#!/usr/bin/env python3
"""
scripts/05_ablation_train.py

Ablation harness for Spike mutation classifier.
- Uses GroupKFold (group by site_SARS2) to avoid leakage across positions.
- Runs several feature-set ablations (embeddings only, +physicochem, +context, +stability if present).
- Computes macro-F1, balanced accuracy, ROC-AUC (OVR), and ECE (expected calibration error).
- Computes permutation importance for the final model per ablation.
- Optionally computes SHAP summary (if shap is installed).

Usage:
    python scripts/05_ablation_train.py --n_splits 5 --random_state 42 --model rf
"""
import os
import argparse
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import (f1_score, balanced_accuracy_score,
                             classification_report, roc_auc_score)
from sklearn.preprocessing import LabelBinarizer
from sklearn.inspection import permutation_importance
from joblib import dump
import warnings
warnings.filterwarnings("ignore")

# ---------- Utilities ----------
def compute_ece(probs, y_true, n_bins=10):
    """
    Simple multiclass ECE approximation:
    - For each sample use the predicted max probability p_max and correctness (1 if predicted == true).
    - Bin by p_max and compute weighted |acc_bin - avg_conf_bin|.
    """
    preds = np.argmax(probs, axis=1)
    p_max = probs.max(axis=1)
    correct = (preds == y_true).astype(int)
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for i in range(n_bins):
        idx = (p_max >= bins[i]) & (p_max < bins[i+1])
        if np.sum(idx) == 0:
            continue
        acc_bin = correct[idx].mean()
        conf_bin = p_max[idx].mean()
        ece += (np.sum(idx) / len(y_true)) * abs(acc_bin - conf_bin)
    return ece

def safe_merge(left, right, on_left, on_right=None):
    if on_right is None:
        on_right = on_left
    return left.merge(right, left_on=on_left, right_on=on_right, how='left')

# ---------- Data loading & merges ----------
def load_and_merge_data(paths):
    """
    Expected files (if present):
      - data/processed/labels.csv      (required)  : has mut_id, site_SARS2, label
      - data/processed/embeddings.csv  (optional)  : mut_id + embedding cols
      - data/processed/mutation_features_enhanced.csv (optional) : physchem cols
      - data/processed/context_features.csv (optional) : rsa, dist_to_ace2, is_contact, is_rbm, sec_struct
      - data/processed/stability_features.csv (optional): ddg or ddg_pred column
      - data/processed/conservation.csv (optional): pos-wise conservation
    """
    labels = pd.read_csv(paths['labels'])
    labels = labels.rename(columns={c: c.strip() for c in labels.columns})
    if 'mut_id' not in labels.columns:
        raise ValueError("labels.csv must contain 'mut_id' column")
    # base frame
    df = labels[['mut_id', 'site_SARS2', 'label']].copy()
    # embeddings
    if os.path.exists(paths['embeddings']):
        emb = pd.read_csv(paths['embeddings'])
        # make sure mut_id column exists in embeddings
        if 'mut_id' not in emb.columns and 'Mutation_ID' in emb.columns:
            emb = emb.rename(columns={'Mutation_ID': 'mut_id'})
        if 'mut_id' not in emb.columns:
            raise ValueError("embeddings.csv must contain 'mut_id' or 'Mutation_ID'")
        df = safe_merge(df, emb, 'mut_id')
    # physchem (mutation_features_enhanced)
    if os.path.exists(paths['physchem']):
        phys = pd.read_csv(paths['physchem'])
        # try common id names
        idcol = None
        for cand in ['mut_id', 'Mutation_ID', 'MutationId', 'MutationId']:
            if cand in phys.columns:
                idcol = cand
                break
        if idcol is None and 'pos' in phys.columns and 'wt' in phys.columns and 'mut' in phys.columns:
            # fallback: build mut_id from wt+pos+mut if possible
            phys['mut_id'] = phys['wt'].astype(str) + phys['pos'].astype(int).astype(str) + phys['mut'].astype(str)
            idcol = 'mut_id'
        if idcol is None:
            raise ValueError("Could not find id column in physchem file")
        df = safe_merge(df, phys, 'mut_id', idcol)
    # context
    if os.path.exists(paths['context']):
        ctx = pd.read_csv(paths['context'])
        # context should have mut_id and / or pos
        if 'mut_id' in ctx.columns:
            df = safe_merge(df, ctx, 'mut_id')
        elif 'pos' in ctx.columns:
            df = df.merge(ctx, left_on='site_SARS2', right_on='pos', how='left')
        else:
            raise ValueError("context_features.csv must have 'mut_id' or 'pos'")
    # stability
    if os.path.exists(paths['stability']):
        stab = pd.read_csv(paths['stability'])
        idcol = 'mut_id' if 'mut_id' in stab.columns else ('Mutation_ID' if 'Mutation_ID' in stab.columns else None)
        if idcol:
            df = safe_merge(df, stab, 'mut_id', idcol)
    # conservation (optional)
    if os.path.exists(paths['conservation']):
        cons = pd.read_csv(paths['conservation'])
        # expect pos, conservation_score
        if 'pos' in cons.columns:
            df = df.merge(cons, left_on='site_SARS2', right_on='pos', how='left')
    # final clean
    df = df.drop_duplicates(subset=['mut_id'])
    return df

# ---------- Feature set builders ----------
def get_feature_columns_for_ablation(df):
    # detect embeddings: everything numeric except site_SARS2,label and some id fields
    non_feat = set(['mut_id','site_SARS2','label','pos','Mutation_ID','wt','mut','wildtype','mutant','mutation','mutation_RBD'])
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    # remove non_feat
    numeric_cols = [c for c in numeric_cols if c not in non_feat]
    # embed columns heuristic: long vector columns (many dims)
    emb_cols = []
    # If embeddings file was present and contains typical names, prefer those by checking count or prefix
    for c in df.columns:
        if c.startswith('esm') or c.startswith('embedding') or c.startswith('emb_') or c.startswith('esm2') or c.startswith('prot_'):
            emb_cols.append(c)
    # fallback: if there are many numeric cols (>50) use those as embeddings (assuming emb dims)
    if len(emb_cols) == 0 and len(numeric_cols) > 50:
        emb_cols = numeric_cols.copy()
    # physchem candidates (explicit known names)
    physchem_cands = ['delta_hydrophobicity','delta_charge','volume_delta','polarity_changed','blosum62',
                      'is_proline_introduced','is_gly_to_x','aromatic_change','is_conservative']
    physchem_cols = [c for c in physchem_cands if c in df.columns]
    # context candidates
    context_cands = ['rsa','dist_to_ace2','is_contact','is_rbm','contact_number','sec_struct']
    ctx_cols = [c for c in context_cands if c in df.columns]
    # stability
    stab_cols = [c for c in ['ddg','ddG','ddg_pred','ddG_pred','ddG_foldx'] if c in df.columns]
    # conservation
    cons_cols = [c for c in ['conservation_score','entropy','shannon'] if c in df.columns]
    return {
        'embeddings': emb_cols,
        'physchem': physchem_cols,
        'context': ctx_cols,
        'stability': stab_cols,
        'conservation': cons_cols
    }

# ---------- Cross-validated training ----------
def run_group_cv(df, feature_sets, n_splits=5, random_state=42, model_type='rf', do_shap=False, outdir='results'):
    os.makedirs(outdir, exist_ok=True)
    groups = df['site_SARS2'].values
    y = df['label'].astype(int).values
    gkf = GroupKFold(n_splits=n_splits)
    results = []
    fold_idx = 0
    for ab_name, use_cols in feature_sets.items():
        print(f"\n=== Ablation: {ab_name} ({len(use_cols)} features) ===")
        # check features exist
        if len(use_cols) == 0:
            print(" -> No features found for this ablation, skipping.")
            continue
        # data matrix
        X = df[use_cols].copy()
        # handle categorical sec_struct if present
        if 'sec_struct' in X.columns:
            X = pd.get_dummies(X, columns=['sec_struct'], dummy_na=True)
        # fill missing
        X = X.fillna(0)
        all_index = np.arange(len(df))
        fold_metrics = []
        fold_preds_proba = np.zeros((len(df), len(np.unique(y))))
        # iterate folds
        for train_idx, test_idx in gkf.split(all_index, y, groups=groups):
            fold_idx += 1
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            # model
            if model_type == 'lgb':
                try:
                    from lightgbm import LGBMClassifier
                    model = LGBMClassifier(random_state=random_state, n_estimators=300)
                except Exception:
                    print("lightgbm not installed, falling back to RandomForest")
                    model = RandomForestClassifier(n_estimators=200, random_state=random_state, n_jobs=-1, class_weight='balanced')
            else:
                model = RandomForestClassifier(n_estimators=200, random_state=random_state, n_jobs=-1, class_weight='balanced')
            model.fit(X_train, y_train)
            # predict
            proba = model.predict_proba(X_test)
            preds = np.argmax(proba, axis=1)
            # store predicted probs for whole dataset
            fold_preds_proba[test_idx] = proba
            # metrics
            f1_macro = f1_score(y_test, preds, average='macro')
            bal_acc = balanced_accuracy_score(y_test, preds)
            # ROC AUC (ovr)
            lb = LabelBinarizer()
            lb.fit(y)
            y_test_bin = lb.transform(y_test)
            try:
                roc_auc = roc_auc_score(y_test_bin, proba, average='macro', multi_class='ovr')
            except Exception:
                roc_auc = np.nan
            ece = compute_ece(proba, y_test, n_bins=10)
            fold_metrics.append({
                'ablation': ab_name,
                'fold': len(fold_metrics) + 1,
                'f1_macro': f1_macro,
                'balanced_accuracy': bal_acc,
                'roc_auc_ovr': roc_auc,
                'ece': ece,
                'n_test': len(test_idx)
            })
            # permutation importance on this fold's test set (store top 20)
            try:
                perm = permutation_importance(model, X_test, y_test, n_repeats=10, random_state=random_state, scoring='f1_macro', n_jobs=-1)
                perm_df = pd.DataFrame({'feature': X.columns, 'imp_mean': perm.importances_mean, 'imp_std': perm.importances_std})
                perm_df = perm_df.sort_values('imp_mean', ascending=False).head(40)
                perm_df.to_csv(os.path.join(outdir, f'perm_importance_{ab_name}_fold{len(fold_metrics)}.csv'), index=False)
            except Exception as e:
                print("Permutation importance failed:", e)
            # optionally save model of last fold
            dump(model, os.path.join(outdir, f'model_{ab_name}_fold{len(fold_metrics)}.joblib'))
            # optional SHAP (compute for last fold only to save time)
            if do_shap:
                try:
                    import shap
                    explainer = shap.TreeExplainer(model)
                    shap_vals = explainer.shap_values(X_test)
                    # produce a small summary plot (requires matplotlib)
                    shap.summary_plot(shap_vals, X_test, show=False)
                    import matplotlib.pyplot as plt
                    plt.tight_layout()
                    plt.savefig(os.path.join(outdir, f'shap_summary_{ab_name}_fold{len(fold_metrics)}.png'))
                    plt.close()
                except Exception as e:
                    print("SHAP failed:", e)
        # aggregate fold metrics
        fold_df = pd.DataFrame(fold_metrics)
        fold_df.to_csv(os.path.join(outdir, f'ablation_{ab_name}_fold_metrics.csv'), index=False)
        # aggregate summary
        summary = fold_df.agg(['mean','std'])[['f1_macro','balanced_accuracy','roc_auc_ovr','ece']].T
        summary.columns = ['mean','std']
        summary['ablation'] = ab_name
        # write to global results list
        for r in fold_metrics:
            results.append(r)
        # Save predicted probs for analysis (optional)
        # Save a CSV with sample-level probs if desired:
        proba_df = pd.DataFrame(fold_preds_proba, columns=[f"class_{c}" for c in sorted(np.unique(y))])
        proba_df['mut_id'] = df['mut_id'].values
        proba_df['y_true'] = y
        proba_df.to_csv(os.path.join(outdir, f'preds_proba_{ab_name}.csv'), index=False)
        print(f"Ablation {ab_name} mean F1_macro: {summary.loc['f1_macro','mean']:.4f} ± {summary.loc['f1_macro','std']:.4f}")
    # save full fold-level results
    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(outdir, 'ablation_fold_results.csv'), index=False)
    print("Saved results to", outdir)
    return results_df

# ---------- Main ----------
if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--n_splits', type=int, default=5)
    p.add_argument('--random_state', type=int, default=42)
    p.add_argument('--model', choices=['rf','lgb'], default='rf')
    p.add_argument('--do_shap', action='store_true')
    args = p.parse_args()

    paths = {
        'labels': 'data/processed/labels.csv',
        'embeddings': 'data/processed/embeddings.csv',
        'physchem': 'data/processed/mutation_features_enhanced.csv',
        'context': 'data/processed/context_features.csv',
        'stability': 'data/processed/stability_features.csv',
        'conservation': 'data/processed/conservation.csv'
    }

    print("Loading and merging data...")
    df = load_and_merge_data(paths)
    print("Rows:", len(df), "Columns:", df.shape[1])

    print("Detecting available feature groups...")
    available = get_feature_columns_for_ablation(df)
    print("Detected features:", {k: len(v) for k,v in available.items()})

    # Build ablation sets in order (change order if you prefer)
    feature_sets = {}
    if len(available['embeddings'])>0:
        feature_sets['embeddings'] = available['embeddings']
    else:
        # if no embeddings, fallback to numeric features (but warn)
        numeric = df.select_dtypes(include=[np.number]).columns.tolist()
        numeric = [c for c in numeric if c not in ['site_SARS2','label']]
        feature_sets['numeric_all'] = numeric

    # embeddings + physchem
    feature_sets['embeddings_plus_physchem'] = feature_sets[list(feature_sets.keys())[0]] + available['physchem']
    # add context
    feature_sets['embeddings_plus_physchem_plus_context'] = feature_sets['embeddings_plus_physchem'] + [c for c in available['context'] if c!='sec_struct']
    # if sec_struct exists add later as one-hot in the pipeline automatically
    # add stability
    feature_sets['all_including_stability'] = feature_sets['embeddings_plus_physchem_plus_context'] + available['stability'] + available['conservation']

    # run CV
    results_df = run_group_cv(df, feature_sets, n_splits=args.n_splits, random_state=args.random_state, model_type=args.model, do_shap=args.do_shap)
    # print summary
    print("\n=== Aggregate summary ===")
    print(results_df.groupby('ablation')[['f1_macro','balanced_accuracy','roc_auc_ovr','ece']].mean())
