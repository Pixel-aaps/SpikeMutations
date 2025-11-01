#!/usr/bin/env python3
"""
06_sanity_check_evaluation.py

Sanity-check evaluation for spike mutation classifier:
 - GroupKFold by site_SARS2 to avoid leakage across residues
 - Compare embeddings-only vs embeddings+physicochemical features
 - Compute macro-F1, balanced accuracy, ROC-AUC (OVR)
 - Bootstrap 95% CIs for macro-F1 and the difference

Usage:
    python scripts/06_sanity_check_evaluation.py --n_splits 5 --n_boot 1000 --random_state 42
"""

import os
import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score, balanced_accuracy_score, accuracy_score, roc_auc_score, classification_report
from sklearn.preprocessing import LabelBinarizer

# ---------- Helpers ----------
def find_id_col(df, candidates=("mut_id","Mutation_ID","MutationId","id")):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def drop_label_proxies(df):
    # remove columns that look like direct measurement proxies for labels
    banned_patterns = ['ace2', 'expr', 'bind', 'ace2_score', 'expr_score', 'ace2score', 'exprscore']
    cols_to_drop = [c for c in df.columns if any(p in c.lower() for p in banned_patterns)]
    # but do not drop site or label or mut_id
    cols_to_drop = [c for c in cols_to_drop if c not in ('mut_id','site_SARS2','label','mut_id')]
    if cols_to_drop:
        print("Dropping label-proxy columns:", cols_to_drop)
        df = df.drop(columns=cols_to_drop)
    return df

def get_bootstrap_ci(vec, alpha=0.05):
    lo = np.percentile(vec, 100*alpha/2)
    hi = np.percentile(vec, 100*(1-alpha/2))
    return lo, hi

# ---------- Core CV to get OOF preds ----------
def get_oof_predictions(X, y, groups, n_splits=5, random_state=42, model_kwargs=None):
    if model_kwargs is None:
        model_kwargs = {"n_estimators":300, "random_state":random_state, "n_jobs":-1, "class_weight":"balanced"}
    classes = np.unique(y)
    n_classes = len(classes)
    n_samples = X.shape[0]
    oof_preds = np.zeros(n_samples, dtype=int)
    oof_proba = np.zeros((n_samples, n_classes), dtype=float)
    oof_models = []
    gkf = GroupKFold(n_splits=n_splits)
    for fold, (train_idx, test_idx) in enumerate(gkf.split(X, y, groups=groups), start=1):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        model = RandomForestClassifier(**model_kwargs)
        model.fit(X_train, y_train)
        # probabilities: align with global 'classes'
        proba = model.predict_proba(X_test)  # shape (n_test, n_classes_in_model)
        # map model.classes_ to indices in global classes
        for i_cls, cls in enumerate(model.classes_):
            global_idx = int(np.where(classes == cls)[0][0])
            oof_proba[test_idx, global_idx] = proba[:, i_cls]
        preds = model.predict(X_test)
        oof_preds[test_idx] = preds
        oof_models.append(model)
        print(f"  Fold {fold}: trained on {len(train_idx)} samples, test {len(test_idx)} samples")
    return oof_preds, oof_proba, oof_models

# ---------- Bootstrap routine ----------
def bootstrap_metric(y_true, y_pred, metric_func, n_boot=1000, random_state=42):
    rng = np.random.RandomState(random_state)
    n = len(y_true)
    stats = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.randint(0, n, size=n)  # sample with replacement
        stats[i] = metric_func(y_true[idx], y_pred[idx])
    return stats

# ---------- Main ----------
def main(args):
    # Paths (adjust if your layout differs)
    labels_path = "data/processed/labels.csv"
    embeddings_path = "data/processed/embeddings.csv"  # produced from ESM2
    physchem_path = "data/processed/seq_features.csv"  # sequence physchem features

    # Load labels
    labels = pd.read_csv(labels_path)
    # Expect 'mut_id', 'site_SARS2', 'label' exist
    if 'mut_id' not in labels.columns:
        raise RuntimeError("labels.csv must contain 'mut_id' column")
    if 'site_SARS2' not in labels.columns:
        raise RuntimeError("labels.csv must contain 'site_SARS2' column for grouping")
    if 'label' not in labels.columns:
        raise RuntimeError("labels.csv must contain 'label' column")

    # Base frame
    df = labels[['mut_id','site_SARS2','label']].copy()

    # Merge embeddings
    if not os.path.exists(embeddings_path):
        raise RuntimeError(f"Embeddings file missing: {embeddings_path}")
    emb = pd.read_csv(embeddings_path)
    emb_id = find_id_col(emb)
    if emb_id is None:
        raise RuntimeError("embeddings.csv must contain an id column (mut_id or Mutation_ID)")
    if emb_id != 'mut_id':
        emb = emb.rename(columns={emb_id: 'mut_id'})
    df = df.merge(emb, on='mut_id', how='left')
    print("Merged embeddings; shape:", df.shape)

    # Merge physchem if present
    physchem_exists = os.path.exists(physchem_path)
    if physchem_exists:
        phys = pd.read_csv(physchem_path)
        # find id col
        phys_id = find_id_col(phys)
        if phys_id is None:
            # try to build mut_id from wt,pos,mut
            if set(['wt','pos','mut']).issubset(set([c.lower() for c in phys.columns])):
                phys['mut_id'] = phys['wt'].astype(str) + phys['pos'].astype(int).astype(str) + phys['mut'].astype(str)
                phys_id = 'mut_id'
            else:
                raise RuntimeError("Could not find id column in physchem file")
        if phys_id != 'mut_id':
            phys = phys.rename(columns={phys_id: 'mut_id'})
        df = df.merge(phys, on='mut_id', how='left', suffixes=('', '_phys'))
        print("Merged physchem; shape:", df.shape)
    else:
        print("Physicochemical file not found; proceeding with embeddings-only comparison.")

    # Drop label proxies
    df = drop_label_proxies(df)

    # Determine embedding columns: everything from embeddings file (except id)
    emb_cols = [c for c in emb.columns if c != 'mut_id']
    # Determine physchem cols: numeric columns from physchem (excluding id and obvious meta)
    phys_cols = []
    if physchem_exists:
        # Keep numeric cols that are not in the non-feature set
        non_feat = set(['mut_id','pos','Mutation_ID','wt','mut','label'])
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        # but numeric_cols includes embeddings too; remove emb_cols
        phys_candidate = [c for c in numeric_cols if c not in emb_cols and c not in non_feat and c not in ('site_SARS2',)]
        # optionally restrict to a known list if present
        known_phys = ['delta_hydrophobicity','delta_charge','volume_delta','polarity_changed','blosum62',
                      'is_proline_introduced','is_gly_to_x','aromatic_change','is_conservative']
        phys_cols = [c for c in known_phys if c in df.columns]
        # if none of the known_phys are found, use the numeric phys_candidate
        if len(phys_cols) == 0:
            phys_cols = phys_candidate
    print("Embedding dims:", len(emb_cols), "; Physchem dims (selected):", len(phys_cols))

    # Build two datasets
    # 1) embeddings-only
    df_embed = df[['mut_id','site_SARS2','label'] + emb_cols].copy()
    # 2) embeddings + physchem
    if phys_cols:
        df_embed_phys = df[['mut_id','site_SARS2','label'] + emb_cols + phys_cols].copy()
    else:
        df_embed_phys = df_embed.copy()

    # Drop rows with missing embeddings
    df_embed = df_embed.dropna(subset=emb_cols)
    df_embed_phys = df_embed_phys.dropna(subset=emb_cols)  # ensure embeddings present

    # Prepare X, y, groups
    def prepare(df_local, use_cols):
        X = df_local[use_cols].copy()
        # if sec_struct exists as string, one-hot encode
        if 'sec_struct' in X.columns:
            X = pd.get_dummies(X, columns=['sec_struct'], dummy_na=True)
        # fill missing with 0 (alternative: median)
        X = X.fillna(0)
        y = df_local['label'].astype(int).values
        groups = df_local['site_SARS2'].values
        return X, y, groups

    X_embed, y_embed, g_embed = prepare(df_embed, emb_cols)
    X_embed_phys, y_embed_phys, g_embed_phys = prepare(df_embed_phys, emb_cols + phys_cols if phys_cols else emb_cols)

    if not (np.array_equal(y_embed, y_embed_phys) and np.array_equal(g_embed, g_embed_phys) and df_embed['mut_id'].tolist() == df_embed_phys['mut_id'].tolist()):
        # align rows: restrict to intersection of mut_id to be fair
        common_ids = set(df_embed['mut_id']).intersection(set(df_embed_phys['mut_id']))
        df_embed = df_embed[df_embed['mut_id'].isin(common_ids)].reset_index(drop=True)
        df_embed_phys = df_embed_phys[df_embed_phys['mut_id'].isin(common_ids)].reset_index(drop=True)
        X_embed, y_embed, g_embed = prepare(df_embed, emb_cols)
        X_embed_phys, y_embed_phys, g_embed_phys = prepare(df_embed_phys, emb_cols + phys_cols if phys_cols else emb_cols)
        print("Aligned both datasets to common mut_id set; n =", len(df_embed))

    # Run GroupKFold OOF predictions for both models
    print("\nRunning GroupKFold CV for embeddings-only model...")
    oof_preds_emb, oof_proba_emb, models_emb = get_oof_predictions(X_embed, y_embed, g_embed, n_splits=args.n_splits, random_state=args.random_state)

    print("\nRunning GroupKFold CV for embeddings + physchem model...")
    oof_preds_ep, oof_proba_ep, models_ep = get_oof_predictions(X_embed_phys, y_embed_phys, g_embed_phys, n_splits=args.n_splits, random_state=args.random_state)

    # Evaluate metrics
    def evaluate(y_true, y_pred, y_proba):
        f1_macro = f1_score(y_true, y_pred, average='macro')
        bal_acc = balanced_accuracy_score(y_true, y_pred)
        acc = accuracy_score(y_true, y_pred)
        lb = LabelBinarizer().fit(y_true)
        try:
            roc = roc_auc_score(lb.transform(y_true), y_proba, average='macro', multi_class='ovr')
        except Exception:
            roc = np.nan
        return {"f1_macro": f1_macro, "balanced_accuracy": bal_acc, "accuracy": acc, "roc_auc_ovr": roc}

    eval_emb = evaluate(y_embed, oof_preds_emb, oof_proba_emb)
    eval_ep  = evaluate(y_embed_phys, oof_preds_ep, oof_proba_ep)

    print("\n--- Embeddings-only metrics (OOF) ---")
    for k,v in eval_emb.items():
        print(f"{k}: {v:.4f}")
    print("\nClassification report (embeddings-only):")
    print(classification_report(y_embed, oof_preds_emb))

    print("\n--- Embeddings + physchem metrics (OOF) ---")
    for k,v in eval_ep.items():
        print(f"{k}: {v:.4f}")
    print("\nClassification report (emb+phys):")
    print(classification_report(y_embed_phys, oof_preds_ep))

    # Bootstrap CIs for macro-F1
    n_boot = args.n_boot
    print(f"\nPerforming bootstrap (n={n_boot}) for macro-F1...")
    stats_emb = bootstrap_metric(y_embed, oof_preds_emb, lambda a,b: f1_score(a,b,average='macro'), n_boot=n_boot, random_state=args.random_state)
    stats_ep  = bootstrap_metric(y_embed_phys, oof_preds_ep, lambda a,b: f1_score(a,b,average='macro'), n_boot=n_boot, random_state=args.random_state+1)

    lo_e, hi_e = get_bootstrap_ci(stats_emb)
    lo_ep, hi_ep = get_bootstrap_ci(stats_ep)
    print(f"Embeddings-only macro-F1: mean={stats_emb.mean():.4f}, 95% CI=({lo_e:.4f}, {hi_e:.4f})")
    print(f"Embeddings+physchem macro-F1: mean={stats_ep.mean():.4f}, 95% CI=({lo_ep:.4f}, {hi_ep:.4f})")

    # Bootstrap difference distribution (ep - emb)
    diff = stats_ep - stats_emb
    lo_d, hi_d = get_bootstrap_ci(diff)
    p_greater = np.mean(diff > 0)  # proportion of bootstraps where ep > emb
    print(f"Difference (emb+phys - emb): mean={diff.mean():.4f}, 95% CI=({lo_d:.4f}, {hi_d:.4f}), P(ep>emb)={p_greater:.3f}")

    # Save OOF predictions for later analysis
    outdir = "results/sanity_check"
    os.makedirs(outdir, exist_ok=True)
    save_df = pd.DataFrame({
        'mut_id': df_embed['mut_id'],
        'site_SARS2': df_embed['site_SARS2'],
        'y_true': y_embed,
        'pred_emb': oof_preds_emb,
        'pred_ep': oof_preds_ep
    })
    # add probs columns
    for i in range(oof_proba_emb.shape[1]):
        save_df[f"prob_emb_class{i}"] = oof_proba_emb[:, i]
    for i in range(oof_proba_ep.shape[1]):
        save_df[f"prob_ep_class{i}"] = oof_proba_ep[:, i]
    save_df.to_csv(os.path.join(outdir, "oof_predictions.csv"), index=False)
    print("Saved OOF predictions to", os.path.join(outdir, "oof_predictions.csv"))

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--n_splits', type=int, default=5, help='GroupKFold splits')
    p.add_argument('--n_boot', type=int, default=1000, help='Bootstrap iterations for CI')
    p.add_argument('--random_state', type=int, default=42)
    args = p.parse_args()
    main(args)
