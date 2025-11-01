#!/usr/bin/env python3
"""
07_sanity_check_prepared_datasets.py

Sanity-check evaluation on:
  - mutation_features.csv  (labels + embeddings)
  - final_dataset_clean.csv (labels + embeddings + physicochemical)

Checks:
  - GroupKFold CV grouped by site_SARS2
  - Macro-F1, balanced accuracy, accuracy
  - Bootstrap confidence intervals for macro-F1 and difference
"""

import os
import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score, balanced_accuracy_score, accuracy_score, classification_report

# ---------- Helpers ----------
def drop_label_proxies(df):
    # remove columns that leak label info
    banned = ["ace2", "expr", "bind", "score"]
    drop_cols = [c for c in df.columns if any(b in c.lower() for b in banned)]
    drop_cols = [c for c in drop_cols if c not in ["mut_id","site_SARS2","label"]]
    if drop_cols:
        print("Dropping proxy cols:", drop_cols)
        df = df.drop(columns=drop_cols)
    return df

def get_oof(X, y, groups, n_splits=5, random_state=42):
    gkf = GroupKFold(n_splits=n_splits)
    preds = np.zeros(len(y), dtype=int)
    for fold,(tr,te) in enumerate(gkf.split(X,y,groups),1):
        model = RandomForestClassifier(
            n_estimators=300, random_state=random_state,
            class_weight="balanced", n_jobs=-1
        )
        model.fit(X.iloc[tr], y[tr])
        preds[te] = model.predict(X.iloc[te])
        print(f"  Fold {fold}: train {len(tr)} test {len(te)}")
    return preds

def bootstrap_metric(y, yhat, metric, n_boot=1000, random_state=42):
    rng = np.random.RandomState(random_state)
    scores = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.randint(0,n,n)
        scores.append(metric(y[idx], yhat[idx]))
    return np.array(scores)

# ---------- Main ----------
def main(args):
    # Paths
    f_mut = "data/processed/mutation_features.csv"
    f_final = "data/processed/final_dataset_clean.csv"

    df_mut = pd.read_csv(f_mut)
    df_final = pd.read_csv(f_final)

    # Ensure required columns
    for df, name in [(df_mut,"mutation_features"), (df_final,"final_dataset_clean")]:
        for col in ["mut_id","site_SARS2","label"]:
            if col not in df.columns:
                raise RuntimeError(f"{name} missing column {col}")

    # Drop proxies
    df_mut = drop_label_proxies(df_mut)
    df_final = drop_label_proxies(df_final)

    # Prepare X,y,groups
    def prepare(df):
        y = df["label"].astype(int).values
        groups = df["site_SARS2"].values
        X = df.drop(columns=["mut_id","site_SARS2","label"])
        # One-hot encode if categorical
        X = pd.get_dummies(X, drop_first=False)
        X = X.fillna(0)
        return X, y, groups

    X_mut, y_mut, g_mut = prepare(df_mut)
    X_final, y_final, g_final = prepare(df_final)

    # Restrict to common rows (same mut_ids) for fair comparison
    common_ids = set(df_mut["mut_id"]).intersection(df_final["mut_id"])
    df_mut = df_mut[df_mut["mut_id"].isin(common_ids)].reset_index(drop=True)
    df_final = df_final[df_final["mut_id"].isin(common_ids)].reset_index(drop=True)
    X_mut, y_mut, g_mut = prepare(df_mut)
    X_final, y_final, g_final = prepare(df_final)

    print("N samples (common):", len(df_mut))

    # Run OOF
    print("\nEmbeddings-only model + labels(ace2Score + exprScore) (mutation_features.csv)")
    pred_mut = get_oof(X_mut, y_mut, g_mut, n_splits=args.n_splits, random_state=args.seed)
    print("Macro-F1:", f1_score(y_mut, pred_mut, average="macro"))
    print(classification_report(y_mut, pred_mut))

    print("\nEmbeddings+physchem model + labels(ace2Score + exprScore) (final_dataset_clean.csv)")
    pred_final = get_oof(X_final, y_final, g_final, n_splits=args.n_splits, random_state=args.seed)
    print("Macro-F1:", f1_score(y_final, pred_final, average="macro"))
    print(classification_report(y_final, pred_final))

    # Bootstrap CIs
    print("\nBootstrap CIs (n=%d)" % args.n_boot)
    m1 = bootstrap_metric(y_mut, pred_mut, lambda a,b:f1_score(a,b,average="macro"),
                          n_boot=args.n_boot, random_state=args.seed)
    m2 = bootstrap_metric(y_final, pred_final, lambda a,b:f1_score(a,b,average="macro"),
                          n_boot=args.n_boot, random_state=args.seed+1)
    diff = m2 - m1
    print("Embeddings-only F1: mean=%.3f, 95%% CI=(%.3f,%.3f)" %
          (m1.mean(), np.percentile(m1,2.5), np.percentile(m1,97.5)))
    print("Emb+physchem F1:   mean=%.3f, 95%% CI=(%.3f,%.3f)" %
          (m2.mean(), np.percentile(m2,2.5), np.percentile(m2,97.5)))
    print("Difference (final - mut): mean=%.3f, 95%% CI=(%.3f,%.3f), P>0=%.3f" %
          (diff.mean(), np.percentile(diff,2.5), np.percentile(diff,97.5),
           np.mean(diff>0)))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_splits", type=int, default=5)
    ap.add_argument("--n_boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    main(args)
    print("-"*50)

