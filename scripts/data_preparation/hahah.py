#!/usr/bin/env python3
"""
09_region_holdout_check.py

Evaluate model generalization with region-based holdout:
 - Train on all residues outside the region
 - Test on residues inside the region
 - Reports macro-F1, accuracy
"""

import os
import argparse
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score, accuracy_score, classification_report

# ----------------------------
def drop_label_proxies(df):
    banned = ["ace2", "expr", "bind", "score"]
    drop_cols = [c for c in df.columns if any(b in c.lower() for b in banned)]
    drop_cols = [c for c in drop_cols if c not in ["mut_id","site_SARS2","label"]]
    if drop_cols:
        print("Dropping proxy cols:", drop_cols)
        df = df.drop(columns=drop_cols)
    return df

def prepare(df):
    y = df["label"].astype(int).values
    X = df.drop(columns=["mut_id","site_SARS2","label"])
    X = pd.get_dummies(X, drop_first=False).fillna(0)
    return X, y

# ----------------------------
def region_holdout(df, region=(438,507), seed=42):
    """
    region = tuple(start, end) inclusive, based on site_SARS2 numbering
    """
    # Split by site_SARS2
    mask_region = df["site_SARS2"].between(region[0], region[1])

    df_train = df[~mask_region].reset_index(drop=True)
    df_test  = df[mask_region].reset_index(drop=True)

    X_train, y_train = prepare(df_train)
    X_test, y_test   = prepare(df_test)

    print(f"Training on {len(df_train)} mutations (outside {region[0]}–{region[1]})")
    print(f"Testing on {len(df_test)} mutations (inside {region[0]}–{region[1]})")

    model = RandomForestClassifier(
        n_estimators=300, random_state=seed,
        class_weight="balanced", n_jobs=-1
    )
    model.fit(X_train, y_train)
    yhat = model.predict(X_test)

    f1 = f1_score(y_test, yhat, average="macro")
    acc = accuracy_score(y_test, yhat)

    print("\n=== Region-holdout results ===")
    print(f"Macro-F1: {f1:.3f}")
    print(f"Accuracy: {acc:.3f}")
    print("\nClassification report (holdout region):")
    print(classification_report(y_test, yhat))

    return f1, acc

# ----------------------------
def main(args):
    df = pd.read_csv(args.dataset)
    print("Loaded", args.dataset, "shape:", df.shape)

    df = drop_label_proxies(df)

    region = (args.start, args.end)
    region_holdout(df, region=region, seed=args.seed)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", type=str,
                    default="data/processed/final_dataset_clean.csv",
                    help="Path to dataset CSV")
    ap.add_argument("--start", type=int, default=438,
                    help="Start of region (inclusive)")
    ap.add_argument("--end", type=int, default=507,
                    help="End of region (inclusive)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    main(args)
