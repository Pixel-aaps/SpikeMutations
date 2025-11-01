# scripts/02_make_labels.py

import os
import pandas as pd
import numpy as np

INTERIM_DIR = "data/interim"
PROCESSED_DIR = "data/processed"
os.makedirs(PROCESSED_DIR, exist_ok=True)

# Load cleaned dataset
df = pd.read_csv(os.path.join(INTERIM_DIR, "dms_core.csv"))

# Thresholds
EXPR_DELETERIOUS = -1.0
ACE2_DELETERIOUS = -1.0
ACE2_BENEFICIAL = 0.5
EXPR_MIN_BENEFICIAL = -0.5

# Classification rules
conditions = [
    (df["expr_score"] < EXPR_DELETERIOUS) | (df["ace2_score"] < ACE2_DELETERIOUS),
    (df["ace2_score"] > ACE2_BENEFICIAL) & (df["expr_score"] > EXPR_MIN_BENEFICIAL),
]
choices = [0, 2]  # 0 = deleterious, 2 = beneficial
df["label"] = np.select(conditions, choices, default=1)  # default = neutral

# # Optional: make it human-readable
# label_map = {0: "deleterious", 1: "neutral", 2: "beneficial"}
# df["label_name"] = df["label"].map(label_map)

# Save label file (keep all columns for now)
out_path = os.path.join(PROCESSED_DIR, "labels.csv")
df.to_csv(out_path, index=False)

# Quick stats
print(f"Saved labeled dataset with {df.shape[0]} mutations → {out_path}")
print(df["label"].value_counts())

print("Deleterious count:", ((df["expr_score"] < EXPR_DELETERIOUS) | (df["ace2_score"] < ACE2_DELETERIOUS)).sum())
print("Beneficial count:", ((df["ace2_score"] > ACE2_BENEFICIAL) & (df["expr_score"] > EXPR_MIN_BENEFICIAL)).sum())

counts = df["label"].value_counts(normalize=True)
print("Class distribution (%):")
print(counts * 100)

print("ACE2 score range:", df["ace2_score"].min(), df["ace2_score"].max())
print("Expr score range:", df["expr_score"].min(), df["expr_score"].max())
