# scripts/004_preprocess_merge.py

import pandas as pd
import numpy as np
import os

# Paths
LABELS_FILE = "data/processed/labels.csv"
EMBEDDINGS_FILE = "data/processed/embeddings.csv"
OUTPUT_FILE = "data/processed/mutation_features.csv"

# Load datasets
labels = pd.read_csv(LABELS_FILE)
emb = pd.read_csv(EMBEDDINGS_FILE)

# Ensure mut_id exists
if "mut_id" not in emb.columns:
    emb = emb.rename(columns={emb.columns[0]: "mut_id"})

# Identify embedding columns
embedding_cols = [c for c in emb.columns if c != "mut_id"]

# 1️⃣ Handle missing embeddings (stop codons or failed embeddings)
rows_allnan = emb[embedding_cols].isna().all(axis=1)
emb.loc[rows_allnan, embedding_cols] = 0.0  # fill missing embeddings with 0
emb["embedding_missing"] = rows_allnan.astype(int)

# 2️⃣ Handle missing DMS / experimental scores
dms_cols = ["bind_lib1", "bind_lib2", "ace2_score", "expr_lib1", "expr_lib2", "expr_score"]

# Fill only NaNs with 0.0 (do NOT overwrite existing 0.0 values)
for col in dms_cols:
    if col in labels.columns:
        labels[col] = labels[col].fillna(0.0)

# Optional: add a flag for DMS score missingness
for col in dms_cols:
    flag_col = col + "_missing"
    labels[flag_col] = labels[col].isna().astype(int)

# 3️⃣ Merge embeddings and labels on mut_id
dataset = pd.merge(labels, emb, on="mut_id", how="inner")

# 4️⃣ Save final ML-ready dataset
dataset.to_csv(OUTPUT_FILE, index=False)
print(f"Final dataset saved to {OUTPUT_FILE}")
print(f"Shape: {dataset.shape}")
print("\nEmbedding missing flag distribution:")
print(dataset["embedding_missing"].value_counts())
