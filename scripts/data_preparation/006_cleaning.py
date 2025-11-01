# preprocess_final_dataset.py
import pandas as pd
import numpy as np
import os

IN_FILE  = "data/processed/final_dataset.csv"
OUT_FILE = "data/processed/final_dataset_clean.csv"

# -----------------------
# Load
# -----------------------
df = pd.read_csv(IN_FILE)
print("Loaded", IN_FILE, "shape:", df.shape)

# -----------------------
# Identify embedding columns
# digits-only column names -> embedding dims produced by ESM
# -----------------------
embed_cols = [c for c in df.columns if c.isdigit()]
print(f"Detected {len(embed_cols)} embedding columns (sample):", embed_cols[:6])

# -----------------------
# Detect stop-codon rows
# We consider a row a stop-codon if:
#   - mutant column contains '*' OR
#   - mut_id contains '*' OR
#   - is_stop == 1 (if present)
# -----------------------
stop_mask = pd.Series(False, index=df.index)

# try common mutant column names
for colname in ("mutant", "mut", "mutation"):
    if colname in df.columns:
        stop_mask = stop_mask | df[colname].astype(str).str.contains(r"\*", regex=True, na=False)

# also check mut_id
if "mut_id" in df.columns:
    stop_mask = stop_mask | df["mut_id"].astype(str).str.contains(r"\*", regex=True, na=False)

# is_stop flag if present
if "is_stop" in df.columns:
    stop_mask = stop_mask | (df["is_stop"] == 1)

n_stop = int(stop_mask.sum())
print("Stop-codon rows detected:", n_stop)

# -----------------------
# Columns to impute ONLY for stop rows
# (only include columns actually present in df)
# -----------------------
exp_cols = [c for c in ["bind_lib1","bind_lib2","ace2_score","expr_lib1","expr_lib2","expr_score"] if c in df.columns]
cols_to_impute = embed_cols + exp_cols
print("Will consider imputation for", len(cols_to_impute), "columns (embeddings + experimental).")

# -----------------------
# Create flags for original NaNs (so you can see what we imputed)
# We create column_name + "_was_nan" for each col we may impute
# -----------------------
for col in cols_to_impute:
    flag = col + "_was_nan"
    # mark 1 where value is NaN (before any fills) and row is a stop row
    df[flag] = 0
    df.loc[stop_mask & df[col].isna(), flag] = 1

# -----------------------
# Fill NaNs with 0.0 ONLY for stop rows
# - do not overwrite non-NaN values
# - this keeps legitimate 0.0 values (like N331N) untouched
# -----------------------
for col in cols_to_impute:
    # Only fill NaNs in the subset of stop rows (safe)
    df.loc[stop_mask, col] = df.loc[stop_mask, col].fillna(0.0)

# -----------------------
# Diagnostics printed for user verification
# -----------------------
print("\nDiagnostics after stop-row imputation:")
# how many embedding dims had NaNs in stop rows originally
embed_flags = [c + "_was_nan" for c in embed_cols if c + "_was_nan" in df.columns]
if embed_flags:
    # rows where any embedding dim was NaN originally
    any_embed_was_nan = (df[embed_flags].sum(axis=1) > 0).sum()
    print("Rows with any embedding-dim originally NaN:", int(any_embed_was_nan))
    all_embed_was_nan = (df[embed_flags].sum(axis=1) == len(embed_flags)).sum()
    print("Rows with ALL embedding-dims originally NaN:", int(all_embed_was_nan))
else:
    print("No embedding flags found (no embedding columns).")

# For experimental columns, show how many stop rows had them NaN originally
for c in exp_cols:
    flag = c + "_was_nan"
    print(f"{c:12s} : was NaN in stop rows = {int(df.loc[stop_mask, flag].sum())}, total NaNs now = {int(df[c].isna().sum())}")

# Show some example rows to confirm (first 10 stop rows)
if n_stop:
    print("\nSample stop rows (first 10) - showing mut_id and the experimental cols and embedding flags:")
    cols_show = ["mut_id"] + exp_cols + (["is_stop"] if "is_stop" in df.columns else []) + embed_flags[:6]
    print(df.loc[stop_mask, cols_show].head(10).to_string(index=False))
else:
    print("No stop rows detected; nothing changed.")

# Confirm N331N (or other known rows) still have 0.0 (don't overwrite legitimate zeros)
for mid in ("N331N","N331E","N331A"):
    if "mut_id" in df.columns and mid in df["mut_id"].values:
        row = df[df["mut_id"]==mid].iloc[0]
        vals = {c: row[c] for c in ["bind_lib1","bind_lib2","ace2_score","expr_lib1","expr_lib2","expr_score"] if c in df.columns}
        print(f"\nSample {mid} values (should be unchanged real values):")
        print(vals)

# -----------------------
# Save cleaned file
# -----------------------
os.makedirs(os.path.dirname(OUT_FILE) or ".", exist_ok=True)
df.to_csv(OUT_FILE, index=False)
print("\nSaved cleaned dataset to", OUT_FILE)
