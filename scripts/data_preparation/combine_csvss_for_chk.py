# scripts/data_preparation/merge_all_features.py
import pandas as pd

# Load datasets
final = pd.read_csv("data/processed/final_dataset_clean.csv")
context = pd.read_csv("data/processed/context_features.csv")

# Filter context to the canonical RBD region
context = context[(context["residue"] >= 331) & (context["residue"] <= 524)]

# Rename for merging
context = context.rename(columns={"residue": "pos", "aa": "pdb_aa"})

# Merge on position
merged = final.merge(context, on="pos", how="left")

print("Before merge:", final.shape)
print("After merge:", merged.shape)

# Save
merged.to_csv("data/processed/final_dataset_with_context.csv", index=False)
print("Saved merged dataset with context features")
