import pandas as pd
import os

# ==== Load files ====
labels = pd.read_csv("data/processed/labels.csv")
seq_features = pd.read_csv("data/processed/seq_features.csv")
embeddings = pd.read_csv("data/processed/embeddings.csv") 

# ==== Merge ====
df = labels.merge(seq_features, on="mut_id", how="inner")
df = df.merge(embeddings, on="mut_id", how="inner")

print("Final dataset shape:", df.shape)

# ==== Save ====
os.makedirs("data/processed", exist_ok=True)
df.to_csv("data/processed/final_dataset.csv", index=False)

print("Saved data/processed/final_dataset.csv")
