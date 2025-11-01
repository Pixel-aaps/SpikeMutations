import pandas as pd
import os
import numpy as np

# --- Paths ---
MUTATIONS_CSV = "data/processed/encoded_sequences.csv"
OUTPUT_CSV    = "data/processed/drug_binding_matrix.csv"

# --- Ligand list (must match names in your combined SDF) ---
ligands = ["apilimod", "camostatmethylate", "cepharanthine", "chloroquine","darunavir", 
           "ebselen", "favipiravir", "hydroxychloroquine", "ivermectin", "lopinavir",
           "molnupiravir", "nelfinavir", "nirmatrelvir", "nitazoxanide", "remdesivir", "ritonavir", 
           "umifenovir"]

# --- Load mutation IDs (only) ---
mut_df = pd.read_csv(MUTATIONS_CSV, usecols=["Mutation_ID"])
mutation_ids = mut_df["Mutation_ID"].tolist()

# --- Prepare empty matrix ---
matrix = pd.DataFrame({"Mutation_ID": mutation_ids})

# --- Fill in placeholder scores (Vina: –12 to –5 kcal/mol) ---
# Replace this block with real score-loading when available
rng = np.random.default_rng(seed=42)
for ligand in ligands:
    matrix[ligand] = rng.uniform(-12.0, -5.0, size=len(matrix))

# --- Save to CSV ---
os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
matrix.to_csv(OUTPUT_CSV, index=False)

print(f"Saved {OUTPUT_CSV} ({len(matrix)} mutations × {len(ligands)} ligands)")
