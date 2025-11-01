import pandas as pd
import numpy as np
import re

# -----------------------------
# Load datasets
# -----------------------------
df_pairs = pd.read_csv("notebooks/ligand_mutation_dataset.csv")
df_wt = pd.read_csv("vina_wt_out/wt_scores.csv")

# -----------------------------
# Normalize LigandName → PubChem CID
# -----------------------------
def extract_cid(name):
    # For dataset (numbers with .0)
    if str(name).replace(".0", "").isdigit():
        return str(int(float(name)))
    # For WT docking names (e.g. Apilimod(10173277)_3D_h)
    m = re.search(r"\((\d+)\)", str(name))
    if m:
        return m.group(1)
    return str(name)

df_pairs["LigandCID"] = df_pairs["LigandName"].apply(extract_cid)
df_wt["LigandCID"] = df_wt["LigandName"].apply(extract_cid)

# -----------------------------
# Merge on CID
# -----------------------------
df = df_pairs.merge(df_wt[["LigandCID","DockingScore_WT"]], on="LigandCID", how="left")

# -----------------------------
# Load WT protein features
# -----------------------------
wt_protein = pd.read_csv("notebooks/protein_features_rich.csv").iloc[0]

# Compute changes
df['Delta_AvgHydro'] = df['AvgHydro'] - wt_protein['AvgHydro']
df['Delta_Charge']   = df['FracCharge'] - wt_protein['FracCharge']

# -----------------------------
# Heuristic correction
# -----------------------------
alpha = -0.7
beta  = 0.5
df['DockingScore_label'] = df['DockingScore_WT'] + alpha*df['Delta_AvgHydro'] + beta*df['Delta_Charge']

# -----------------------------
# Save
# -----------------------------
df.to_csv("ligand_mutation_with_labels_heuristic.csv", index=False)
print("✅ Saved ligand_mutation_with_labels_heuristic.csv with", len(df), "rows")
