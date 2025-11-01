import pandas as pd
import numpy as np

# thresholds (tunable)
EXPR_DELETERIOUS = -1.0
ACE2_DELETERIOUS = -1.0
ACE2_BENEFICIAL = 0.5
EXPR_MIN_BENEFICIAL = -0.5
RSA_BURIED = 0.2        # relative solvent accessibility threshold
DDG_DESTABILIZING = 1.0 # kcal/mol
CONSERVED = 0.8         # high conservation score

# known Variants of Concern (override as beneficial)
VOC_MUTS = {"N501Y", "E484K", "L452R", "K417N", "K417T"}
labels = pd.read_csv("data/labels.csv")
seq_features = pd.read_csv("data/seq_features.csv")
context = pd.read_csv("data/processed/context_features.csv")

# optional features (if available later)
try:
    conservation = pd.read_csv("data/processed/conservation.csv")
except FileNotFoundError:
    conservation = pd.DataFrame(columns=["mut_id", "conservation"])

try:
    stability = pd.read_csv("data/processed/stability_features.csv")
except FileNotFoundError:
    stability = pd.DataFrame(columns=["mut_id", "ddG"])

# merge on mut_id
df = labels.merge(seq_features, on="mut_id", how="left")
df = df.merge(context, left_on="site_RBD", right_on="residue", how="left")
df = df.merge(conservation, on="mut_id", how="left")
df = df.merge(stability, on="mut_id", how="left")
def classify(row):
    mut = row["mutation_RBD"]

    # Rule 1: deleterious if low expr or low ace2
    if (row["expr_score"] < EXPR_DELETERIOUS) or (row["ace2_score"] < ACE2_DELETERIOUS):
        return 0  # deleterious

    # Rule 2: beneficial if strong ace2 + decent expr
    if (row["ace2_score"] > ACE2_BENEFICIAL) and (row["expr_score"] > EXPR_MIN_BENEFICIAL):
        return 2  # beneficial

    # Rule 3: stability penalty
    if not pd.isna(row.get("ddG")) and row["ddG"] > DDG_DESTABILIZING:
        return 0  # deleterious

    # Rule 4: conservation penalty
    if not pd.isna(row.get("conservation")) and row["conservation"] > CONSERVED:
        return 0  # deleterious

    # Rule 5: buried mutations (low RSA) more likely deleterious
    if not pd.isna(row.get("rsa")) and row["rsa"] < RSA_BURIED:
        return 0  # deleterious

    # Rule 6: known Variants of Concern
    if mut in VOC_MUTS:
        return 2  # beneficial

    # otherwise neutral
    return 1
df["refined_label"] = df.apply(classify, axis=1)
df.to_csv("data/processed/mutation_features_refined.csv", index=False)

print(df["refined_label"].value_counts())
