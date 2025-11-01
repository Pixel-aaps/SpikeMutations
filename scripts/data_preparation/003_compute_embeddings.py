# scripts/03_compute_embeddings_gpu.py

import os
import torch
import pandas as pd
from tqdm import tqdm
import esm
import numpy as np

# ==============================
# Paths
# ==============================

labels_file = "/content/labels.csv"

# Load labeled mutations
df = pd.read_csv(labels_file)

# ==============================
# Load pretrained ESM-2 model
# ==============================
print("Loading ESM-2 650M model...")
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()

# Send model to GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")
model = model.to(device)

# ==============================
# Define WT Spike RBD (319–541)
# ==============================
WT_RBD = (
    "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYAD"
    "SFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGST"
    "PCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"
)

# ==============================
# Mutation function
# ==============================
def apply_mutation(seq, mut):
    """Apply a single mutation like 'N501Y'. Returns (mut_seq, pos)."""
    wt = mut[0]
    mt = mut[-1]
    pos = int(mut[1:-1]) - 319  # adjust to RBD numbering

    if mt == "*":
        raise ValueError(f"Stop codon {mut} encountered")

    if seq[pos] != wt:
        raise AssertionError(f"Mismatch at {mut}: expected {wt}, got {seq[pos]}")

    return seq[:pos] + mt + seq[pos+1:], pos

# ==============================
# Embedding loop
# ==============================
embeddings = {}

print("Encoding mutations with ESM-2 650M...")
for _, row in tqdm(df.iterrows(), total=len(df)):
    mut_id = row["mut_id"]

    try:
        seq, pos = apply_mutation(WT_RBD, mut_id)
    except (AssertionError, ValueError) as e:
        print(f"Skipping {mut_id}: {e}")
        # Save NaN embedding for skipped mutation
        embeddings[mut_id] = np.full((1280,), np.nan)  # 1280 dims for ESM-2 650M
        continue

    # Convert to tokens
    batch_labels, batch_strs, batch_tokens = batch_converter([("mut", seq)])
    batch_tokens = batch_tokens.to(device)

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=False)

    token_representations = results["representations"][33]

    # Save embedding for mutated residue
    emb = token_representations[0, pos+1, :].cpu().numpy()  # +1 for [CLS]
    embeddings[mut_id] = emb

# ==============================
# Save embeddings
# ==============================
emb_df = pd.DataFrame.from_dict(embeddings, orient="index")
emb_df.index.name = "mut_id"
emb_csv_path = "/content/embeddings.csv"
emb_df.to_csv(emb_csv_path)

print(f"Saved embeddings for {len(emb_df)} mutations → {emb_csv_path}")