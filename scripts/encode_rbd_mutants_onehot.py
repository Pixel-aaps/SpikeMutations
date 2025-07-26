import os
import pandas as pd
from Bio import SeqIO
import numpy as np

fasta_path = "data/interim/mutated_variants.fasta"
output_path = "data/processed/encoded_sequences.csv"
rbd_start = 319 - 1
rbd_end = 541

# One-Hot Encoding Map
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
aa_to_index = {aa: i for i, aa in enumerate(amino_acids)}

def one_hot_encode(seq):
    encoded = np.zeros((len(seq), 20))
    for i, aa in enumerate(seq):
        if aa in aa_to_index:
            encoded[i, aa_to_index[aa]] = 1
        else:
            pass  # Unknown AA → leave as zero
    return encoded.flatten()

# Read mutated sequences
records = list(SeqIO.parse(fasta_path, "fasta"))
encoded_data = []
ids = []

for rec in records:
    full_seq = str(rec.seq)
    rbd_seq = full_seq[rbd_start:rbd_end]
    encoded_vector = one_hot_encode(rbd_seq)
    encoded_data.append(encoded_vector)
    ids.append(rec.id)

# Create DataFrame
column_names = [f"{pos}_{aa}" for pos in range(1, 224) for aa in amino_acids]
df = pd.DataFrame(encoded_data, columns=column_names)
df.insert(0, "Mutation_ID", ids)

# Save CSV
os.makedirs(os.path.dirname(output_path), exist_ok=True)
df.to_csv(output_path, index=False)

print(f"Saved encoded features to {output_path}")
