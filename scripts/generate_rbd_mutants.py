from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

fasta_path = "data/raw/spike_protein.fasta"
output_path = "data/interim/mutated_variants.fasta"
rbd_start = 319 - 1  # 0-based indexing
rbd_end = 541        # inclusive end position
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")  # 20 standard AAs

# Load spike protein 
record = SeqIO.read(fasta_path, "fasta")
full_seq = str(record.seq)
rbd_seq = full_seq[rbd_start:rbd_end]

# Create mutated sequences 
mutant_records = []
for i, original_aa in enumerate(rbd_seq):
    position = rbd_start + i + 1  # convert back to 1-based position
    for mutant_aa in amino_acids:
        if mutant_aa == original_aa:
            continue
        mutated_seq = list(full_seq)
        mutated_seq[rbd_start + i] = mutant_aa
        mutated_seq_str = ''.join(mutated_seq)

        mutant_id = f"mut_pos{position}_{original_aa}->{mutant_aa}"
        rec = SeqRecord(Seq(mutated_seq_str), id=mutant_id, description="")
        mutant_records.append(rec)

# Save to FASTA
os.makedirs(os.path.dirname(output_path), exist_ok=True)
SeqIO.write(mutant_records, output_path, "fasta")

print(f"Saved {len(mutant_records)} RBD mutant sequences to {output_path}")
