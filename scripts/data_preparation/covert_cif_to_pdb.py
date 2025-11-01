from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename="data/external/6M0J.cif")
PDBFile.writeFile(fixer.topology, fixer.positions, open("data/external/6M0J.pdb", "w"))

# import os

# pdb_file = "data/external/6M0J.pdb"
# clean_pdb_file = "data/external/6M0J_clean.pdb"

# # Read the original PDB and remove problematic lines (HETATM for metals like ZN)
# with open(pdb_file, "r") as f:
#     lines = []
#     for l in f:
#         if l.startswith("HETATM") and "ZN" in l:
#             continue  # skip zinc
#         lines.append(l)

# # Save the cleaned PDB
# os.makedirs(os.path.dirname(clean_pdb_file), exist_ok=True)
# with open(clean_pdb_file, "w") as f:
#     f.writelines(lines)

# print(f"Cleaned PDB saved → {clean_pdb_file}")

import os

pdb_file = "data/external/6M0J_clean.pdb"
fully_clean_pdb = "data/external/6M0J_clean_no_het.pdb"

with open(pdb_file, "r") as f:
    lines = []
    for l in f:
        # Keep only standard protein residues (ATOM) and ignore HETATM (Zn, Cl, water, etc.)
        if l.startswith("HETATM"):
            continue
        if l.startswith("ATOM"):
            lines.append(l)
        if l.startswith("TER") or l.startswith("END"):
            lines.append(l)

os.makedirs(os.path.dirname(fully_clean_pdb), exist_ok=True)
with open(fully_clean_pdb, "w") as f:
    f.writelines(lines)

print(f"Fully cleaned PDB (no HETATM) saved → {fully_clean_pdb}")
