# # make_mutants_quick.py
# from Bio.PDB import PDBParser, PDBIO
# import copy
# import os
# import pandas as pd

# one2three = {
#  'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU',
#  'G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE',
#  'P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'
# }

# # Input PDB and mutation list
# parser = PDBParser(QUIET=True)
# structure = parser.get_structure("wt", "notebooks/6M0J_clean_h.pdb")

# df = pd.read_csv("data/processed/mutation_list.csv")
# mutations = df["mutation"].dropna().tolist()

# # Output folder (absolute or relative)
# output_dir = "data/processed/mutants_quick"
# os.makedirs(output_dir, exist_ok=True)

# # Generate mutants
# for mut in mutations:
#     wt = mut[0]; new = mut[-1]; pos = int(mut[1:-1])
#     st_copy = copy.deepcopy(structure)  # fresh copy for each mutation
#     for model in st_copy:
#         for chain in model:
#             for res in chain:
#                 if res.id[1] == pos:
#                     res.resname = one2three[new]
#     outp = os.path.join(output_dir, "6m0j_" + mut + ".pdb")
#     io = PDBIO()
#     io.set_structure(st_copy)
#     io.save(outp)
#     print("Saved", outp)
import os
import subprocess

# Paths
receptor = "6m0j_wt.pdbqt"
ligand_dir = "notebooks/ligands_ready_docking"
output_dir = "vina_wt_out2"

# Grid center and size
CX, CY, CZ = -20.475, 13.629, -25.369
SX, SY, SZ = 24, 24, 24
exhaustiveness = 8

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Get all pdbqt ligands
ligands = [f for f in os.listdir(ligand_dir) if f.endswith(".pdbqt")]

for lig in ligands:
    lig_path = os.path.join(ligand_dir, lig)
    name = os.path.splitext(lig)[0]
    out_file = os.path.join(output_dir, "{}_out.pdbqt".format(name))
    log_file = os.path.join(output_dir, "{}_wt.log".format(name))  # log file path

    cmd = [
        "vina",
        "--receptor", receptor,
        "--ligand", lig_path,
        "--center_x", str(CX),
        "--center_y", str(CY),
        "--center_z", str(CZ),
        "--size_x", str(SX),
        "--size_y", str(SY),
        "--size_z", str(SZ),
        "--exhaustiveness", str(exhaustiveness),
        "--out", out_file
    ]

    print("Running docking for {}...".format(name))
    
    # Redirect stdout and stderr to the log file
    with open(log_file, "w") as logfile:
        result = subprocess.run(cmd, stdout=logfile, stderr=logfile, text=True)

    if result.returncode != 0:
        print("❌ Docking failed for {}, check {}".format(name, log_file))
    else:
        print("✅ Docking finished for {}, results saved to {}".format(name, out_file))
