# # def fix_pdb(input_pdb, output_pdb):
# #     with open(input_pdb, "r") as f:
# #         lines = f.readlines()

# #     fixed_lines = ["HEADER    FOLDX REPAIRED STRUCTURE\n"]
# #     atom_counter = 1

# #     for line in lines:
# #         if line.startswith(("ATOM", "HETATM")):
# #             # overwrite atom serial (cols 7-11)
# #             fixed_line = f"{line[:6]}{atom_counter:5d}{line[11:]}"
# #             fixed_lines.append(fixed_line)
# #             atom_counter += 1
# #         else:
# #             fixed_lines.append(line)

# #     with open(output_pdb, "w") as f:
# #         f.writelines(fixed_lines)

# #     print(f"Fixed PDB written to {output_pdb}")


# # # Usage
# # fix_pdb("data/external/fold_x/6m0j_Repair.pdb",
# #         "data/external/fold_x/6m0j_Repair_fixed.pdb")

# # from pdbfixer import PDBFixer
# # from openmm.app import PDBFile

# # fixer = PDBFixer(filename='notebooks/6M0J_clean.pdb')
# # fixer.findMissingResidues()
# # fixer.findMissingAtoms()
# # fixer.addMissingHydrogens(7.0)  # pH 7.0

# # PDBFile.writeFile(fixer.topology, fixer.positions, open('6M0J_clean_h.pdb', 'w'))
# # print("Hydrogens added.")
# # import os
# # import subprocess

# # input_folder = "data/raw/3D_temp_sdf"
# # output_folder = "ligands_ready_docking"
# # os.makedirs(output_folder, exist_ok=True)

# # for file in os.listdir(input_folder):
# #     if file.endswith(".sdf"):
# #         input_path = os.path.join(input_folder, file)
# #         base_name = os.path.splitext(file)[0]

# #         output_h_sdf = os.path.join(output_folder, f"{base_name}_3D_h.sdf")
# #         output_pdbqt = os.path.join(output_folder, f"{base_name}_3D_h.pdbqt")

# #         # Step 1: Add hydrogens at pH 7.4
# #         subprocess.run([
# #             "obabel", "-i", "sdf", input_path,
# #             "-o", "sdf", "-O", output_h_sdf, "-p", "7.4"
# #         ], check=True)

# #         # Step 2: Convert to PDBQT
# #         subprocess.run([
# #             "obabel", "-i", "sdf", output_h_sdf,
# #             "-o", "pdbqt", "-O", output_pdbqt
# #         ], check=True)

# #         print(f"Processed {file}: hydrogenated SDF and PDBQT created")

# from Bio.PDB import PDBParser
# import numpy as np
# import pandas as pd
# import freesasa  # Optional: for SASA calculation, install via pip

# # Load protein
# parser = PDBParser(QUIET=True)
# structure = parser.get_structure("spike", "/mnt/c/Users/avant/OneDrive/Desktop/VSCode/C++/SpikeMutations/notebooks/6M0J_clean_h.pdb")  # cleaned & hydrogenated

# # Define pocket residues (example RBD: 331-528)
# pocket_residues = [res for res in structure.get_residues() if 331 <= res.id[1] <= 528]

# # Residue properties
# hydrophobicity_scale = {'A': 1.8,'R': -4.5,'N': -3.5,'D': -3.5,'C': 2.5,'Q': -3.5,
#                         'E': -3.5,'G': -0.4,'H': -3.2,'I': 4.5,'L': 3.8,'K': -3.9,
#                         'M': 1.9,'F': 2.8,'P': -1.6,'S': -0.8,'T': -0.7,'W': -0.9,
#                         'Y': -1.3,'V': 4.2}

# polar_residues = {'R','K','H','D','E','N','Q','S','T','Y'}
# positive_residues = {'R','K','H'}
# negative_residues = {'D','E'}
# aromatic_residues = {'F','W','Y','H'}

# # Compute residue-level features
# hydrophobicity = [hydrophobicity_scale.get(res.get_resname()[0], 0) for res in pocket_residues]
# avg_hydro = np.mean(hydrophobicity)
# std_hydro = np.std(hydrophobicity)

# num_polar = sum(1 for res in pocket_residues if res.get_resname()[0] in polar_residues)
# frac_polar = num_polar / len(pocket_residues)

# num_positive = sum(1 for res in pocket_residues if res.get_resname()[0] in positive_residues)
# num_negative = sum(1 for res in pocket_residues if res.get_resname()[0] in negative_residues)
# frac_charge = (num_positive + num_negative) / len(pocket_residues)

# num_aromatic = sum(1 for res in pocket_residues if res.get_resname()[0] in aromatic_residues)
# frac_aromatic = num_aromatic / len(pocket_residues)

# pocket_size = len(pocket_residues)

# # Optional: compute SASA using freesasa
# try:
#     structure_str = "/mnt/c/Users/avant/OneDrive/Desktop/VSCode/C++/SpikeMutations/notebooks/6M0J_clean_h.pdb"
#     result = freesasa.Structure(structure_str)
#     sasa_result = freesasa.calc(result)
#     sasa = sasa_result.totalArea()
# except:
#     sasa = np.nan

# # Combine features
# protein_features = [avg_hydro, std_hydro, frac_polar, frac_charge, frac_aromatic, pocket_size, sasa]

# # Save to CSV
# columns = ["AvgHydro","StdHydro","FracPolar","FracCharge","FracAromatic","PocketSize","SASA"]
# df_protein = pd.DataFrame([protein_features], columns=columns)
# df_protein.to_csv("protein_features_rich.csv", index=False)

# print("Protein pocket features extracted (rich set).")

import difflib

file1 = r"C:\Users\avant\OneDrive\Desktop\VSCode\C++\SpikeMutations\data\processed\mutants_pdbqt\6m0j_A344C.pdbqt"
file2 = r"C:\Users\avant\OneDrive\Desktop\VSCode\C++\SpikeMutations\data\processed\mutants_pdbqt\6m0j_A344D.pdbqt"

# Read files
with open(file1, "r") as f:
    lines1 = f.readlines()

with open(file2, "r") as f:
    lines2 = f.readlines()

# Compare using difflib
diff = difflib.unified_diff(lines1, lines2, fromfile='A344C', tofile='A344D', lineterm='')

# Print differences
for line in diff:
    print(line)
