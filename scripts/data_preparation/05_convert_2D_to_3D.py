from rdkit import Chem
from rdkit.Chem import AllChem
import os

input_files = [
    "data/raw/2D_temp_sdf/2D_CamostatMethylate(5284360).sdf",
    "data/raw/2D_temp_sdf/2D_Ebselen(3194).sdf",
    "data/raw/2D_temp_sdf/2D_Ivermectin(6321424).sdf",
    "data/raw/2D_temp_sdf/2D_Lopinavir(92727).sdf",
    "data/raw/2D_temp_sdf/2D_Ritonavir(392622).sdf",
]

output_files = [
    "data/raw/3D_temp_sdf/CamostatMethylate.sdf",
    "data/raw/3D_temp_sdf/Ebselen.sdf",
    "data/raw/3D_temp_sdf/Ivermectin.sdf",
    "data/raw/3D_temp_sdf/Lopinavir.sdf",
    "data/raw/3D_temp_sdf/Ritonavir.sdf"
]

for input_file, output_file in zip(input_files, output_files):
    mol_supplier = Chem.SDMolSupplier(input_file)
    mol_writer = Chem.SDWriter(output_file)

    for mol in mol_supplier:
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) == 0:
            AllChem.UFFOptimizeMolecule(mol)
            mol_writer.write(mol)
            print(f"Converted 2D to 3D: {os.path.basename(output_file)}")
        else:
            print(f"Failed to embed: {os.path.basename(input_file)}")

    mol_writer.close()
