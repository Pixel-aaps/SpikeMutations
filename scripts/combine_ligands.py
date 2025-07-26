from rdkit import Chem
import os

# Path where your individual SDFs are stored
input_folder = "data/raw/3D_temp_sdf/"
output_file = "data/raw/drug_ligands.sdf"

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Prepare writer for combined SDF
writer = Chem.SDWriter(output_file)

# Loop through all .sdf files in the folder
for filename in os.listdir(input_folder):
    if filename.endswith(".sdf"):
        filepath = os.path.join(input_folder, filename)
        suppl = Chem.SDMolSupplier(filepath)

        for mol in suppl:
            if mol:
                # Set name to filename (without extension)
                mol.SetProp("_Name", os.path.splitext(filename)[0])
                writer.write(mol)
                print(f"Added: {filename}")
            else:
                print(f"Skipped: {filename} (invalid structure)")

writer.close()
print(f"\n🎉 Combined SDF saved to: {output_file}")
