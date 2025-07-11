from rdkit import Chem
from rdkit.Chem import AllChem
from vina import Vina
import os

def smiles_to_pdbqt(smiles, pdbqt_path):
    # Generate 3D conformation from SMILES
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    
    # Save as MOL file (temporary)
    mol_file = pdbqt_path.replace('.pdbqt', '.mol')
    Chem.MolToMolFile(mol, mol_file)

    # Convert to PDBQT using Open Babel
    os.system(f'obabel {mol_file} -O {pdbqt_path}')
    os.remove(mol_file)

def prepare_receptor_pdbqt(pdb_path, pdbqt_path):
    # Convert receptor PDB to PDBQT using Open Babel
    os.system(f'obabel {pdb_path} -O {pdbqt_path}')
    
def run_vina(receptor_pdbqt, ligand_pdbqt, center, size):
    v = Vina(sf_name='vina')

    # Set receptor and ligand
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)

    # Set box
    v.compute_vina_maps(center=center, box_size=size)

    # Run docking
    v.dock(exhaustiveness=8, n_poses=1)
    affinity = v.score()
    print(affinity)
    print(f'Predicted binding affinity: {affinity:.2f} kcal/mol')

# === USER INPUT ===
smiles = "CCO"  # replace with your SMILES
receptor_pdb = "receptor.pdb"  # replace with your receptor PDB path

# === FILE PREP ===
ligand_pdbqt = "ligand.pdbqt"
receptor_pdbqt = "receptor.pdbqt"

smiles_to_pdbqt(smiles, ligand_pdbqt)
prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt)

# === BOX SETUP ===
# You must choose this based on known site or guess. Example:
center = [0, 0, 0]
size = [20, 20, 20]

run_vina(receptor_pdbqt, ligand_pdbqt, center=center, size=size)
