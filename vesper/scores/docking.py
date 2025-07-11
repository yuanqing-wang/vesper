from rdkit import Chem
from rdkit.Chem import AllChem
from vina import Vina
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import tempfile

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
    os.system(f'prepare_receptor -r {pdb_path} -o {pdbqt_path}')
    
def get_pdb(pdb_code, pdb_path):
    # Download PDB file if not exists
    if not os.path.exists(pdb_path):
        os.system(f'wget -O {pdb_path} https://files.rcsb.org/download/{pdb_code}.pdb')
        
def clean_pdb_with_pdbfixer(input_pdb, output_pdb):
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()
    # Remove water and ions
    chains_to_remove = [chain.index for chain in fixer.topology.chains() if any(res.name in ['HOH', 'WAT', 'NA', 'K', 'CL', 'CA', 'MG', 'ZN', 'SO4', 'PO4'] for res in chain.residues())]
    fixer.removeChains(chains_to_remove)
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    
def run_vina(pdb_code, smiles, center=[0.0, 0.0, 0.0], size=[20.0, 20.0, 20.0]):
    # Download PDB if needed and clean it up
    # get_pdb(pdb_code, pdb_code + ".pdb")
    # clean_pdb_with_pdbfixer(pdb_code + ".pdb", "clean_" + pdb_code + ".pdb")

    with tempfile.TemporaryDirectory() as tmpdir:
        receptor_pdbqt = os.path.join(tmpdir, "receptor.pdbqt")
        ligand_pdbqt = os.path.join(tmpdir, "ligand.pdbqt")
        prepare_receptor_pdbqt("clean_" + pdb_code + ".pdb", receptor_pdbqt)
        smiles_to_pdbqt(smiles, ligand_pdbqt)

        v = Vina(sf_name='vina')
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(center=center, box_size=size)
        v.dock(exhaustiveness=8, n_poses=1)
        affinity = v.score()[0]
        print(affinity)

if __name__ == "__main__":
    # Example usage
    pdb_code = "1a4w"  # Replace with your PDB code
    smiles = "CCO"     # Replace with your ligand SMILES
    run_vina(pdb_code, smiles)
