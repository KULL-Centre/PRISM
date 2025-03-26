from Bio import PDB
import os
from ptm_dict import modres

def clean_pdb(input_pdb, output_pdb, chains_to_keep='AB', keep_ligands=False, ligands_to_keep=None, ptm_mode='keep'):
    global modres
    """
    Cleans and restructures a PDB file, renumbers residues, and reorders ligands at the end.

    Parameters:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to save the cleaned PDB file.
        chains_to_keep (list): List of chain IDs to retain.
        keep_ligands (bool): Whether to keep ligands (small molecules).
        ligands_to_keep (list): List of ligand residue names to keep (if keep_ligands=True).
        ptm_mode (str): Mode for handling PTMs ('keep', 'skip', 'reverse').
    """

    # A dict to convert three letter code to one letter code
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    modres_new = modres.copy()
    
    for aa in d3to1.keys():
        modres_new.pop(aa)
        
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    
    clean_structure = PDB.Structure.Structure("cleaned_protein")
    model = structure[0]  # Most PDBs have a single model
    clean_model = PDB.Model.Model(0)

    new_atom_serial = 1  # Continuous atom numbering
    new_residue_id = 1   # Continuous residue numbering
    protein_chains = []
    ligand_residues = []
    fasta_sequences = {}  # Store sequences per chain

    if ligands_to_keep is None:
        ligands_to_keep = []
    if modres_new is None:
        modres_new = {}

    for chain in model:
        if chain.id in chains_to_keep:
            clean_chain = PDB.Chain.Chain(chain.id)
            sequence = list()
            for residue in chain:
                resname = residue.get_resname()                
                het_flag = residue.id[0]  # Identifies if it's a HETATM
                
                # Remove water molecules
                if het_flag == "W":
                    continue  

                # Handle PTMs based on user selection
                if het_flag != " " and resname in modres_new:
                    if ptm_mode == 'keep':
                        het_flag = " "  # Treat PTMs as standard residues
                        sequence.append(d3to1[modres_new[resname]]) # Append FASTA sequence
                    elif ptm_mode == 'skip':
                        continue  # Remove PTMs
                    elif ptm_mode == 'reverse':
                        het_flag = " "  # Convert PTM to standard residue
                        residue.resname = modres_new[resname]  # Rename residue
                        sequence.append(d3to1[modres_new[resname]]) # Append FASTA sequence
                elif het_flag == " ":
                    sequence.append(d3to1[resname]) # Append FASTA sequence

                # Separate ligands for reordering
                if het_flag != " ":
                    if not keep_ligands or resname not in ligands_to_keep:
                        continue  # Skip unwanted ligands
                    ligand_residues.append(residue.copy())  # Store for later
                    continue

                # Renumber residue ID
                residue.id = (" ", new_residue_id, " ")
                new_residue_id += 1

                # Renumber atoms
                for atom in residue:
                    atom.serial_number = new_atom_serial
                    new_atom_serial += 1

                clean_chain.add(residue.copy())  # Add residue to clean chain

            protein_chains.append(clean_chain)  # Store cleaned chain

            # Create fasta file
            case_id = output_pdb.split('/')[-1].split('.')[0]
            fasta_path = os.path.join('/'.join(output_pdb.split('/')[0:-1]), f'{case_id}_{chain.id}.fasta')
            with open(fasta_path, 'w') as fp:
                fp.write(f'>input_{chain.id}\n')
                fp.write(''.join(sequence))

    # Renumber and append ligands at the end
    for residue in ligand_residues:
        residue.id = ("HETATM", new_residue_id, " ")
        new_residue_id += 1

        for atom in residue:
            atom.serial_number = new_atom_serial
            new_atom_serial += 1

    # Add chains and ligands to the model
    for chain in protein_chains:
        clean_model.add(chain)
    if ligand_residues:
        ligand_chain = PDB.Chain.Chain("L")  # Create a separate ligand chain
        for ligand in ligand_residues:
            ligand_chain.add(ligand)
        clean_model.add(ligand_chain)

    clean_structure.add(clean_model)

    # Save cleaned structure
    io = PDB.PDBIO()
    io.set_structure(clean_structure)
    io.save(output_pdb)

    print(f"Cleaned PDB file saved as: {output_pdb}")

# Example usage
if __name__ == "__main__":
    input_pdb = '/lustre/hpc/sbinlab/panf/Cartesian_pipeline_runs/test_Beatriz/1KNE.pdb'
    output_pdb = "1KNE_clean.pdb"
    chains_to_keep = ["A", "P"]  # Chains to retain
    keep_ligands = True  # Keep ligands?
    ligands_to_keep = ["017", "PO4"]  # List of ligands to keep
    ptm_mode = "reverse"  # Change to 'keep', 'skip', or 'reverse'

    clean_pdb(input_pdb, output_pdb, chains_to_keep, keep_ligands, ligands_to_keep, ptm_mode)