from Bio.PDB import PDBParser, PDBIO
import string
import warnings
from Bio import BiopythonWarning

def make_residues_index_non_negative(input_pdb, output_pdb):
    """ 
    Shifts all residue numbers in the input PDB to be non-negative and saves the modified structure.
    """

    # Suppress Biopython warnings about ID collisions
    warnings.filterwarnings("ignore", category=BiopythonWarning)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', input_pdb)
    
    # Find the minimum residue number to determine the shift needed
    for model in structure:
        for chain in model:
            neg_res = [res.id[1] for res in chain if res.id[1] < 0]
            if not neg_res:
                continue
            min_res = min(neg_res)
            print(min_res)
            shift = abs(min_res) + 1
            print(f"Shifting residues in chain {chain.id} by {shift}.")
            
            for residue in chain:
                # Get the current ID components
                res_id = list(residue.id)
                # Update the sequence number (the second element in the tuple)
                res_id[1] = residue.id[1] + shift
                # Apply the new ID as a tuple
                residue.id = tuple(res_id)

    # Save the modified structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb) 

def clean_pdb(input_pdb, output_pdb, peptide_chain):
    """ 
    Cleans the input PDB file by reordering chains (peptide chain that requires diffusion are placed first) 
    and renumbering residues to non-negative integers.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to save the cleaned PDB file.
        peptide_chain (str): Chain ID of the peptide that requires diffusion.
    """

    # Suppress Biopython warnings about ID collisions
    warnings.filterwarnings("ignore", category=BiopythonWarning)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', input_pdb)

    new_pep_chains_mapping = {}
    
    for model in structure:
        chains = list(model.get_chains())
        
        sorted_chains = sorted(
            chains, 
            key=lambda c: (0, '') if c.id == peptide_chain else (1, c.id)
        )

        # print(f"chains: {[c.id for c in sorted_chains]}")
        
        # Reorder and Rename Chains in the Model
        # Chains must be detached before re-adding to maintain a clean hierarchy
        for c in chains:
            model.detach_child(c.id)
        
        # available chain labels (the same as in RFdiffusion)
        available_chains = sorted(list(set(string.ascii_uppercase)))

        for i, chain in enumerate(sorted_chains):
            # Assign new sequential chain ID (A, B, C...)
            new_id = available_chains[i] if i < len(available_chains) else f"C{i}"
            # save the new peptide chain IDs 
            if sorted_chains[i].id == peptide_chain:
                new_pep_chains_mapping[sorted_chains[i].id] = new_id
            
            chain.id = new_id
            model.add(chain)
            
            # # Renumber Residues to Non-Negative (Two-Pass Method)
            # # Find the minimum sequence number to calculate the shift
            # all_res_ids = [res.id[1] for res in chain]
            # min_val = min(all_res_ids)
            
            # if min_val < 1:
            #     shift = abs(min_val) + 1
            #     # Pass 1: Shift to a temporary high range to prevent sibling ID conflicts
            #     for res in chain:
            #         res.id = (res.id[0], res.id[1] + 10000, res.id[2])
            #     # Pass 2: Apply the final non-negative sequential ID
            #     for res in chain:
            #         res.id = (res.id[0], res.id[1] - 10000 + shift, res.id[2])
            
            # Reassign residue indices to remove insertion codes (icodes)
            # This ensures residues with insertion codes get proper sequential indices
            residue_list = list(chain.get_residues())
            for idx, res in enumerate(residue_list, start=1):
                res.id = (res.id[0], idx, ' ')

    # Save the modified structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    print(f"Cleaned PDB saved to {output_pdb}")

    return new_pep_chains_mapping



