from Bio.PDB import PDBParser, PDBIO, MMCIFParser
import string
import warnings
from Bio import BiopythonWarning
import gemmi


PROTEIN_POLYMERS = [
    "PeptideL",
    "PeptideD",
]

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

def extract_structure_from_cif(cif_path: str, assembly_idx=0):
    """
    Extract PDB structure from mmCIF

    Args:
        script_path: Path to the PDB files.
        assembly_idx: Index of the assembly to transform to, default is 0 (first assembly).
    Returns:
        tuple: A tuple containing the PDB ID and the structure object if it exists.
    """
    try:
        doc = gemmi.cif.read_file(cif_path)
        block = doc.sole_block()
        structure = gemmi.make_structure_from_block(block)

        if structure.assemblies:
            # if the structure has assemblies, by default, transform to the first assembly
            assembly_name = structure.assemblies[assembly_idx].name
            how = gemmi.HowToNameCopiedChain.AddNumber
            structure.transform_to_assembly(assembly_name, how)
            
        return structure

    except Exception as e:
        print(f"Error processing {cif_path}: {e}")
        return None

def parse_cif_and_clean(input_cif, output_pdb, peptide_chain):
    """ 
    Cleans the input CIF file by reordering chains (peptide chain that requires diffusion are placed first) 
    and renumbering residues to non-negative integers.

    Args:
        input_cif (str): Path to the input CIF file.
        output_pdb (str): Path to save the cleaned PDB file.
        peptide_chain (str): Chain ID of the peptide that requires diffusion.
    """

    # Suppress Biopython warnings about ID collisions
    warnings.filterwarnings("ignore", category=BiopythonWarning)

    print(f"Parsing CIF and cleaning: {input_cif} -> {output_pdb} with peptide chain {peptide_chain}")

    old_to_new_chains_mapping = {}

    structure_gemmi = extract_structure_from_cif(input_cif)

    chain_label_id_to_auth_id = {}
    for chain in structure_gemmi[0]:  # Assuming we are working with the first model
        for sub in chain.subchains():
            chain_label_id_to_auth_id[sub.subchain_id()] = chain.name

    print(f"Chain label ID to Author ID mapping: {chain_label_id_to_auth_id}")

    # Build entity mapping to identify which chains belong to the same entity 
    entities = {}
    k = 1
    for entity in structure_gemmi.entities:
        print(f"Processing entity: {entity.name}, polymer type: {entity.polymer_type.name}, subchains: {entity.subchains}")
        # maps chain name to entity object, skips over non-protein entities
        if entity.polymer_type.name not in PROTEIN_POLYMERS:
            continue
        k_entity = None
        # if this entity is the same as the chain of interest's entity
        for sc in entity.subchains:
            if chain_label_id_to_auth_id[sc].upper() == peptide_chain.upper():
                print(f"Found peptide chain {peptide_chain} belongs to entity {entity.name}")
                k_entity = 0
                break
        if k_entity is None:
            k_entity = k
            k += 1
        for chain_id_ in entity.subchains:
            entities[chain_id_] = k_entity # k maps to entity_idx

    print(f"Entity mapping keys: {entities.keys()}")
    
    # Now reorder chains (place peptide chain first, then sort the rest by original chain ID) and renumber residues
    parser = MMCIFParser(auth_chains=False) 
    structure = parser.get_structure("struct", filename=input_cif) # biopython structure for easier manipulation of chains and residues

    for model in structure:
        chains = list(model.get_chains())
        
        sorted_chains = sorted(
            chains,  
            # when peptide chain extracted in ./curate_dataset/curate_dataset.py, auth_asym_id was used
            # if the chain's auth id matches the peptide chain, and it is a peptide chain (exists in entity mapping), then place it first, otherwise sort by original chain ID
            # chains that are not in the entity mapping (e.g. water, ligand) are placed at the end of the sorted list, and sort them by original chain ID as well
            key=lambda c: (0 if chain_label_id_to_auth_id[c.id] == peptide_chain and c.id in entities else 1 if c.id in entities else 2, c.id)
        )
        print(f"Old sorted chain label ids: {[c.id for c in sorted_chains]}")

        # Reorder and Rename Chains in the Model
        # Chains must be detached before re-adding to maintain a clean hierarchy
        for c in chains:
            model.detach_child(c.id)
        
        # available chain labels (the same as in RFdiffusion)
        available_chains = sorted(list(set(string.ascii_uppercase)))

        for i, chain in enumerate(sorted_chains):
            if (sorted_chains[i].id) not in entities:
                continue # skip non-protein chains that are not in the entity mapping (e.g. water, ligand)
            
            # Assign new sequential chain ID (A, B, C...)
            if i < len(available_chains):
                new_id = available_chains[i]
            else:
                return None, None
            
            old_to_new_chains_mapping[sorted_chains[i].id] = new_id
            
            chain.id = new_id
            model.add(chain)
            
            # Reassign residue indices to remove insertion codes (icodes)
            # This ensures residues with insertion codes get proper sequential indices
            residue_list = list(chain.get_residues())
            for idx, res in enumerate(residue_list, start=1):
                res.id = (res.id[0], idx, ' ')

    new_entities = {}
    # check that mappings are consistent with entity mapping
    for old_chain_id, new_chain_id in old_to_new_chains_mapping.items():
        if old_chain_id not in entities:
            print(f"Warning: Old chain ID {old_chain_id} not found in entity mapping.")
        else:
            # switch the old chain ID to the new chain ID in the entity mapping
            k_entity = entities[old_chain_id]
            new_entities[new_chain_id] = k_entity

    # Save the modified structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    print(f"Cleaned PDB saved to {output_pdb}")

    return old_to_new_chains_mapping, new_entities



def count_chains_in_cif(cif_path: str) -> int:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('struct', cif_path)
    first_model = next(structure.get_models())
    return sum(1 for _ in first_model.get_chains())