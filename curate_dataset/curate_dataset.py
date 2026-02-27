import os
import requests
import pandas as pd
import gemmi
import freesasa
from Bio.PDB import PDBList
from utils import get_chain_seq

BSA_THRESHOLD = 400.0
CANDIDATE_ID_CSV = 'protein-peptide_complexes_candidates.csv'
OUTPUT_CSV = 'protein_peptide_dataset.csv'
SAVE_DIR = './dataset_cif_files'

def get_pdb_complexes():
    """Queries RCSB for high-res protein-peptide complexes."""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal", "value": 2.5}},
                {"type": "terminal", "service": "text", "parameters": {"attribute": "entity_poly.rcsb_sample_sequence_length", "operator": "greater", "value": 30}},
                {"type": "terminal", "service": "text", "parameters": {"attribute": "entity_poly.rcsb_sample_sequence_length", "operator": "less_or_equal", "value": 25}},
                {"type": "terminal", "service": "text", "parameters": {"attribute": "entity_poly.rcsb_entity_polymer_type", "operator": "exact_match", "value": "Protein"}}
            ]
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }
    response = requests.post(url, json=query)
    if response.status_code == 200:
        return [hit['identifier'] for hit in response.json().get('result_set', [])]
    return []

def identify_peptide_chains(model):
    unique_peptides = []
    p_ent = set()

    for chain in model:
        poly = chain.get_polymer()
        if len(poly) == 0 or len(poly) > 25:
            continue
        p_type = poly.check_polymer_type()
        if not(p_type == gemmi.PolymerType.PeptideL or p_type == gemmi.PolymerType.PeptideD):
            continue

        # if same entity as encountered before, skip
        ent_id = chain[0].entity_id
        if ent_id in p_ent:
            continue

        # get sequence (non-standard residue is marked as 'X')
        seq = poly.make_one_letter_sequence()
        seq = "".join('X' if char.islower() else char for char in seq)

        standard_res_count = sum(1 for char in seq if char != 'X')

        if standard_res_count >= 5:
            p_ent.add(ent_id)
            unique_peptides.append(chain)
        
    return unique_peptides

def calculate_interface_data(cif_path):
    pdb_id = os.path.basename(cif_path).split('-')[0]
    structure = gemmi.read_structure(cif_path)

    # remove water and hydrogen
    structure.remove_waters()
    structure.remove_hydrogens() 
    model = structure[0]

    classifier = freesasa.Classifier.getStandardClassifier('protor')
    
    peptides = identify_peptide_chains(model)
    
    if not peptides: 
        print(f"\nNo valid peptide found for {pdb_id}.")
        return []
    
    print(f"\n{pdb_id}: Peptides: {[c.name for c in peptides]}")

    found_pairs = []

    all_coords = []
    all_radii = []

    peptide_chain_start_ends = {}

    for chain in model:
        if chain in peptides:
            start = len(all_radii)
        for res in chain:
            for atom in res:
                all_coords.extend(atom.pos.tolist())
                all_radii.append(classifier.radius(res.name, atom.name))

        if chain in peptides:
            peptide_chain_start_ends[chain.name] = (start, len(all_radii))
    
    if not all_coords:
        print(f"No valid atoms found in structure {pdb_id}.")
        return []
    
    # Full complex SASA
    full_complex_result = freesasa.calcCoord(all_coords, all_radii)

    for p_chain in peptides:

        # Extract isolated peptide coords/radii
        start_pos, end_pos = peptide_chain_start_ends[p_chain.name]
        
        if start_pos == end_pos:
            print(f"No coordinates found for peptide chain {p_chain.name} in {pdb_id}.")
            continue
        
        p_coords = all_coords[start_pos*3:end_pos*3]
        p_radii = all_radii[start_pos:end_pos]
        # isolated peptide SASA
        sasa_p_iso = freesasa.calcCoord(p_coords, p_radii).totalArea()

        # Extract all other chains
        other_coords = all_coords[:start_pos*3] + all_coords[end_pos*3:]
        other_radii = all_radii[:start_pos] + all_radii[end_pos:]
        if not other_coords:
            print(f"No other chains coordinates found for {pdb_id} for peptide {p_chain.name}.")
            continue
        # isolated other chains SASA
        sasa_other_iso = freesasa.calcCoord(other_coords, other_radii).totalArea()

        bsa = (sasa_other_iso + sasa_p_iso) - full_complex_result.totalArea()

        print(f"SASA peptide {p_chain.name} (isolated): {sasa_p_iso}")
        print(f"SASA other chains (isolated): {sasa_other_iso}")
        print(f"SASA complex: {full_complex_result.totalArea()}")
        print(f"Calculated BSA for {p_chain.name}: {bsa}")

        if bsa >= BSA_THRESHOLD:
            found_pairs.append({
                'pdb_id': pdb_id.upper(),
                'Peptide': p_chain.name,
                'Buried surface area': round(bsa, 3),
                'Sequence': get_chain_seq(p_chain),
            })
    return found_pairs

if __name__ == "__main__":
    if os.path.exists(CANDIDATE_ID_CSV):
        ids = pd.read_csv(CANDIDATE_ID_CSV)['pdb_id'].tolist()
        print(f"Loaded {len(ids)} candidate IDs from {CANDIDATE_ID_CSV}.")
    else:
        ids = get_pdb_complexes()
        print(f"Found {len(ids)} candidates. Filtering...")

        # save csv 
        df = pd.DataFrame(ids, columns=['pdb_id'])
        df.to_csv(CANDIDATE_ID_CSV, index=False)

    # download files
    os.makedirs(SAVE_DIR, exist_ok=True)
    pdbl = PDBList()
    for pdb_id in ids:
        try:
            if os.path.exists(os.path.join(SAVE_DIR, f"{pdb_id.lower()}-assembly1.cif")):
                print(f"Skipping {pdb_id} (already exists)")
                continue

            # Using BioPython to download the Biological Assembly (Assembly 1)
            downloaded_path = pdbl.retrieve_assembly_file(
                pdb_id, pdir=SAVE_DIR, file_format="mmCif", assembly_num=1
            )

        except Exception as e:
            print(f"Error downloading {pdb_id}: {e}")

    # process each CIF file
    final_list = []
    for pdb_id in ids:
        cif_path = os.path.join(SAVE_DIR, f"{pdb_id.lower()}-assembly1.cif")
        if not os.path.exists(cif_path):
            print(f"File not found: {cif_path}, skipping.")
            continue

        data = calculate_interface_data(cif_path)

        if data:
            final_list.extend(data)
        else:
            # delete the cif files that do not qualify for the dataset
            os.remove(cif_path)
            print(f"Removed {cif_path} as it does not meet criteria.")

    # save dataset info to csv
    df = pd.DataFrame(final_list)
    df.to_csv(OUTPUT_CSV, index=False)

    print(f"\nDone! Saved {len(df)} entries to {OUTPUT_CSV}.") 
    