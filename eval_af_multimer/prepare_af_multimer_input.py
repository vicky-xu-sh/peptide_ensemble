import gemmi
import os

import argparse
import csv

def extract_structure_from_pdb(pdb_path: str, assembly_idx=0):
    """
    Extract PDB structure from pdb

    Args:
        pdb_path: Path to the PDB files in PDB format.
        assembly_idx: Index of the assembly to transform to, default is 0 (first assembly).
    Returns:
        tuple: A tuple containing the PDB ID and the structure object if it exists.
    """
    try:
        structure = gemmi.read_pdb(pdb_path)
        
        if structure.assemblies:
            # if the structure has assemblies, by default, transform to the first assembly
            assembly_name = structure.assemblies[assembly_idx].name
            how = gemmi.HowToNameCopiedChain.AddNumber
            structure.transform_to_assembly(assembly_name, how)
            return (os.path.basename(pdb_path).split(".pdb")[0], structure)
        else:
            return (os.path.basename(pdb_path).split(".pdb")[0], structure)

    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return (os.path.basename(pdb_path).split(".pdb")[0], None)

def write_fasta(path, name, seqs):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        f.write(":".join(seqs) + "\n")

def load_peptide_ensemble_info(path: str) -> dict[str, tuple[str, str]]:
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        field_map = {name.lower(): name for name in reader.fieldnames}
        if ("pdb_id" not in field_map) or ("peptide_chain_id" not in field_map) or ("peptide_sequence" not in field_map):
            raise ValueError(
                "CSV must contain 'PDB_ID', 'PEPTIDE_CHAIN_ID', and 'PEPTIDE_SEQUENCE' columns. "
                f"Found: {reader.fieldnames}"
            )
        pdb_col = field_map["pdb_id"]
        pep_chain_col = field_map["peptide_chain_id"]
        seq_col = field_map["peptide_sequence"]
        result: dict[str, tuple[str, str]] = {}
        for row in reader:
            pdb_id = row.get(pdb_col, "").strip().lower()
            seq = row.get(seq_col, "").strip().replace(" ", "").replace("\n", "")
            pep_chain = row.get(pep_chain_col, "").strip()
            if pdb_id in result:
                raise ValueError(f"Duplicate PDB ID found: {pdb_id}")
            result[pdb_id] = (pep_chain, seq)
        return result

def get_polymer_residues(chain):
    residue_span = []
    for res in chain:
        # skip water 
        if res.name in ('HOH', 'WAT'):
            continue

        # look up residue one-letter code, remove spaces
        code = gemmi.find_tabulated_residue(res.name).one_letter_code.strip()
        
        if code:
            residue_span.append(res)
            
    return residue_span

def write_full_fasta_sequence(
    pdb_id: str,
    chain_of_interest: str,
    structure: gemmi.Structure,
    peptide_designed_seqs: list[str],
    output_pdb_dir: str,
) -> None:
    """
    Extracts complete chain sequences of interest for a PDB structure, and replace designed sequences.
    
    """
    complete_seqs = []

    polymer_alignment_idx = None
    sequence_alignment_idx = None

    for chain in structure[0]:
        residue_span = get_polymer_residues(chain)

        chain_of_interest_flag = (
            chain.name.rstrip("1234567890").upper()
            == chain_of_interest.rstrip("1234567890").upper()
        )
        # print(f"chain {chain.name} chain of interest flag: {chain_of_interest_flag}")

        if chain_of_interest_flag:
            complete_seqs.append("*") # mark the chain of interest with a special character
            continue

        sequence = [
            gemmi.find_tabulated_residue(res.name).one_letter_code
            for res in residue_span
        ]
        raw_sequence = "".join(code for code in sequence)
        # print(f"raw_sequence: {raw_sequence}.")
        sequence = "".join((code if code.isupper() else "X") for code in sequence)

        if not sequence or set(sequence) == {"X"}:
            print(f"Warning: Chain {chain.name} has no recognizable residues. Skip.")
            continue

        complete_seqs.append(sequence)

    # print(f"complete_seqs: {complete_seqs}")

    # write fasta
    for idx, designed_seq in enumerate(peptide_designed_seqs):
        final_seqs = []
        for seq in complete_seqs:
            if seq == "*":
                final_seqs.append(designed_seq)
            else:
                final_seqs.append(seq)

        outname = f"{pdb_id}_{idx+1}.fasta"
        folder_path = os.path.join(output_pdb_dir,f"{pdb_id}_{idx+1}")
        os.makedirs(folder_path, exist_ok=True)
        fasta_path = os.path.join(folder_path, outname)
        write_fasta(fasta_path, pdb_id, final_seqs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--output_dir",
        type=str,
        default=".",
        help="Output dirrectory",
    )
    parser.add_argument(
        "--designed_sequences_fasta",
        type=str,
        required=True,
        help="Path to a fasta file containing designed sequences for the chain of interest.",
    )
    parser.add_argument(
        "--orig_pdb_dir", 
        type=str, 
        default="./subset_cleaned_pdb_files",
        help="Directory containing the original PDB file."
    )
    parser.add_argument(
        "--ensemble_pdb_dir",
        type=str,
        default="./subset_peptide_ensembles",
        help="Directory containing subdirectries of the synthetic ensembles.",
    )
    parser.add_argument(
        "--dataset_info_csv",
        type=str,
        default="protein_peptide_rep_subset_updated_info.csv",
        help="CSV file containing the dataset information, including the PDB IDs and peptide chain of interest for each structure.",
    )
    parser.add_argument(
        "--dataset",
        type=str,
        default="peptide_ensemble",
    )
    args = parser.parse_args()

    # read designed sequences 
    # exmaple: 
    # >1lvm_dup1_1 subset_peptide_ensemble_inference_results
    # KVKELKF
    designed_seqs = {}
    with open(args.designed_sequences_fasta, "r") as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:].strip().split(" ")[0]
            elif current_id is not None:
                
                current_id = "_".join(current_id.split('_')[:-1]) # remove the _# suffix 
                if current_id in designed_seqs:
                    designed_seqs[current_id].append(line.strip())
                else:
                    designed_seqs[current_id] = [line.strip()]
    
    # get the pdb ids and peptide chain of interest 
    pdb_id_to_peptide_info = load_peptide_ensemble_info(args.dataset_info_csv)

    ### parse ###
    for pdb_id in pdb_id_to_peptide_info.keys():
        print(f"\nProcessing {pdb_id}...")
        pep_chain, ref_seq = pdb_id_to_peptide_info[pdb_id]

        
        if args.dataset == "peptide_baseline":
            # load the original structure
            pdb_file_path = os.path.join(args.orig_pdb_dir, f"{pdb_id}.pdb")
            if not os.path.exists(pdb_file_path):
                print(f"Warning! Missing original PDB: {pdb_file_path}")
                continue

        elif args.dataset == "peptide_ensemble":
            # load ensemble pdb
            pdb_dir = os.path.join(args.ensemble_pdb_dir, pdb_id)
            pdb_file_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
            
            if not os.path.exists(pdb_file_path):
                print(f"Warning! Missing topology PDB: {pdb_file_path}")
                continue

        if designed_seqs.get(pdb_id) is None:
            print(f"Warning! No designed sequence found for {pdb_id} in the designed sequences fasta file.")
            continue

        pdb_file_name, struct = extract_structure_from_pdb(pdb_file_path)
        pdb_output_dir = os.path.join(args.output_dir, pdb_id)
        os.makedirs(pdb_output_dir, exist_ok=True)
        write_full_fasta_sequence(pdb_id, pep_chain, struct, designed_seqs[pdb_id], pdb_output_dir)


