import argparse
import csv
import os
from typing import Dict, List, Optional, Tuple
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
import pandas as pd
from utils import clean_pdb


def num_residues_in_chain(chain) -> int:
    seen_indices = set()
    
    for res in chain.get_residues():
        # RFDiffusion logic: 
        # 1. Must be ATOM (handled by res.id[0] == " ")
        # 2. Must have a CA
        # 3. Must have a unique (Chain, ResNum)
        if res.id[0] == " " and "CA" in res:
            res_num = res.id[1]
            seen_indices.add(res_num)
            
    return seen_indices

def read_csv_column(path: str, column: str) -> List[str]:
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        field_map = {name.lower(): name for name in reader.fieldnames}
        if column.lower() not in field_map:
            raise ValueError(
                f"CSV must contain a '{column}' column. Found: {reader.fieldnames}"
            )
        col = field_map[column.lower()]
        return [row[col].strip() for row in reader if row.get(col)]


def load_peptide_info(path: str) -> Dict[str, List[str]]:
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        field_map = {name.lower(): name for name in reader.fieldnames}
        if ("pdb_id" not in field_map) or ("peptide" not in field_map) or ("sequence" not in field_map):
            raise ValueError(
                "CSV must contain 'pdb_id', 'peptide', and 'sequence' columns. "
                f"Found: {reader.fieldnames}"
            )
        pdb_col = field_map["pdb_id"]
        pep_chain_col = field_map["peptide"]
        seq_col = field_map["sequence"]
        result: Dict[str, List[Tuple[str, str]]] = {}
        for row in reader:
            pdb_id = row.get(pdb_col, "").strip().upper()
            seq = row.get(seq_col, "").strip().replace(" ", "").replace("\n", "")
            pep_chain = row.get(pep_chain_col, "").strip()
            if pdb_id in result:
                result[pdb_id].append((pep_chain, seq))
            else:
                result[pdb_id] = [(pep_chain, seq)]
        return result

def process_one_pdb(pdb_path: str, peptide_chain_id: str) -> Tuple[int, str]:
    structure = PDBParser(QUIET=True).get_structure("structure", pdb_path)
    model = next(structure.get_models())

    total_residues = 0
    chain_res_counts = []
    fix_chain_strs = []

    for chain in model:
        res_indices = num_residues_in_chain(chain)
        if not res_indices:
            continue

        res_nums = len(res_indices)
        sorted_res = sorted(res_indices)

        if chain.id.strip() != peptide_chain_id:
            # print(f"Processing chain {chain.id.strip()} with residue numbers: {sorted_res}")
            chain_str = ""

            # iterate the indices in order, find contiguous ranges
            i = 0
            while i < res_nums:
                res_num_start = sorted_res[i]
                res_num_end = res_num_start
                # extend the end as long as next index is contiguous
                while (res_num_end + 1) in res_indices:
                    res_num_end += 1
                    i += 1
                # print(f"Found contiguous range: {res_num_start} to {res_num_end}")

                # move to next index
                i += 1

                chain_str += f"{chain.id.strip()}{res_num_start}-{res_num_end}"

                if res_num_end != sorted_res[-1]:
                    chain_str += "/"

                # # Hacky way to prevent chain reordering in RFdiffusion (not used)
                # if res_num_end == sorted_res[-1]:
                #     # if we reached the end of the list,
                #     # have the last residue to be diffused (hack to prevent chain reording)
                #     res_num_end = res_num_end-1
                #     chain_str += f"{chain.id.strip()}{res_num_start}-{res_num_end}/1-1" # 1-1 at the end to have the last residue diffusable
                #     break
                # else:
                #     chain_str += f"{chain.id.strip()}{res_num_start}-{res_num_end}/"

            fix_chain_strs.append(chain_str)

        else:
            # peptide chain
            fix_chain_strs.append(f"{res_nums}-{res_nums}")
        
        total_residues += res_nums
        chain_res_counts.append(res_nums)

    chain_res_str = "[" + "/0 ".join(f"{c}-{c}" for c in chain_res_counts) + "]"
    fixed_chain_res_str = "[" + "/0 ".join(fix_chain_strs) + "]"

    return total_residues, chain_res_str, fixed_chain_res_str


def main():
    parser = argparse.ArgumentParser(
        description="Add total residues and peptide ranges to subset CSV."
    )
    parser.add_argument(
        "--subset-csv",
        default="protein_peptide_rep_subset.csv",
        help="CSV containing pdb_id column for the subset, which will be parsed and updated with new info",
    )
    parser.add_argument(
        "--peptide-csv",
        default="protein_peptide_dataset.csv",
        help="CSV containing PDB/Sequence columns for the complete dataset, used to look up peptide sequences for the subset",
    )
    parser.add_argument(
        "--pdb-dir",
        default="subset_pdb_files",
        help="Directory with the PDB files",
    )
    parser.add_argument(
        "--cleaned-pdb-dir",
        default="subset_cleaned_pdb_files",
        help="Directory for output cleaned PDB files",
    )
    parser.add_argument(
        "--out",
        default="protein_peptide_rep_subset_updated_info.csv",
        help="Output CSV path",
    )
    args = parser.parse_args()

    print(f"Reading peptide info from {args.peptide_csv}...")
    print(f"Reading set of PDB IDs from {args.subset_csv}...")
    pdb_ids = read_csv_column(args.subset_csv, "pdb_id")
    peptide_info = load_peptide_info(args.peptide_csv)

    os.makedirs(args.cleaned_pdb_dir, exist_ok=True)

    results = []

    for pdb_id in pdb_ids:
        pdb_path = os.path.join(args.pdb_dir, f"{pdb_id.lower()}.pdb")
        if not os.path.exists(pdb_path):
            print(f"\nWarning: PDB file not found for {pdb_id} at {pdb_path}")
            continue

        print(f"\nProcessing {pdb_id}...")

        peptide_seqs = peptide_info.get(pdb_id.upper(), [("", "")])
        if len(peptide_seqs) > 1:
            print(f"Multiple peptide chains found for {pdb_id}, processing all.")

        for i, peptide_seq_tuple in enumerate(peptide_seqs):
            pdb_id = pdb_id.split("_dup")[0]  # remove any existing _dup suffix
            pdb_path = os.path.join(args.pdb_dir, f"{pdb_id.lower()}.pdb") # the original pdb file
            
            if i > 0:
                pdb_id = f"{pdb_id}_dup{i}"
            
            old_peptide_chain_id, peptide_seq = peptide_seq_tuple
            # Gemmi probably shortened the chain ID
            if len(old_peptide_chain_id) > 1:
                old_peptide_chain_id = old_peptide_chain_id[0]

            cleaned_pdb_path = os.path.join(args.cleaned_pdb_dir, f"{pdb_id.lower()}.pdb")

            new_pep_chains_mapping = clean_pdb(pdb_path, cleaned_pdb_path, old_peptide_chain_id)
            if not new_pep_chains_mapping or old_peptide_chain_id not in new_pep_chains_mapping:
                print(f"Error: No new peptide chain IDs found for {pdb_id} after cleaning. Skipping this PDB.")
                continue
            
            # this is the new/cleaned pdb file with reordered chains and non-negative residue numbers
            pdb_path = os.path.join(args.cleaned_pdb_dir, f"{pdb_id.lower()}.pdb")
            print(f"New peptide chains for {pdb_id}: {new_pep_chains_mapping}")

            pep_chain = new_pep_chains_mapping.get(old_peptide_chain_id)
            assert pep_chain == 'A', f"Expected new peptide chain ID to be 'A', but got {pep_chain} for {pdb_id}"
            total_residues, chain_res_str, fixed_chain_res_str = process_one_pdb(pdb_path, pep_chain)
            print(f"{pdb_id}: Total residues={total_residues}, Chain residue string={chain_res_str}, Fixed chain residue string={fixed_chain_res_str}")

            # add to result
            results.append({
                "pdb_id": pdb_id,
                "total_residues": total_residues,
                "peptide_chain_id": pep_chain,
                "peptide_sequence": peptide_seq,
                "chain_res_str": chain_res_str,
                "fixed_chain_res_str": fixed_chain_res_str,
            })

    # save result to CSV
    result_df = pd.DataFrame(results)
    result_df.to_csv(args.out, index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())