"""Extract receptor (non-peptide) protein sequences from mmCIF files.

Example:
  python extract_receptor_sequences.py \
    --csv protein_peptide_dataset.csv \
    --cif-dir dataset_cif_files \
    --output receptor_sequences.fasta
"""

import argparse
import csv
import os
from typing import Dict, List, Set, Tuple
import gemmi
from utils import get_chain_seq


def read_peptide_chains(csv_path: str) -> Dict[str, Set[str]]:
    pdb_to_peptides: Dict[str, Set[str]] = {}
    with open(csv_path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pdb_id = row["pdb_id"].strip().upper()
            peptide_chain = row["Peptide"].strip()
            if not pdb_id or not peptide_chain:
                continue
            pdb_to_peptides.setdefault(pdb_id, set()).add(peptide_chain)
    return pdb_to_peptides

def format_header(
    pdb_id: str,
    peptide_chains: Set[str],
    chain_id: str | None = None,
    chain_ids: List[str] | None = None,
) -> str:
    peptide_str = ",".join(sorted(peptide_chains))
    if chain_ids is not None:
        chains_str = ",".join(chain_ids)
        return f">{pdb_id}|chains={chains_str}|peptide={peptide_str}"
    if chain_id is None:
        raise ValueError("Provide chain_id or chain_ids")
    return f">{pdb_id}|chain={chain_id}|peptide={peptide_str}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract receptor (non-peptide) protein sequences from mmCIF files."
    )
    parser.add_argument(
        "--csv",
        default="protein_peptide_dataset.csv",
        help="Path to protein_peptide_dataset.csv",
    )
    parser.add_argument(
        "--cif-dir",
        default="dataset_cif_files",
        help="Directory with mmCIF files (e.g., 1a1m-assembly1.cif)",
    )
    parser.add_argument(
        "--output",
        default="receptor_sequences.fasta",
        help="Output FASTA path",
    )
    parser.add_argument(
        "--concat",
        action="store_true",
        help="Concatenate all receptor chains into a single FASTA record per PDB",
    )
    parser.add_argument(
        "--delimiter",
        default="",
        help="Delimiter inserted between concatenated chains when using --concat",
    )

    args = parser.parse_args()

    pdb_to_peptides = read_peptide_chains(args.csv)
    if not pdb_to_peptides:
        raise SystemExit(f"No PDB entries found in {args.csv}")

    missing_files = 0
    kept_chains = 0
    fasta_lines: List[str] = []

    for pdb_id in sorted(pdb_to_peptides.keys()):
        peptide_chains = pdb_to_peptides[pdb_id]
        cif_name = f"{pdb_id.lower()}-assembly1.cif"
        cif_path = os.path.join(args.cif_dir, cif_name)
        if not os.path.isfile(cif_path):
            missing_files += 1
            continue

        structure = gemmi.read_structure(cif_path)
        model = structure[0]

        receptor_sequences: List[Tuple[str, str]] = []
        for chain in model:
            chain_id = chain.name
            if chain_id in peptide_chains:
                continue

            seq = get_chain_seq(chain)

            if len(seq) < 30:  # not a protein chain
                continue

            receptor_sequences.append((chain_id, seq))

        if not receptor_sequences:
            continue

        if args.concat:
            chain_ids = [cid for cid, _ in receptor_sequences]
            header = format_header(
                pdb_id,
                peptide_chains,
                chain_ids=chain_ids,
            )
            concat_seq = args.delimiter.join([seq for _, seq in receptor_sequences])
            fasta_lines.append(header)
            fasta_lines.append(concat_seq)
            kept_chains += len(receptor_sequences)
        else:
            for chain_id, seq in receptor_sequences:
                header = format_header(
                    pdb_id,
                    peptide_chains,
                    chain_id=chain_id,
                )
                fasta_lines.append(header)
                fasta_lines.append(seq)
                kept_chains += 1

    with open(args.output, "w") as out_handle:
        if fasta_lines:
            out_handle.write("\n".join(fasta_lines) + "\n")

    print(
        f"Wrote {kept_chains} receptor chain sequences to {args.output}. "
        f"Missing CIF files: {missing_files}."
    )


if __name__ == "__main__":
    main()
