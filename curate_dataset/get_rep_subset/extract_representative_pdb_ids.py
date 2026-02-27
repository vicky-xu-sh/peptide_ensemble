"""Extract representative PDB IDs from MMseqs2 clustering results.

Reads the cluster TSV file where the first column contains centroid/representative
PDB IDs and saves the unique representative IDs to a CSV file.
"""

import argparse
import csv
from pathlib import Path


def extract_pdb_id(full_id: str) -> str:
    """Extract the PDB ID from the full identifier.
    
    Args:
        full_id: Full identifier like "1EG4|chain=A|peptide=P"
    
    Returns:
        Just the PDB ID, e.g., "1EG4"
    """
    return full_id.split("|")[0]


def main():
    parser = argparse.ArgumentParser(
        description="Extract representative PDB IDs from MMseqs2 clustering"
    )
    parser.add_argument(
        "--tsv",
        default="cluster/cluster_receptor_seqs_cluster.tsv",
        help="Path to MMseqs2 cluster TSV file",
    )
    parser.add_argument(
        "--output",
        default="protein_peptide_rep_subset.csv",
        help="Output CSV file for representative PDB IDs",
    )
    args = parser.parse_args()

    # Read TSV and collect unique representative IDs from first column
    representatives = set()
    
    with open(args.tsv, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                centroid = parts[0]
                pdb_id = extract_pdb_id(centroid)
                representatives.add(pdb_id)
    
    # Sort for consistent output
    representatives = sorted(representatives)
    
    # Write to CSV
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdb_id"])
        for pdb_id in representatives:
            writer.writerow([pdb_id])
    
    print(f"Extracted {len(representatives)} representative PDB IDs")
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
