"""
Compute rmsd of between structure and reference pdb. Parse results (pLDDT, pTM).

"""

import mdtraj as md
from pathlib import Path
import numpy as np
import csv
import re
import pandas as pd
import argparse
from typing import Optional
import os
import tempfile
from Bio.PDB import MMCIFParser, PDBIO

VERBOSE = False

def load_metadata(csv_path):
    """Load and preprocess multiconf metadata for alignment indices."""

    df = pd.read_csv(csv_path)
    metadata = {}

    for _, row in df.iterrows():
        key = row["pdb_id"].strip()
        metadata[key] = {
            "alignment_idx": [int(x) for x in str(row["alignment_idx"]).split(",") if x != "" and x != "nan"],
        }

    return metadata


def find_model_pdbs(out_dir: Path) -> Path:
    """Find the 'rank_001' model pdb files under a ColabFold output directory."""
    pats = list(out_dir.glob("**/*.pdb"))
    if not pats:
        return None

    # Return only best-ranked model: filename contains 'rank_001'
    rank1 = [p for p in pats if "rank_001" in p.name.lower()]
    if rank1:
        return rank1[0]

    # If none found, return empty so caller can skip (strict behavior requested)
    return None


def process_out_root(metadata: dict, out_root: Path, summary_csv: Path, is_localunfolding=False) -> int:
    """Go through alphafold output root and compute RMSDs, parse pLDDT info, writing summary CSV."""
   
    rows = []

    # find all *af_multimer_results directories
    for base_out in sorted(out_root.rglob("*af_multimer_results")):
        if not any(f.endswith(".done.txt") for f in os.listdir(base_out)):
            print(f"Skipping {base_out} because no .done.txt file found")
            continue
        # base_out is e.g., out_root/pdb_id/pdb_id_#/af_multimer_results
        row_entry_name = base_out.parent.name

        # extract rank_001 metrics from log.txt if present
        plddt_val = None
        ptm_val = None
        iptm_val = None
        log_path = base_out / "log.txt"
        if log_path.exists():
            try:
                txt = log_path.read_text()
                m = re.search(r"rank_001[^\n]*pLDDT=([0-9.]+)[^\n]*pTM=([0-9.]+)[^\n]*ipTM=([0-9.]+)", txt)
                if m:
                    plddt_val = float(m.group(1))
                    ptm_val = float(m.group(2))
                    iptm_val = float(m.group(3))
            except Exception:
                pass

        
        # calculate peptide plddt
        pep_plddt_val = None
        # get the plddt values for the CA atoms of the peptide chain(s) and average them
        model_pdb = find_model_pdbs(base_out)
        if model_pdb is not None:
            try:
                traj = md.load(model_pdb)
                # get CA atoms of peptide chain (chain A)
                peptide_chains = [c for c in traj.topology.chains if c.name == "A"]
                if peptide_chains:
                    ca_atoms = [a for c in peptide_chains for a in c.atoms if a.name == "CA"]
                    if ca_atoms:
                        pep_plddt_vals = [a.bfactor for a in ca_atoms]
                        pep_plddt_val = np.mean(pep_plddt_vals)
            except Exception:
                pass


        rows.append({
            "id": row_entry_name,
            "plddt": plddt_val,
            "ptm": ptm_val,
            "iptm": iptm_val,
            "peptide plddt": pep_plddt_val
        })
    

    # write CSV
    fieldnames = ["id", "plddt", "ptm", "iptm", "peptide plddt"]
    with summary_csv.open("w", newline="") as wf:
        w = csv.DictWriter(wf, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote summary to {summary_csv} ({len(rows)} rows)")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute RMSD and pLDDT metrics from ColabFold af-multimer outputs")
    parser.add_argument("--out_root", 
        required=True,
        help="Root directory containing <seq_id>_results directories (or dataset/<seq_id>_results)",
    )
    parser.add_argument(
        "--dataset_name",
        type=str,
        default="peptide_ensemble",
        help="Name of the dataset to process (default: peptide_ensemble)",
    )
    parser.add_argument(
        "--processed_data_location",
        type=str,
        default="./subset_processed_data/",
        help="Directory storing processed data",
    )
    parser.add_argument("--output_csv_dir", 
        type=str,
        default="./",
        help="Path to output CSV metrics summary",
    )
    parser.add_argument("--csv-name",
        type=str,
        default="af-multimer_results_summary.csv",
        help="Name of the output CSV file",
    )

    args = parser.parse_args()

    out_root = Path(args.out_root).expanduser().resolve()
    if not out_root.exists():
        print(f"ERROR: out_root not found: {out_root}")
        return 2

    if args.output_csv_dir:
        summary = Path(args.output_csv_dir).expanduser().resolve() / args.csv_name
    else:
        summary = out_root / args.csv_name

    dataset_mapping = {
        "peptide_ensemble": "synthetic_ensemble",
        "peptide_baseline": "orig_structures",
    }

    ### load benchmark metadata ###
    metadata_path = os.path.join(
        args.processed_data_location,
        dataset_mapping[args.dataset_name],
        "multiconf_metadata.csv",
    )

    metadata = load_metadata(metadata_path)
    return process_out_root(metadata, out_root, summary)

if __name__ == "__main__":
    raise SystemExit(main())
