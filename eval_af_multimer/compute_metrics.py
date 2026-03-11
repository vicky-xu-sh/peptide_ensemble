from pathlib import Path
import numpy as np
import csv
import re
import argparse
import os
from Bio.PDB import PDBParser


def find_model_pdbs(out_dir: Path) -> Path:
    """Find the 'rank_001' model pdb files under a ColabFold output directory."""
    pats = list(out_dir.glob("*.pdb"))
    if not pats:
        print(f"WARNING: no pdb files in {out_dir}")
        return None

    # Return only best-ranked model: filename contains 'rank_001'
    rank1 = [p for p in pats if "rank_001" in p.name.lower()]
    if rank1:
        return rank1[0]

    # If none found, return empty so caller can skip (strict behavior requested)
    print(f"WARNING: No rank_001 model pdb found under {out_dir}, skipping")
    return None


def process_out_root(out_root: Path, summary_csv: Path, is_localunfolding=False) -> int:
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
                structure = PDBParser(QUIET=True).get_structure("structure", model_pdb)
                model = next(structure.get_models())
                for chain in model:
                    if chain.id == "A":  # assuming peptide chain is labeled 'A'
                        ca_atoms = [res["CA"] for res in chain if "CA" in res]
                        if ca_atoms:
                            pep_plddt_vals = [a.bfactor for a in ca_atoms]
                            pep_plddt_val = np.mean(pep_plddt_vals)
                        else:
                            print(f"No CA atoms found in chain A of {model_pdb}, skipping peptide pLDDT calculation")
                        break
                    else:
                        print(f"No CA atoms found in chain A of {model_pdb}, skipping peptide pLDDT calculation")
            except Exception:
                pass


        rows.append({
            "id": row_entry_name,
            "plddt": plddt_val,
            "ptm": ptm_val,
            "iptm": iptm_val,
            "peptide_plddt": pep_plddt_val
        })
    

    # write CSV
    fieldnames = ["id", "plddt", "ptm", "iptm", "peptide_plddt"]
    with summary_csv.open("w", newline="") as wf:
        w = csv.DictWriter(wf, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote summary to {summary_csv} ({len(rows)} rows)")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Parse and compute pLDDT metrics from ColabFold af-multimer outputs")
    parser.add_argument("--out_root", 
        required=True,
        help="Root directory containing <seq_id>_results directories (or dataset/<seq_id>_results)",
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

    return process_out_root(out_root, summary)

if __name__ == "__main__":
    raise SystemExit(main())
