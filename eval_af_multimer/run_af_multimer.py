from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import os
import re
import subprocess
import sys
import shutil
from pathlib import Path
from typing import Tuple


def run_colabfold(
    colabfold_bin: Path,
    fasta_file: Path,
    out_dir: Path,
    model_name: str,
    dry_run: bool,
) -> None:

    # check if results already exist
    if out_dir.exists():
        if any(f.endswith(".done.txt") for f in os.listdir(out_dir)):
            print(f"Skipping {fasta_file.name} because results already exist in {out_dir}")
            return

    os.makedirs(out_dir, exist_ok=True)

    cmd = [
        str(colabfold_bin),
        str(fasta_file),
        str(out_dir),
        "--model-type",
        model_name,
    ]

    print(">>>", " ".join(cmd))
    if dry_run:
        return

    subprocess.run(cmd, check=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root_dir", 
        required=True, 
        help="Input root directory containing sub directories of PDBs, with sub-folders containing FASTA files. \
            FASTA header format: >seq_id. Each folder corresponding to one test sequence, and the FASTA file containing the designed sequence for that test sequence.")
    ap.add_argument(
        "--colabfold-bin",
        default="~/colabfold_venv/bin/colabfold_batch",
        help="Path to colabfold_batch binary",
    )
    ap.add_argument(
        "--model-name",
        default="alphafold2_multimer_v3",
        help="ColabFold model name to use (e.g., alphafold2_multimer_v3, alphafold2_multimer_v2, etc.)",
    )
    ap.add_argument(
        "--max-workers",
        type=int,
        default=1,
        help="Number of ColabFold jobs to run concurrently. Use 1 for sequential execution.",
    )
    ap.add_argument("--dry-run", action="store_true", help="Print commands without running them")
    args = ap.parse_args()

    root_dir = Path(args.root_dir).expanduser().resolve()
    colabfold_bin = Path(args.colabfold_bin).expanduser().resolve()

    if not root_dir.exists():
        print(f"ERROR: input not found: {root_dir}", file=sys.stderr)
        return 2
    if not colabfold_bin.exists():
        print(f"ERROR: colabfold_batch not found: {colabfold_bin}", file=sys.stderr)
        return 2
    if args.max_workers < 1:
        print("ERROR: --max-workers must be >= 1", file=sys.stderr)
        return 2

    jobs = []

    # Iterate over subdirectories in root_dir (subdir is each PDB)
    for subdir in sorted(root_dir.iterdir()):
        if not subdir.is_dir():
            continue
        
        pdb_id = subdir.name
        print(f"\n=== Processing {pdb_id} ===")

        # Find the FASTA file in the subfolders of this subdir
        for folder in sorted(subdir.iterdir()):
            if not folder.is_dir():
                continue
            fasta_candidates = sorted(folder.glob("*.fasta"))
            if not fasta_candidates:
                print(f"WARNING: no FASTA file found in {folder}; skipping", file=sys.stderr)
                continue
            if len(fasta_candidates) > 1:
                print(f"WARNING: multiple FASTA files found in {folder}; using {fasta_candidates[0].name}", file=sys.stderr)

            fasta_file = fasta_candidates[0]
            print(f"Found FASTA file: {fasta_file}")

            out_path = folder / "af_multimer_results"
            jobs.append((pdb_id, fasta_file, out_path))

    if args.max_workers == 1:
        for pdb_id, fasta_file, out_path in jobs:
            try:
                run_colabfold(colabfold_bin, fasta_file, out_path, args.model_name, args.dry_run)
            except subprocess.CalledProcessError as e:
                print(f"ERROR: colabfold failed for {pdb_id} with exit code {e.returncode}", file=sys.stderr)
                return e.returncode
    else:
        print(f"\nRunning {len(jobs)} jobs with up to {args.max_workers} workers")
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            future_to_job = {
                executor.submit(run_colabfold, colabfold_bin, fasta_file, out_path, args.model_name, args.dry_run): (pdb_id, fasta_file)
                for pdb_id, fasta_file, out_path in jobs
            }

            for future in as_completed(future_to_job):
                pdb_id, fasta_file = future_to_job[future]
                try:
                    future.result()
                except subprocess.CalledProcessError as e:
                    print(
                        f"ERROR: colabfold failed for {pdb_id} ({fasta_file.name}) with exit code {e.returncode}",
                        file=sys.stderr,
                    )
                    return e.returncode
                except Exception as e:
                    print(
                        f"ERROR: unexpected failure for {pdb_id} ({fasta_file.name}): {e}",
                        file=sys.stderr,
                    )
                    return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())