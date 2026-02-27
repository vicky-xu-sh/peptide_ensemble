"""Convert mmCIF files to PDB files for PDB IDs listed in a CSV.

- Reads a CSV with a `pdb_id` column (case-insensitive).
- Expects CIF files in a directory named `dataset_cif_files` with names like
  "<pdb_id>-assembly1.cif" (lowercase PDB ID).
- Writes PDB files into a `pdbs` directory.
- Remaps chain IDs that are invalid for PDB format (must be single char).
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import string
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Set, Tuple
from Bio.PDB import MMCIFParser, PDBIO

ALLOWED_CHAIN_IDS = string.ascii_uppercase + string.ascii_lowercase + string.digits


def read_pdb_ids(csv_path: str) -> List[str]:
    with open(csv_path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {csv_path}")
        field_map = {name.lower(): name for name in reader.fieldnames}
        if "pdb_id" not in field_map:
            raise ValueError(
                f"CSV must contain a 'pdb_id' column. Found: {reader.fieldnames}"
            )
        col = field_map["pdb_id"]
        ids = [row[col].strip() for row in reader if row.get(col)]
    return [pid for pid in ids if pid]


def _get_next_available(used_ids: Set[str]) -> Optional[str]:
    for c in ALLOWED_CHAIN_IDS:
        if c not in used_ids:
            return c
    return None


def remap_chain_ids_biopython(structure) -> Dict[str, str]:
    """Remap chain IDs for Biopython structures.

    Only remaps when needed (invalid or duplicate IDs within a model).
    Returns a mapping of original -> new IDs (only for changed IDs).
    """
    remap: Dict[str, str] = {}

    for model in structure:
        used: Set[str] = set()
        to_remap = []

        for chain in model:
            original = chain.id
            if len(original) == 1 and original in ALLOWED_CHAIN_IDS and original not in used:
                used.add(original)
            else:
                to_remap.append(chain)

        for chain in to_remap:
            original = chain.id
            new_id = _get_next_available(used)
            if new_id is None:
                raise ValueError("No available PDB chain IDs left for remapping.")
            chain.id = new_id
            used.add(new_id)
            remap[original] = new_id

    return remap


def count_chains_in_cif(cif_path: str) -> int:
    parser = MMCIFParser(QUIET=True)
    structure_id = os.path.basename(cif_path).split(".")[0]
    structure = parser.get_structure(structure_id, cif_path)
    first_model = next(structure.get_models())
    return sum(1 for _ in first_model.get_chains())


def convert_one(
    cif_path: str, pdb_path: str) -> Dict[str, str]:
    parser = MMCIFParser(QUIET=True)
    structure_id = os.path.basename(cif_path).split(".")[0]
    structure = parser.get_structure(structure_id, cif_path)
    remap = remap_chain_ids_biopython(structure)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)
    return remap


def _convert_with_gemmi(
    gemmi_bin: str, cif_path: str, pdb_path: str, gemmi_shorten: bool
) -> None:
    cmd = [gemmi_bin, "convert", cif_path, pdb_path]
    if gemmi_shorten:
        cmd.insert(2, "--shorten")
    subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert CIFs to PDBs from a CSV list.")
    parser.add_argument(
        "--csv",
        default="protein_peptide_dataset.csv",
        help="CSV file with pdb_id column",
    )
    parser.add_argument(
        "--cif-dir",
        default="dataset_cif_files",
        help="Directory with CIF files",
    )
    parser.add_argument(
        "--pdb-dir",
        default="dataset_pdb_files",
        help="Output directory for PDB files",
    )
    parser.add_argument(
        "--engine",
        choices=["auto", "gemmi", "biopython"],
        default="auto",
        help="Conversion engine (auto uses gemmi CLI if available)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel workers for gemmi conversions (>=1)",
    )
    parser.add_argument(
        "--max-chains",
        type=int,
        default=5,
        help="Skip conversion when CIF has more than this many chains",
    )
    parser.add_argument(
        "--no-gemmi-shorten",
        action="store_true",
        help="Disable gemmi --shorten when converting to PDB",
    )
    args = parser.parse_args()

    pdb_ids = read_pdb_ids(args.csv)
    if not pdb_ids:
        raise ValueError("No PDB IDs found in CSV.")

    os.makedirs(args.pdb_dir, exist_ok=True)

    missing: List[str] = []
    failed: List[Tuple[str, str]] = []
    skipped: List[Tuple[str, int]] = []

    jobs: List[Tuple[str, str, str]] = []

    for pdb_id in pdb_ids:
        pdb_id_lower = pdb_id.lower()
        cif_name = f"{pdb_id_lower}-assembly1.cif"
        cif_path = os.path.join(args.cif_dir, cif_name)
        pdb_path = os.path.join(args.pdb_dir, f"{pdb_id_lower}.pdb")

        if not os.path.exists(cif_path):
            missing.append(pdb_id)
            continue

        try:
            chain_count = count_chains_in_cif(cif_path)
        except Exception as exc:
            failed.append((pdb_id, f"Failed to read CIF for chain count: {exc}"))
            continue

        if chain_count > args.max_chains:
            skipped.append((pdb_id, chain_count))
            continue

        jobs.append((pdb_id, cif_path, pdb_path))

    gemmi_bin = None
    if args.engine in {"auto", "gemmi"}:
        gemmi_bin = shutil.which("gemmi")
        if not gemmi_bin and args.engine == "gemmi":
            raise RuntimeError("gemmi executable not found in PATH.")

    gemmi_shorten = not args.no_gemmi_shorten

    if gemmi_bin and args.engine in {"auto", "gemmi"}:
        workers = max(1, args.workers)
        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_map = {
                executor.submit(
                    _convert_with_gemmi, gemmi_bin, cif_path, pdb_path, gemmi_shorten
                ): pdb_id
                for pdb_id, cif_path, pdb_path in jobs
            }
            for future in as_completed(future_map):
                pdb_id = future_map[future]
                try:
                    future.result()
                except Exception as exc:
                    failed.append((pdb_id, str(exc)))
    else:
        for pdb_id, cif_path, pdb_path in jobs:
            try:
                remap = convert_one(cif_path, pdb_path)
                if remap:
                    print(f"{pdb_id}: remapped chains {remap}")
            except Exception as exc:
                failed.append((pdb_id, str(exc)))

    if missing:
        print(f"Missing CIFs ({len(missing)}): {', '.join(missing)}")
    if skipped:
        print("Skipped due to chain count:")
        for pdb_id, chain_count in skipped:
            print(f"  {pdb_id}: {chain_count} chains")
    if failed:
        print("Failed conversions:")
        for pdb_id, err in failed:
            print(f"  {pdb_id}: {err}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


