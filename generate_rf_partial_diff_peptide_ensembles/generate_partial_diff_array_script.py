import argparse
import csv
import os
from typing import List, Tuple


def load_rows(csv_path: str) -> List[Tuple[str, str]]:
    rows: List[Tuple[str, str]] = []
    with open(csv_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        required_cols = {"pdb_id", "fixed_chain_res_str"}
        missing = required_cols - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV is missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            pdb_id = (row.get("pdb_id") or "").strip().lower()
            fixed_chain_res_str = (row.get("fixed_chain_res_str") or "").strip()
            if not pdb_id or not fixed_chain_res_str:
                continue
            rows.append((pdb_id, fixed_chain_res_str))

    return rows


def write_array_script(
    path: str,
    csv_path: str,
    pdb_dir: str,
    output_root: str,
    num_designs: int,
    pd_steps: int,
    array_max_index: int,
) -> None:
    script = f"""#!/bin/bash
#SBATCH --job-name=rf_partial_diff
#SBATCH --output=partial_diff_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=01:30:00
#SBATCH --array=0-{array_max_index}

export INPUT_CSV={csv_path}
export PDB_DIR={pdb_dir}
export OUTPUT_ROOT={output_root}
export NUM_DESIGNS={num_designs}
export PD_STEPS={pd_steps}

line=$(python - <<'PY'
import csv
import os

csv_path = os.environ.get("INPUT_CSV")
if not csv_path:
    raise SystemExit("Missing INPUT_CSV")

idx_env = os.environ.get("SLURM_ARRAY_TASK_ID", "0")
idx = int(idx_env)

with open(csv_path, "r", encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    for i, row in enumerate(reader):
        if i == idx:
            pdb_id = (row.get("pdb_id") or "").strip().lower()
            contig = (row.get("fixed_chain_res_str") or "").strip()
            if pdb_id and contig:
                print(pdb_id + "\t" + contig)
            break
PY
)

if [ -z "$line" ]; then
    echo "No valid row for index $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

pdb_id=$(echo "$line" | cut -f1)
contig_str=$(echo "$line" | cut -f2)

OUTPUT_DIR="${{OUTPUT_ROOT}}/${{pdb_id}}"
OUTPUT_PREFIX="${{OUTPUT_DIR}}/${{pdb_id}}_design"
INPUT_PDB="${{PDB_DIR}}/${{pdb_id}}.pdb"
CONTIG_STR="${{contig_str}}"

echo "Running partial diffusion for PDB ID: $pdb_id with contig: $CONTIG_STR"
echo "Output will be saved to: $OUTPUT_PREFIX"
echo "Input PDB file: $INPUT_PDB"
echo "Number of designs: $NUM_DESIGNS. Partial diffusion steps: $PD_STEPS"

# Load Apptainer
module load apptainer

# Run the container
apptainer run --nv -C -W $SLURM_TMPDIR -B $HOME/scratch \\
    $HOME/scratch/RFdiffusion/rfd.sif \\
    inference.output_prefix=$OUTPUT_PREFIX \\
    inference.model_directory_path=/app/RFdiffusion/models \\
    inference.schedule_directory_path=schedules \\
    inference.input_pdb=$INPUT_PDB \\
    inference.num_designs=$NUM_DESIGNS \\
    contigmap.contigs="$CONTIG_STR" \\
    diffuser.partial_T=$PD_STEPS
"""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(script)
    os.chmod(path, 0o755)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a Slurm job array script and input list from a CSV."
    )
    parser.add_argument(
        "--csv",
        default="protein_peptide_rep_subset_updated_info.csv",
        help="Path to CSV with pdb_id and contig string columns",
    )
    parser.add_argument(
        "--pdb_dir",
        default="../subset_cleaned_pdb_files",
        help="Root input PDB directory",
    )
    parser.add_argument(
        "--output_root",
        required=True,
        help="Root output directory for RFdiffusion runs",
    )
    parser.add_argument(
        "--array_script",
        default="subset_partial_diff_array.sh",
        help="Output Slurm array sbatch script",
    )
    parser.add_argument(
        "--num_designs",
        type=int,
        default=10,
        help="Value for NUM_DESIGNS (default: 10)",
    )
    parser.add_argument(
        "--pd_steps",
        type=int,
        default=10,
        help="Value for PD_STEPS (default: 10)",
    )

    args = parser.parse_args()
    rows = load_rows(args.csv)
    if not rows:
        raise ValueError("No valid rows found in CSV.")

    write_array_script(
        args.array_script,
        os.path.abspath(args.csv),
        os.path.abspath(args.pdb_dir),
        os.path.abspath(args.output_root),
        args.num_designs,
        args.pd_steps,
        len(rows) - 1,
    )


if __name__ == "__main__":
    main()
