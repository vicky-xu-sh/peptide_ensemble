import argparse
import csv
import os
import re
from typing import Dict


def replace_vars(template_text: str, values: Dict[str, str]) -> str:
    def repl(match: re.Match) -> str:
        key = match.group(1)
        if key in values:
            return f"{key}={values[key]}"
        return match.group(0)

    pattern = re.compile(r"^(OUTPUT_PREFIX|INPUT_PDB|NUM_DESIGNS|PD_STEPS|CONTIG_STR)=.*$", re.MULTILINE)
    return pattern.sub(repl, template_text)


def main():
    parser = argparse.ArgumentParser(
        description="Generate partial_diff sh scripts from a CSV and a template."
    )
    parser.add_argument(
        "--csv", 
        default="protein_peptide_rep_subset_updated_info.csv", 
        help="Path to CSV with pdb_id and contig string columns"
    )
    parser.add_argument(
        "--template", 
        default="template_complex_partial_diff.sh", 
        help="Path to template sh script"
    )
    parser.add_argument(
        "--pdb_dir", 
        default="../subset_cleaned_pdb_files",
        help="Root input PDB directory"
    )
    parser.add_argument(
        "--scripts_dir",
        default="subset_partial_diff_scripts",
        help="Directory to write generated sh scripts",
    )
    parser.add_argument("--num_designs", type=int, default=10, help="Value for NUM_DESIGNS (default: 10)")
    parser.add_argument("--pd_steps", type=int, default=10, help="Value for PD_STEPS (default: 10)")
    parser.add_argument("--output_root", required=True, help="Root output directory for RFdiffusion runs")

    args = parser.parse_args()

    with open(args.template, "r", encoding="utf-8") as f:
        template_text = f.read()

    os.makedirs(args.scripts_dir, exist_ok=True)

    with open(args.csv, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        required_cols = {"pdb_id", "fixed_chain_res_str"}
        missing = required_cols - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV is missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            pdb_id = (row.get("pdb_id") or "").strip().lower()
            contig_str = (row.get("fixed_chain_res_str") or "").strip()
            if not pdb_id or not contig_str:
                continue

            output_dir = os.path.join(os.path.abspath(args.output_root), pdb_id)
            output_prefix = output_dir + f"/{pdb_id}_design"
            input_pdb = os.path.join(os.path.abspath(args.pdb_dir), f"{pdb_id}.pdb")

            values = {
                "OUTPUT_PREFIX": output_prefix,
                "INPUT_PDB": input_pdb,
                "NUM_DESIGNS": str(args.num_designs),
                "PD_STEPS": str(args.pd_steps),
                "CONTIG_STR": "'"+ contig_str + "'",
            }

            script_text = replace_vars(template_text, values)
            script_path = os.path.join(args.scripts_dir, f"{pdb_id}_partial_diff.sh")

            with open(script_path, "w", encoding="utf-8") as out_f:
                out_f.write(script_text)

            os.chmod(script_path, 0o755)


if __name__ == "__main__":
    main()
