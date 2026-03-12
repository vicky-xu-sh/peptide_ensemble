# Evaluate with AlphaFold-Multimer

Scripts for preparing AlphaFold-Multimer inputs, running structure predictions, parsing results, and comparing metrics between the baseline and ensemble-based sequence generation approaches.

---

## Scripts

### `prepare_af_multimer_input.py`

Prepares FASTA input files for ColabFold from designed peptide sequences. For each complex, it takes the receptor chain sequences from the original or ensemble PDB structure and substitutes the designed peptide sequence in place of the original peptide chain, writing a multi-chain FASTA file ready for ColabFold.

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `--designed_sequences_fasta` | *(required)* | FASTA file of designed peptide sequences from the sequence generation model |
| `--dataset_info_csv` | `protein_peptide_rep_subset_updated_info.csv` | CSV with PDB IDs and peptide chain IDs for each complex |
| `--dataset` | `peptide_ensemble` | `peptide_ensemble` to use ensemble PDBs as the receptor scaffold, or `peptide_baseline` to use the original experimental PDBs |
| `--orig_pdb_dir` | `./subset_cleaned_pdb_files` | Directory containing original PDB files (used when `--dataset peptide_baseline`) |
| `--ensemble_pdb_dir` | `./subset_peptide_ensembles` | Directory containing RFdiffusion ensemble PDB subdirectories (used when `--dataset peptide_ensemble`) |
| `--output_dir` | `.` | Output directory for generated FASTA files |

**Output structure:**
```
output_dir/
└── <pdb_id>/
    └── <pdb_id>_<n>/
        └── <pdb_id>_<n>.fasta   # one FASTA per designed sequence
```

---

### `run_af_multimer.py`

Runs `colabfold_batch` on all FASTA files prepared by `prepare_af_multimer_input.py`. Iterates over the output directory structure, skipping any jobs that have already completed (detected by the presence of a `.done.txt` file). Supports parallel execution via `--max-workers`.

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `--root_dir` | *(required)* | Root directory produced by `prepare_af_multimer_input.py`, containing per-PDB subdirectories with FASTA files |
| `--colabfold-bin` | `~/colabfold_venv/bin/colabfold_batch` | Path to the `colabfold_batch` binary |
| `--model-name` | `alphafold2_multimer_v3` | ColabFold model to use |
| `--max-workers` | `1` | Number of ColabFold jobs to run concurrently |
| `--dry-run` | — | Print commands without executing them |

**Results** are written to `af_multimer_results/` inside each FASTA's subfolder.

---

### `parse_af_multimer_results.py`

Parses ColabFold output directories and compiles a summary CSV of per-complex structural confidence metrics. For each completed run, it extracts pLDDT, pTM, and ipTM from `log.txt`, and computes the peptide chain pLDDT separately by averaging the B-factor (pLDDT) values of Cα atoms in chain `A` of the top-ranked model PDB.

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `--out_root` | *(required)* | Root directory containing `*af_multimer_results` subdirectories |
| `--output_csv_dir` | `./` | Directory to write the summary CSV |
| `--csv-name` | `af-multimer_results_summary.csv` | Output CSV filename |

**Output CSV columns:** `id`, `plddt`, `ptm`, `iptm`, `peptide_plddt`

> **Note:** Only runs with a `.done.txt` file are included. Runs without it are skipped.

---

### `plot_comparison.py`

Compares AlphaFold-Multimer metrics between the baseline (single original structure) and ensemble approaches. Loads the two summary CSVs produced by `parse_af_multimer_results.py` and generates box plots and scatter plots for each metric (`plddt`, `ptm`, `iptm`, `peptide_plddt`). Two comparisons are made: using the **mean** metric across sequences per complex, and using the **best** (max) metric.

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `--csv_dir` | `./` | Directory containing the two summary CSVs |
| `--baseline_csv_name` | `af-multimer_results_summary_baseline.csv` | Filename of the baseline results CSV |
| `--ensemble_csv_name` | `af-multimer_results_summary_ensemble.csv` | Filename of the ensemble results CSV |
| `--output_dir` | `../plots` | Directory to save output plots |

**Output:** Box plots and scatter plots saved as PNGs, named `<metric>_box_plot_<mean|best>.png` and `<metric>_scatter_plot_<mean|best>.png`.

---

## Typical Usage

### 1. Prepare FASTA inputs

For the ensemble condition:
```bash
python prepare_af_multimer_input.py \
    --designed_sequences_fasta designed_sequences.fasta \
    --dataset peptide_ensemble \
    --ensemble_pdb_dir /path/to/subset_peptide_ensembles \
    --dataset_info_csv protein_peptide_rep_subset_updated_info.csv \
    --output_dir ./af_inputs_ensemble
```

For the baseline condition:
```bash
python prepare_af_multimer_input.py \
    --designed_sequences_fasta designed_sequences.fasta \
    --dataset peptide_baseline \
    --orig_pdb_dir /path/to/subset_cleaned_pdb_files \
    --dataset_info_csv protein_peptide_rep_subset_updated_info.csv \
    --output_dir ./af_inputs_baseline
```

### 2. Run ColabFold

```bash
python run_af_multimer.py \
    --root_dir ./af_inputs_ensemble \
    --colabfold-bin ~/colabfold_venv/bin/colabfold_batch
```

### 3. Parse results

```bash
python parse_af_multimer_results.py \
    --out_root ./af_inputs_ensemble \
    --output_csv_dir ./metrics \
    --csv-name af-multimer_results_summary_ensemble.csv

python parse_af_multimer_results.py \
    --out_root ./af_inputs_baseline \
    --output_csv_dir ./metrics \
    --csv-name af-multimer_results_summary_baseline.csv
```

### 4. Plot comparison

```bash
python plot_comparison.py \
    --csv_dir ./metrics \
    --output_dir ./plots
```

---

## Dependencies

- Python 3.11
- [Biopython](https://biopython.org/) — Cα pLDDT extraction from PDB files
- [Gemmi](https://gemmi.readthedocs.io/) — PDB parsing and assembly handling
- [ColabFold](https://github.com/sokrypton/ColabFold) (`colabfold_batch`) — AlphaFold-Multimer predictions
- NumPy, Pandas, Matplotlib, Seaborn