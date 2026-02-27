# Curate Dataset

This directory contains scripts for curating a high-quality protein-peptide complex dataset from the Protein Data Bank (PDB).

## Overview

The workflow identifies protein-peptide complexes with significant binding interfaces, downloads their structures, and prepares them for downstream analysis. The dataset is filtered based on structural quality, peptide length, and buried surface area at the interface.

## Files

- **`curate_dataset.py`**: Main script for identifying and downloading protein-peptide complexes
- **`convert_cif_to_pdbs.py`**: Converts mmCIF format files to PDB format with proper chain ID handling
- **`utils.py`**: Shared utility functions (e.g., reading PDB IDs, extracting sequences)
- **`get_rep_subset/`**: Scripts for generating a non-redundant representative subset (see subdirectory README)

## Workflow

### 1. Curate the Dataset

```bash
python curate_dataset.py
```

This script:
1. **Queries RCSB PDB** for protein-peptide complexes with:
   - Resolution ≤ 2.5 Å
   - At least one peptide chain (5-25 residues with ≥5 standard amino acids)
   - At least one protein chain (>30 residues)

2. **Downloads structures** as biological assembly 1 (mmCIF format) to `dataset_cif_files/`

3. **Calculates buried surface area (BSA)** using FreeSASA for each protein-peptide pair

4. **Filters complexes** by BSA ≥ 400 Å² threshold

5. **Saves results** to:
   - `protein-peptide_complexes_candidates.csv`: Initial candidate PDB IDs
   - `protein_peptide_dataset.csv`: Final filtered dataset with columns:
     - `pdb_id`: PDB identifier
     - `Peptide`: Chain ID of the peptide
     - `Buried surface area`: BSA in Ų
     - `Sequence`: Peptide amino acid sequence

### 2. Generate Representative Subset (Optional)

See [get_rep_subset/README.md](get_rep_subset/README.md) for clustering proteins by sequence similarity to create a non-redundant subset.

## Output Files

- `protein_peptide_dataset.csv`: Complete curated dataset
- `dataset_cif_files/`: Downloaded mmCIF structure files
- `dataset_pdb_files/`: PDB format structures (if converted)
- `protein_peptide_rep_subset.csv`: Non-redundant subset at 30% sequence identity (if generated)

### 3. Convert to PDB Format (Optional)

```bash
python convert_cif_to_pdbs.py \
    --csv protein_peptide_dataset.csv \
    --cif-dir dataset_cif_files \
    --pdb-dir dataset_pdb_files
```

Converts mmCIF files to PDB format with automatic chain ID remapping when needed (e.g., multi-character chain IDs). Structures with too many chains can be skipped to avoid conversion issues.

**Options:**
- `--engine {auto,gemmi,biopython}`: Conversion method (default: auto, uses gemmi if available)
- `--workers N`: Parallel workers for gemmi conversions (default: 1)
- `--max-chains N`: Skip structures with more than N chains (default: 5)
- `--no-gemmi-shorten`: Disable gemmi's chain ID shortening

If only want to convert a small subset to PDB format (use for early testing purposes):

```bash
python convert_cif_to_pdbs.py \
    --csv protein_peptide_rep_subset.csv \
    --cif-dir dataset_cif_files \
    --pdb-dir subset_pdb_files
    --max-chains 5
```

## Requirements

- Python packages: `pandas`, `requests`, `gemmi`, `freesasa`, `biopython`
- Optional: `gemmi` CLI tool for faster CIF→PDB conversion
- Optional: `mmseqs2` for sequence clustering (representative subset only)
