# Generate RFdiffusion Partial Diffusion Inputs

This directory contains scripts to prepare PDB structures and generate SLURM job scripts for running RFdiffusion partial diffusion on protein-peptide complexes.

## Overview

The workflow processes protein-peptide complex structures to:
1. Clean and reorder chains (peptide chain first)
2. Renumber residues to non-negative indices
3. Generate contig strings specifying which regions to keep fixed during diffusion
4. Create SLURM array scripts for high-throughput partial diffusion runs

## Files

- **`parse_pdbs.py`**: Processes PDB files, cleans them, and generates contig strings
- **`generate_partial_diff_array_script.py`**: Generates SLURM array job scripts for batch processing
- **`generate_partial_diff_individual_scripts.py`**: Generates individual job scripts per PDB
- **`utils.py`**: Utility functions for PDB cleaning and chain reordering
- **`template_complex_partial_diff.sh`**: Template for individual SLURM job scripts

## Workflow

### 1. Parse and Clean PDB Files

```bash
python parse_pdbs.py \
    --subset-csv ../protein_peptide_rep_subset.csv \
    --peptide-csv ../protein_peptide_dataset.csv \
    --pdb-dir <dir_with_input_pdb_files>/ \
    --cleaned-pdb-dir <dir_for_output_cleaned_pdb_files>/ \
    --out protein_peptide_rep_subset_updated_info.csv
```

This script:
- Reads PDB IDs from the subset CSV
- Looks up peptide chain information from the complete dataset CSV
- **Cleans each PDB structure**:
  - Reorders chains to place the peptide chain first (becomes chain A)
  - Renumbers all residues to non-negative sequential indices
  - Handles structures with multiple peptide chains (creates `_dup` variants)
- **Generates contig strings** that specify:
  - Which residues to keep fixed during diffusion (protein receptor regions)
  - Which residues to diffuse (peptide, specified as final residue range)
- Saves cleaned PDBs to `cleaned-pdb-dir`
- Outputs updated CSV with columns:
  - `pdb_id`: PDB identifier (with `_dup` suffix for multiple peptides)
  - `total_residues`: Total residue count
  - `peptide_chain_id`: New peptide chain ID (always 'A')
  - `peptide_sequence`: Peptide amino acid sequence
  - `chain_res_str`: Chain residue breakdown
  - `fixed_chain_res_str`: Contig string for RFdiffusion

**Example contig string**: `[5-5/0 B16-284]`
- Fix receptor chains B (residues 16-284)
- Diffuse chain A residues (the peptide chain)

### 2. Generate SLURM Scripts

You have two options for running RFdiffusion:

#### Option A: Array Script (Recommended for batch processing)

```bash
python generate_partial_diff_array_script.py \
    --csv protein_peptide_rep_subset_updated_info.csv \
    --pdb_dir <dir_with_cleaned_pdb_files> \
    --output_root <output_root_dir> \
    --array_script subset_partial_diff_array.sh \
    --num_designs 10 \
    --pd_steps 10
```

Creates a single SLURM job array script that runs RFdiffusion partial diffusion for all structures in the CSV. This is the most efficient approach for processing many structures.

**Options:**
- `--csv`: Input CSV with PDB IDs and contig strings
- `--pdb_dir`: Directory containing cleaned PDB files
- `--output_root`: Root directory for RFdiffusion outputs
- `--num_designs`: Number of designs to generate per structure (default: 10)
- `--pd_steps`: Number of partial diffusion timesteps (default: 10)

**Submit:**
```bash
sbatch subset_partial_diff_array.sh
```

#### Option B: Individual Scripts (For selective processing)

```bash
python generate_partial_diff_individual_scripts.py \
    --csv protein_peptide_rep_subset_updated_info.csv \
    --template template_complex_partial_diff.sh \
    --pdb_dir <dir_with_cleaned_pdb_files> \
    --output_root <output_root_dir> \
    --scripts_dir <dir_for_output_generated_scripts> \
    --num_designs 10 \
    --pd_steps 10
```

Generates individual SLURM job scripts for each PDB structure using the provided template. This approach is useful when you need to:
- Run only specific structures
- Customize parameters for individual structures
- Debug or rerun failed jobs

**Options:**
- `--template`: Path to template script (default: `template_complex_partial_diff.sh`)
- `--scripts_dir`: Directory for generated scripts (default: `subset_partial_diff_scripts`)
- Other options same as Option A

**Submit individual jobs:**
```bash
sbatch <dir_with_generated_scripts>/1abc_partial_diff.sh
```

Or submit all at once:
```bash
for script in <dir_with_generated_scripts>/*.sh; do sbatch "$script"; done
```

## Job Execution

Both script generation options will:
- Run one job per PDB structure
- Load the RFdiffusion Apptainer container
- Execute partial diffusion with the specified contig constraints
- Save outputs to `{output_root}/{pdb_id}/`

After the jobs finish, scan SLURM logs for failures:

```bash
sh grep_sbatch_errors.sh <dir_with_sbatch_outputs>
```

This script searches the job output directory for common error patterns and prints a summary of problematic jobs.

### 3. Unmerge Chains from Design Output

After RFdiffusion completes (if using an old RFdiffusion version), the generated PDBs contain only the designed peptide chain merged into the receptor. To restore the original chain IDs and generate trajectory ensembles:

```bash
python unmerge_chains.py \
    --orig_pdb_dir <dir_with_cleaned_pdb_files> \
    --generated_pdbs_dir <output_root_dir> \
    --output_dir <dir_for_output_ensembles>
```

This script:
- **Restores original chain IDs**: RFdiffusion merges all chains into chain A. This script maps residues back to their original chains (A, B, C, etc.) using the cleaned original PDB as reference.
- **Matches residues between original and design**: Uses CA atoms and residue identifiers (resSeq) to track positions across files.
- **Handles non-standard residues**: Skips non-standard amino acids (those not in AA3_TO_1) to ensure accurate residue matching.
- **Generates trajectory ensembles**: Combines all design outputs (design_0, design_1, ..., design_N) for each PDB into:
  - `{pdb_id}.pdb`: First design as single structure
  - `{pdb_id}.xtc`: Full trajectory ensemble (mdtraj format)

**Options:**
- `--orig_pdb_dir`: Directory containing the cleaned original PDB files (input to RFdiffusion)
- `--generated_pdbs_dir`: Root directory containing RFdiffusion outputs (subdirectories by pdb_id)
- `--output_dir`: Directory to save unmerged PDB and XTC files

**Note**: Both the original and design files must have the same number of residues (after filtering non-standard residues) for proper chain restoration.

## Output Structure

```
<output_root_dir>/
  {pdb_id}/
    {pdb_id}_design_0.pdb
    {pdb_id}_design_1.pdb
    ...
```

## Requirements

- Python packages: `biopython`, `pandas`
- RFdiffusion (via Apptainer/Singularity container)
- SLURM workload manager (for HPC clusters)

## Notes

- The peptide chain is always reordered to be chain A in cleaned PDBs
- Structures with multiple peptide chains are duplicated with `_dup` suffixes
- Non-standard residues and negative residue numbers are handled automatically
- The contig string ensures only the peptide region is redesigned while keeping the receptor fixed