# Peptide Ensemble

A pipeline for curating a protein–peptide dataset from the PDB, generating synthetic peptide conformational ensembles via RFdiffusion partial diffusion, and evaluating whether ensemble-based inputs improve peptide sequence generation compared to single static structures.

---

## Overview

This project supports four sequential stages of a computational peptide design and evaluation workflow:

1. **Dataset curation** — Parse CIF/PDB structural files from the Protein Data Bank, filter and preprocess protein–peptide complexes at the chain and complex level, and produce a curated dataset ready for downstream use.
2. **Ensemble generation** — Prepare input files and run RFdiffusion partial diffusion to generate diverse synthetic peptide conformational ensembles for a representative subset of complexes.
3. **Sequence generation** *(external)* — Use a sequence generation model of your choice to generate peptide sequences conditioned on the ensemble structures. For baseline comparison, sequences are also generated from the original single static structures.
4. **Evaluation** — Score generated sequences by running AlphaFold-Multimer predictions and computing structural confidence metrics.

---

## Repository Structure

```
peptide_ensemble/
├── curate_dataset/                          # Stage 1: dataset curation scripts
├── generate_rf_partial_diff_peptide_ensembles/  # Stage 2: RFdiffusion ensemble generation
├── eval_af_multimer/                        # Stage 4: AlphaFold-Multimer evaluation
├── protein_peptide_dataset.csv              # Full curated protein–peptide dataset
├── protein-peptide_complexes_candidates.csv # Filtered complex-level candidates
├── protein_peptide_rep_subset.csv           # Representative subset for ensemble generation
└── .gitignore
```

---

## Workflow

### Stage 1 — Curate Dataset (`curate_dataset/`)

Processes raw PDB/CIF structural files to identify and extract protein–peptide complexes. Applies chain-level filters (e.g., peptide length thresholds, resolution cutoffs) and complex-level filters (e.g., Buried Surface Area (BSA)) to produce a clean dataset.

**Output files:**
- `protein_peptide_dataset.csv` — Full dataset of curated complexes
- `protein-peptide_complexes_candidates.csv` — Subset passing complex-level quality filters
- `protein_peptide_rep_subset.csv` — Representative, deduplicated subset for downstream use

### Stage 2 — Generate RF Partial Diffusion Ensembles (`generate_rf_partial_diff_peptide_ensembles/`)

Prepares contig maps and input PDB files for [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) partial diffusion runs. Partial diffusion noises and re-denoises the peptide chain while keeping the receptor fixed, generating a diverse ensemble of plausible peptide conformations for each complex.

### Stage 3 — Generate Peptide Sequences Using a Sequence Generation Model (outside of this repo)

Peptides are highly flexible molecules, and existing sequence generation models are typically conditioned on a single static structure — which fails to capture this conformational diversity and can limit sequence quality. This stage investigates whether providing a structural ensemble as input leads to better sequence generation.

Use a sequence generation model of your choice (e.g., ProteinMPNN) to generate peptide sequences conditioned on the RFdiffusion ensemble structures from Stage 2. For baseline comparison, also generate sequences using the original experimental PDB structures as input.

### Stage 4 — Evaluate with AlphaFold-Multimer (`eval_af_multimer/`)

Evaluates the generated sequences by running AlphaFold-Multimer predictions (using colabfold_batch) and computing structural confidence metrics (e.g., pTM, ipTM, pLDDT, peptide chain pLDDT).

---

## Dependencies

- Python 3.11
- [Biopython](https://biopython.org/) — PDB/CIF parsing
- [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) — partial diffusion inference
- [AlphaFold-Multimer](https://github.com/google-deepmind/alphafold) — structure prediction and evaluation
- NumPy, Pandas

RFdiffusion and AlphaFold-Multimer require separate installation. Please refer to their respective repositories for setup instructions.

---

## Notes

1. Check the README in each subdirectory for more details on each stage.
2. **Known limitation:** `parse_pdbs.py` in `generate_rf_partial_diff_peptide_ensembles/` currently parses PDB files converted from CIF, rather than the original CIF files directly. Parsing from CIF is preferable as it preserves entity-level chain information, which may be useful for sequence generation models. A fix is available in the `fix_parse_pdbs` branch, but has not yet been tested or validated on the downstream task.
3. **Chain ID note:** The peptide chain IDs in `protein_peptide_dataset.csv` currently use `auth_asym_id`. The `fix_parse_pdbs` branch also updates this so that cleaned PDBs use `label_asym_id` instead — this matters because a single `auth_asym_id` chain can span multiple entity subchains.
4. **RFdiffusion version limitation:** This pipeline was run using an older version of RFdiffusion inside a container on Compute Canada. That version has a known bug where chains are merged together after diffusion if the partially diffused chain is not placed first in the input PDB. See [RFdiffusion PR#348](https://github.com/RosettaCommons/RFdiffusion/pull/348) for details and the fix. Therefore, `parse_pdbs.py` in `generate_rf_partial_diff_peptide_ensembles/` reorders the chains to place peptide chain as chain 'A'. If using a newer version, this reordering is not necessary.

---

## Data

| File | Description |
|---|---|
| `protein_peptide_dataset.csv` | Full curated set of protein–peptide complexes from PDB |
| `protein-peptide_complexes_candidates.csv` | Complexes passing chain- and complex-level quality filters |
| `protein_peptide_rep_subset.csv` | Representative subset used for ensemble generation |

---

## Citation

If you use this pipeline in your research, please cite the relevant tools:

- **RFdiffusion:** Watson et al., *Nature* (2023). [doi:10.1038/s41586-023-06415-8](https://doi.org/10.1038/s41586-023-06415-8)
- **AlphaFold-Multimer:** Evans et al., *bioRxiv* (2021). [doi:10.1101/2021.10.04.463034](https://doi.org/10.1101/2021.10.04.463034)
- **RCSB PDB:** Berman et al., *Nucleic Acids Research* (2000).

---

## License

This project is for research use. Please check the licenses of RFdiffusion and AlphaFold before deploying in any commercial context.