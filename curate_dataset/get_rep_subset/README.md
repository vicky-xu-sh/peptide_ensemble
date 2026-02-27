# Get Representative Subset

This directory contains scripts to generate a representative subset of the protein-peptide dataset by clustering receptor proteins based on sequence similarity.

## Overview

The workflow reduces redundancy in the curated protein-peptide dataset by clustering receptor sequences and selecting representative structures (cluster centroids).

## Workflow

### 1. Extract Receptor Sequences
```bash
./extract_receptor_seq.sh
```
This generates `receptor_sequences.fasta` containing all receptor protein sequences from the dataset.

### 2. Cluster Sequences
```bash
./cluster_receptor_seqs.sh
```
This clusters the receptor sequences using MMseqs2 with a 30% sequence identity threshold. The clustering identifies groups of similar proteins and selects centroid representatives for each cluster.

### 3. Extract Representative PDB IDs
```bash
python extract_representative_pdb_ids.py
```
The MMseqs2 clustering produces a `.tsv` file where the first column contains the centroid PDB IDs (representative structures). This script extracts the unique centroid PDB IDs and saves them to `protein_peptide_rep_subset.csv`.

## Result

The final `protein_peptide_rep_subset.csv` contains a non-redundant subset of the protein-peptide dataset, with one representative structure per cluster at 30% sequence identity.