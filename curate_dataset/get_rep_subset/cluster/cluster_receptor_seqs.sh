#!/bin/bash
#SBATCH --job-name=cluster_test
#SBATCH --output=cluster_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=05:00:00

module load mmseqs2

# --- DEFINE FILENAMES ---
INPUT_FASTA="receptor_sequences.fasta"
OUTPUT="cluster_receptor_seqs"
TMP_FILES="tmp"

echo "Starting Clustering..."
mmseqs easy-cluster $INPUT_FASTA $OUTPUT $TMP_FILES \
    --threads $SLURM_CPUS_PER_TASK \
    --min-seq-id 0.3

echo "Job completed successfully."


