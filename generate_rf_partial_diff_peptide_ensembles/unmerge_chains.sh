#!/bin/bash

#SBATCH --job-name=pdb_unmerge           
#SBATCH --output=unmerge_%j.log          
#SBATCH --nodes=1                       
#SBATCH --cpus-per-task=16               
#SBATCH --mem=64G                        
#SBATCH --time=04:00:00                 


source ~/python_3_env/bin/activate  

ORIG_DIR="/home/vickyxu/scratch/peptide_ensemble/subset_cleaned_pdb_files"
GEN_DIR="/home/vickyxu/scratch/peptide_ensemble/subset_partial_diff_results"
OUT_DIR="/home/vickyxu/scratch/peptide_ensemble/subset_peptide_ensembles"

python /home/vickyxu/scratch/peptide_ensemble/generate_rf_partial_diff_peptide_ensembles/unmerge_chains.py \
    --orig_pdb_dir "$ORIG_DIR" \
    --generated_pdbs_dir "$GEN_DIR" \
    --output_dir "$OUT_DIR" \
    --num_workers $SLURM_CPUS_PER_TASK

echo "Job completed at $(date)"