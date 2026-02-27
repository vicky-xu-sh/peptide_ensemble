#!/bin/bash
#SBATCH --job-name=extract_seq
#SBATCH --output=extract_seq_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=05:00:00

source ~/python_3_11_env/bin/activate

echo "Start extracting receptor sequences..."

python extract_receptor_sequences.py 

echo "Job completed successfully."