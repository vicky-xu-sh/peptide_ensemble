#!/bin/bash
#SBATCH --job-name=rf_partial_diff
#SBATCH --output=partial_diff_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1   
#SBATCH --time=01:30:00

OUTPUT_PREFIX=
INPUT_PDB=
NUM_DESIGNS=
PD_STEPS=
CONTIG_STR=

# Load Apptainer
module load apptainer

# Run the container
apptainer run --nv -C -W $SLURM_TMPDIR -B $HOME/scratch \
    $HOME/scratch/RFdiffusion/rfd.sif \
    inference.output_prefix=$OUTPUT_PREFIX \
    inference.model_directory_path=/app/RFdiffusion/models \
    inference.schedule_directory_path=schedules \
    inference.input_pdb=$INPUT_PDB \
    inference.num_designs=$NUM_DESIGNS \
    contigmap.contigs="$CONTIG_STR" \
    diffuser.partial_T=$PD_STEPS 