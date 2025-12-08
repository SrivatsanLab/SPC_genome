#!/bin/bash
#SBATCH -J coalescence_array
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH -t 6:00:00
#SBATCH -p short
#SBATCH --array=0-10000

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <combinations_file> <adata_file> <output_dir> <num_tasks>"
    exit 1
fi

# Assign positional arguments to variables
COMBOS="$1"
ADATAPATH="$2"
OUTPUT_DIR="$3"
NUM_TASKS="$4"

# Run the Python script with the task ID
python scripts/coalescence.py -c $COMBOS -d $ADATAPATH -o $OUTPUT_DIR -t $SLURM_ARRAY_TASK_ID -n $NUM_TASKS

