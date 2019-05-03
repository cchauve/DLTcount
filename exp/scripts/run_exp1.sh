#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=8000M
#SBATCH --account=rrg-chauvec
#SBATCH --output /home/chauvec/wg-anoph/COUNTING_TREES/DLT-Histories/DLTcount/exp/traces/exp_1.trace
#SBATCH --error  /home/chauvec/wg-anoph/COUNTING_TREES/DLT-Histories/DLTcount/exp/traces/exp_1.error

cd /home/chauvec/wg-anoph/COUNTING_TREES/DLT-Histories/DLTcount/exp/scripts
python exp1_generate_trees.py
