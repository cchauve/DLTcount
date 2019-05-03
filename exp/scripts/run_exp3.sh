#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=8000M
#SBATCH --account=rrg-chauvec
#SBATCH --output /home/chauvec/wg-anoph/COUNTING_TREES/DLT-Histories/DLTcount/exp/traces/exp_3.trace
#SBATCH --error  /home/chauvec/wg-anoph/COUNTING_TREES/DLT-Histories/DLTcount/exp/traces/exp_3.error

cd /home/chauvec/wg-anoph/COUNTING_TREES/DLT-Histories/DLTcount/exp/scripts
python exp3_count_histories_ranked_trees.py

