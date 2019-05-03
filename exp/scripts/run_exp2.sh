#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=8000M
#SBATCH --account=rrg-chauvec
#SBATCH --output /home/chauvec/wg-anoph/COUNTING_TREES/DLT-Histories/DLTcount/exp/traces/exp_2.trace
#SBATCH --error  /home/chauvec/wg-anoph/COUNTING_TREES/DLT-Histories/DLTcount/exp/traces/exp_2.error

cd /home/chauvec/wg-anoph/COUNTING_TREES/DLT-Histories/DLTcount/exp/scripts
python exp2_count_histories_unranked_trees.py
