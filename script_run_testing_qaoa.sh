#!/bin/bash
#
#    Modify this to your needs
#
#SBATCH -c 64
#SBATCH -t 24:0:0
#SBATCH --mem-per-cpu=3G

module load cesga/2020 miniconda3/22.11.1-1

conda activate phylo-tree

python qaoa_trees.py -n 64