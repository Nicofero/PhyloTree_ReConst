#!/bin/bash
#
#    Modify this to your needs
#
#SBATCH -c 64
#SBATCH -t 6:0:0
#SBATCH --mem-per-cpu=1G

module load cesga/2020 miniconda3/22.11.1-1

conda activate phylo-tree

echo "Programa para los arboles pequeños con 64"

python qaoa_trees.py -n 64