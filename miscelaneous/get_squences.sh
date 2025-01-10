#!/bin/bash

for i in {20..70}
do
    ./org_script.py generate --num-msas 10 --filter-outliers -q "NUM_TAXA = $i AND BRANCH_LENGTH_MEAN >= 0.25 AND BRANCH_LENGTH_MEAN <= 0.75 AND BRANCH_LENGTH_VARIANCE < 0.25"
    ./org_script.py generate --num-msas 10 --filter-outliers -q "NUM_TAXA = $i AND BRANCH_LENGTH_MEAN >= 0.25 AND BRANCH_LENGTH_MEAN <= 0.75 AND BRANCH_LENGTH_VARIANCE < 0.25"
done

# Iterate through the directory and copy the file assembled_sequences.fasta to a new directory with a new name

for dir in out/*
do
    dirname=$(basename "$dir")
    cp "$dir/assembled_sequences.fasta" "sequences/assembled_sequences_${dirname}.fasta"
    cp "$dir/tree_best.newick" "trees/tree_best_${dirname}.newick"
done