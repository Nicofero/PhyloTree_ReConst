from qa_functions import *
import re
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-f',"--folder", type=str, help="Folder name containing the matrices")
args = parser.parse_args()

folder = args.folder
completed = []

with open(f'./benchmarking_results/{folder}/timer_hybrid_{folder}.csv','r') as fp:
    for i,line in enumerate(fp):
        if i == 0:
            continue
        aux = line.split(',')[1].strip()[:-len("_matrix.txt")]
        completed.append(aux)

files = os.listdir(f'./benchmarking matrices/{folder}')

with open('token.txt','r') as f:
    TOKEN = f.read()
sampler = LeapHybridBQMSampler(token = TOKEN, cache = False)

for i,file in enumerate(files):
    # if i == 0:
    #     continue  # skip first file if needed
    if not file.endswith("_matrix.txt"):
        continue  # skip non-matrix files
    file_ext = file[:-len("_matrix.txt")]
    
    if file_ext in completed:
        print(f'Skipping {file} as it is already completed.\n---------------------------------------------------------\n')
        continue  # skip already completed files
    
    # if folder == "Bos_taurus" and int(re.search(r'\d{4}', file_ext)[0]) < 6123:
    #     continue  # skip files with size less than 6123 for Bos_taurus
    # elif folder == "Penicillium_digitatum" and int(re.search(r'\d{4}', file_ext)[0]) < 6365:
    #     continue  # skip files with size less than 6365 for Penicillium_digitatum
    # elif folder == "Schizosaccharomyces_pombe" and int(re.search(r'\d{4}', file_ext)[0]) < 4988:
    #     continue  # skip files with size less than 4988 for Schizosaccharomyces_pombe
    # elif folder == "Puccinia_graminis" and int(re.search(r'\d{4}', file_ext)[0]) < 7112:
    #     continue  # skip files with size less than 7112 for Puccinia_graminis
    
    # Load distance matrix
    distance_matrix = np.loadtxt(f'./benchmarking matrices/{folder}/{file}')
    print(f'Processing {file} of size {distance_matrix.shape[0]} for {folder}...\n')
        
    # Begin NMcutQA
    timer = Timer(0.0)
    try:
        tree_qa = hybrid_phylo_tree(distance_matrix,sampler = sampler, timer=timer,name=folder)   # --- For advantage 1 ---

        with open(f'./benchmarking_results/{folder}/timer_hybrid_{folder}.csv','a') as fp:
            fp.write(f'{distance_matrix.shape[0]},{file},NMcutHy,{timer.value}\n')

        with open(f'./benchmarking matrices/{folder}/{file_ext }_species.txt','r') as f:
            species = [line.strip() for line in f.readlines()]

        new_file = f'{file_ext}_nmcuthy_tree.txt'
        tree_qa.create_newick_file(f'./benchmarking_results/{folder}/{new_file}',labels=species)
        print(f'Finished {file} in {timer.value}ms\n---------------------------------------------------------\n')
        del distance_matrix
        del tree_qa
        del species
    except Exception as e:
        print(f'Error processing {file}: {e}')
        print(f'Moving to next file...\n---------------------------------------------------------\n')
        del distance_matrix
