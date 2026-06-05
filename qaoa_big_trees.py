from qaoa_functions import *
import re
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-f',"--folder", type=str, help="Folder name containing the matrices")
parser.add_argument('-m',"--method", type=str, choices=['qaoa','qrao','pce','sparse_pce', 'qrao_esu2'], help="Method for the experiment")
# parser.add_argument('--ansatz', type=str, choices=['real_amplitudes', 'efficientSU2'], default='real_amplitudes', help="Ansatz for the QRAO method")
args = parser.parse_args()

folder = args.folder

if args.method == 'qaoa':
    method_name = "nmcutqaoa"
    timer_name = "timer_qaoa"
    
elif args.method == 'qrao':    
    method_name = "nmcutqrao"
    timer_name = "timer_qrao"

elif args.method == 'sparse_pce':
    method_name = "nmcutspce"
    timer_name = "timer_spce"

elif args.method == 'qrao_esu2':
    method_name = "nmcutqrao_esu2"
    timer_name = "timer_qrao_esu2"
else:
    method_name = "nmcutpce"
    timer_name = "timer_pce"

with open(f'./benchmarking_results/{folder}/{timer_name}_{folder}.csv','w') as fp:
    fp.write(f'size,name,method,time\n')

files = os.listdir(f'./benchmarking matrices/{folder}')
for file in files:
    if not file.endswith("_matrix.txt"):
        continue  # skip non-matrix files
    file_ext = file[:-len("_matrix.txt")]
    
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
    
    # For QAOA
    if args.method in ['qaoa','qrao','qrao_esu2']:
        if distance_matrix.shape[0] <26:
            
            print(f'PROCESSING {file} of size {distance_matrix.shape[0]} for {folder}...\n')
        else:
            print(f'SKIPPING {file} because size {distance_matrix.shape[0]} is too big for {args.method}...\n')
            continue
    
    # Define parallel AerSimulator
    backend = AerSimulator(
        method="statevector",
        max_parallel_threads=64,        # use all CPU cores
        max_parallel_experiments=0,    # parallelize circuits
        runtime_parameter_bind_enable=True  
    )
    
    estimator = AerEstimator(
        options={
                                "backend_options": {
                                    "method": "statevector",   # exact simulation — same as StatevectorEstimator
                                    "device": "CPU",           # swap to "GPU" if you have a CUDA-capable card
                                    "max_parallel_threads": 0, # 0 = use all available cores
                                    "max_parallel_experiments":0
                                }
                            }
    )
        
    # Begin NMcutQAOA
    timer = Timer(0.0)
    try:
        # tree_qa = qaoa_phylo_tree_qiskit(distance_matrix,timer=timer,layers=3,backend=backend)   # --For QAOA testing--
        if args.method == 'qaoa':
            tree_qa = qaoa_phylo_tree_qiskit(distance_matrix, timer=timer, layers = 3, backend=backend)

        elif args.method == 'qrao':
            tree_qa = qrao_phylo_tree_qiskit(distance_matrix, timer=timer)

        elif args.method == 'qrao_esu2':
                tree_qa = qrao_phylo_tree_qiskit(distance_matrix, timer=timer, ansatz="efficientSU2")

        elif args.method == 'pce':
            tree_qa = pce_phylo_tree_qiskit(distance_matrix, timer=timer, estimator= estimator)
        else:
            tree_qa = sparse_pce_phylo_tree_qiskit(distance_matrix, timer=timer, estimator= estimator)

        with open(f'./benchmarking_results/{folder}/{timer_name}_{folder}.csv','a') as fp:
            # fp.write(f'{distance_matrix.shape[0]},{file},NMcutQAOA,{timer.value}\n')
            fp.write(f'{distance_matrix.shape[0]},{file},{method_name},{timer.value}\n')

        with open(f'./benchmarking matrices/{folder}/{file_ext}_species.txt','r') as f:
            species = [line.strip() for line in f.readlines()]

        new_file = f'{file_ext}_{method_name}_tree.txt'
        tree_qa.create_newick_file(f'./benchmarking_results/{folder}/{new_file}',labels=species)
        print(f'Finished {file} in {timer.value}ms\n---------------------------------------------------------\n')
        del distance_matrix
        del tree_qa
        del species
    except Exception as e:
        print(f'Error processing {file}: {e}')
        print(f'Moving to next file...\n---------------------------------------------------------\n')
        del distance_matrix

    
    