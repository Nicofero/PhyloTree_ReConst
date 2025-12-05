from qa_functions import *
import re
import os

file_names = os.listdir('./alignments')

# Remove the folder test
file_names.remove('test')

for folder in file_names:
    if int(folder) > 99 and int(folder) < 130:
        files = os.listdir(f'./alignments/{folder}')
        for file in files:
            file_ext = re.search(r"_([0-9]+)\.",file).group(1)
            try:
                
                matrix_file = f'./matrices/{folder}/matrix_{file_ext}.npy'
                distance_matrix = np.load(matrix_file)
                
            except:
                print(f'alignments/{folder}/{file}')
                # Load the sequences and create similarity matrix
                sequences,labels = load_sequences(f'./alignments/{folder}/{file}',format='fasta')
                distance_matrix = compute_distance_matrix(sequences)

                # Save the distance matrix
                matrix_file = new_file = f'matrix_{file_ext}'
                np.save(f'./matrices/{folder}/{matrix_file}',distance_matrix)
                print(f'Alignment file generated, generating the phylogenetic tree...')
                # Erase memory
                del sequences
                del labels
                
            # Begin NMcutQA
            timer = Timer(0.0)
            tree_qa = phylo_tree(distance_matrix,timer=timer)

            with open('./metrics/timer_advantage2.csv','a') as fp:
                fp.write(f'{folder},{file},{timer.value}\n')

            new_file = f'qa2_tree_{file_ext}.newick'
            tree_qa.create_newick_file(f'./trees/{folder}/{new_file}')
            print(f'Finished {file} in {timer.value}ms\n---------------------------------------------------------\n')

            
            del distance_matrix
            del tree_qa