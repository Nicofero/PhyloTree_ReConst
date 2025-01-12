from qa_functions import *
from qaoa_functions import *
import time
import re
import dask
from distributed import Client
import os
import argparse
import multiprocessing

def main():
    parser = argparse.ArgumentParser(description="Input the size of the cluster.")
    multiprocessing.set_start_method('spawn', force=True)
    # Add arguments
    parser.add_argument("-n","--number", type=int, choices=range(1,65), required=True, help="Number of cores of the cluster (An integer from 1 to 64)")
    args = parser.parse_args()

    num_workers = args.number

    print(f'Selected {num_workers} workers')
    dask.config.set(scheduler='threads')
    client = Client(n_workers=num_workers)
    # print(client.scheduler_info())
    # print(client.dashboard_link)
    time.sleep(10)

    file_names = os.listdir('matrices')

    for folder in file_names:
        files = os.listdir(f'matrices/{folder}')
        for file in files:
            # Load the sequences
            distance_matrix = np.load(f'./matrices/{folder}/{file}')

            print(f'Generating tree from: {file}')
            # Create the tree
            timer = Timer(0.0)
            tree_qaoa = qaoa_phylo_tree(distance_matrix,client=client,timer=timer)

            with open('timer.csv','a') as fp:
                fp.write(f'{folder},{file},{timer.value}\n')

            file_ext = re.search(r'_([0-9]+)\.',file).group(1)
            new_file = f'qaoa_tree_{file_ext}.newick'
            tree_qaoa.create_newick_file(f'./trees/{folder}/{new_file}')
            print(f'Finished {file} in {timer.value}ms')

            del tree_qaoa
            del distance_matrix


if __name__ == '__main__':
    main()