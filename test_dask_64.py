import dask.config
import dask.config
from qa_functions import *
from qaoa_functions import *
import time
import re
from qiskit.visualization import plot_histogram
import multiprocessing
import dask
from distributed import Client
import argparse
from multiprocessing import freeze_support

def main():
    # Initialize parser
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
    print(client.dashboard_link)
    time.sleep(10)
    # print(dask.config.get('scheduler'))    
    timer = Timer(0.0)
    matrix = np.array( [[ 0, 92, 73, 78, 92],
                        [92,  0, 21, 49, 34],
                        [73, 21,  0, 35, 63],
                        [78, 49, 35,  0, 29],
                        [92, 34, 63, 29,  0]])
    tree = qaoa_phylo_tree(matrix,client=client,timer=timer)
    print(f'El tiempo utilizado es: {timer}')
    tree.display_tree()
    client.close()
    
if __name__ == '__main__':
    # freeze_support()
    main()