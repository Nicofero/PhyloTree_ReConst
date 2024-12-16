import numpy as np
import dimod
from dimod import BinaryQuadraticModel, BINARY
from typing import Optional
import time
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random
import re
from qaoa_functions import *
from qa_functions import TreeNode,min_cut_c,n_cut,Timer
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
from qiskit.quantum_info import Pauli,SparsePauliOp,Statevector
from typing import Union
from colorama import Fore
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import SamplerV2
from qiskit.visualization import plot_histogram
from scipy.optimize import minimize
import dask
from distributed import Client
import dask.array as da
import argparse
from multiprocessing import freeze_support

def main():
    # Initialize parser
    parser = argparse.ArgumentParser(description="Input the size of the cluster.")

    # Add arguments
    parser.add_argument("-n","--number", type=int, choices=range(1,65), required=True, help="Number of cores of the cluster (An integer from 1 to 64)")
    args = parser.parse_args()

    num_workers = args.number

    print(f'Selected {num_workers} workers')

    with Client(n_workers=num_workers) as client:


        matrix = np.array( [[ 0, 92, 73, 78, 92],
                            [92,  0, 21, 49, 34],
                            [73, 21,  0, 35, 63],
                            [78, 49, 35,  0, 29],
                            [92, 34, 63, 29,  0]])

        problem = prepare_exp(matrix,c=1)
        size = max([int(x) for x in re.findall(r'Z([0-9]+)',problem)])+1

        result = [None for i in range(num_workers)]
        instances = [QAOA(problem,size) for _ in range(num_workers)]
        start = time.time()
        for i in range(num_workers):
            result[i]=instances[i].get_min()            
        cond_list = dask.compute(result)
        cond_list = cond_list[0]
        end = time.time()
        
    best = min(cond_list)
    print('Time used:',end-start)
    print('Solution:',best)
    best = best[1]
    qc = create_ansatz(problem,size,1,[best[0]], [best[1]])
    sim = AerSimulator()

    # Transpile the circuit for the simulator or real QPU
    qc.measure_all()
    qc = transpile(qc,sim)

    # Run the circuit and collect results
    sampler = SamplerV2()
    job = sampler.run([qc],shots=2048)
    job_result = job.result()
    counts=job_result[0].data.meas.get_counts()
    n_counts = combine_inverse_keys(counts)
    plot_histogram(n_counts)
    plt.savefig('test.png', bbox_inches='tight')
    
if __name__ == '__main__':
    freeze_support()
    main()