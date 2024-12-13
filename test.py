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
import multiprocessing

# Standalone function for calling get_min
def call_get_min(instance):
    return instance.get_min()

# External parallelization function
def parallel_minimization(instances:QAOA):
    with multiprocessing.Pool(processes=len(instances)) as pool:
        # Use pool to call get_min on each instance
        results = pool.map(call_get_min, instances)
    return results

if __name__ == '__main__':
    
    
    matrix = np.array( [[ 0, 92, 73, 78, 92],
                    [92,  0, 21, 49, 34],
                    [73, 21,  0, 35, 63],
                    [78, 49, 35,  0, 29],
                    [92, 34, 63, 29,  0]])
    problem = prepare_exp(matrix,c=1)
    size = max([int(x) for x in re.findall(r'Z([0-9]+)',problem)])+1
    
    poolsize = multiprocessing.cpu_count()
    instances = [QAOA(problem,size) for _ in range(10)]
    start = time.time()
    try:
        results = parallel_minimization(instances)
    except Exception as e:
        print(e)
    end = time.time()
    print('Results:',results)
    print('time:',end-start)

    # act_min = np.inf
    # start = time.time()
    # for i in range(10):
    #     initial_guess = np.random.random(2)
    #     qaoa = QAOA(problem,size,x0=initial_guess)
    #     qaoa.get_opt_circ()
    #     if qaoa.min< act_min:
    #         act_min = qaoa.min
    #         arg = qaoa.param
    #         j=i
    #         print(f'El nuevo minimo es {act_min} con {arg}')
    # end = time.time()
    # print(f'El minimo es: {act_min}, y se da con los parametros: {arg}')
    # print('time:',end-start)

    best = min(results)[1]
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
    plt.savefig('test.pdf', bbox_inches='tight')