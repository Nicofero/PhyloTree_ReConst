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