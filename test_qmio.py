from qiskit import transpile, QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import QAOA
from qmiotools.integrations.qiskitqmio import QmioBackend,FakeQmio
from qiskit.primitives import BackendSampler
from qiskit_algorithms.optimizers import COBYLA
from qa_functions import *
from qaoa_functions import *
import re
import os
import logging

def to_pauli_string(expr: str,size: int) -> SparsePauliOp:
    """Convert a string expression to a SparsePauliOp.
    
    Args:
        expr (str): The expression in the form of '+coeff*Z0Z1+coeff*Z1Z2+...'.
        size (int): The number of qubits in the system.
    
    Returns:
        SparsePauliOp: The corresponding SparsePauliOp representation.
    """
    # Handle the negative coefficients and ensure the expression is well-formed
    expr = re.sub(r'(?<=[0-9])(-)(?=[0-9])','+-',expr)
    if expr[0] == '+':
        expr = expr[1:]
    # Handle exponential like E+10 and E-10
    expr = re.sub(r'([+-]?\d*\.?\d+)([Ee][+-]?\d+)', r'\1\2', expr)
    expr = expr.replace(' ', '')  # Remove any spaces for easier parsing
    terms = re.split(r'(?<=[0-9])(?:\+)(?=[0-9]|-|Z)', expr)
    pauli_terms = []
    for term in terms:
        coeffs = term.split('Z')
        coeffs[0] = float(coeffs[0]) if coeffs[0] else 1.0
        pauli_string = 'I' * size  # Assuming size qubits
        if coeffs[1:]:
            for qubit in coeffs[1:]:
                index = int(qubit)
                pauli_string = pauli_string[:index] + 'Z' + pauli_string[index + 1:]
        pauli_terms.append((pauli_string, coeffs[0]))
    return SparsePauliOp.from_list(pauli_terms)


def qaoa_phylo_tree_qmio(matrix:np.ndarray,tags=[],**kwargs):
    r"""
    Recursive function that uses QAOA to create the Phylogenetic tree using Ncut
    
    Args:
        `matrix`: The matrix defining the graph.
        `tags`: Tags defining the names of the nodes, used for recursivity. **MUST BE AN INT LIST**
    Returns:
        The `TreeNode` containing the full tree. 
    """
    ncuts = []
    
    if not tags:
        sub_mat = matrix
        tags = list(range(matrix.shape[0]))
    else:
        sub_mat = matrix[np.ix_(tags, tags)]
        
    rows = sub_mat.shape[0]
    
    var = int(np.floor(rows/2.0))+1
    
    sampler = BackendSampler(backend=FakeQmio())
    result = '0'*rows
    # Repeat until a cut is found
    while not ncuts:
        
        n_graph_0 = []
        n_graph_1 = []
        # Run min_cut for each configuration
        for i in range(1,var):
            # print(f'Corte con {i}')
            if 'timer' in kwargs:
                start = time.time_ns()/1000000
                # Prepare the expression and run the QAOA    
            problem = prepare_exp(sub_mat,c=i)
            pstring = to_pauli_string(problem,rows)
            while result == '0'*rows or result == '1'*rows:
                qaoa = QAOA(sampler=sampler, optimizer=COBYLA(maxiter=20), reps=1)

                # Compute the minimum eigenvalue (i.e., approximate ground state)
                rest = qaoa.compute_minimum_eigenvalue(pstring)

                # Extract the measurement distribution
                result = rest.best_measurement['bitstring']
                minim = rest.eigenvalue.real
                    
            # Time measurement
            if 'timer' in kwargs:
                end = time.time_ns()/1000000
                kwargs['timer'].update(end-start)
                
            n_graph_0.append([tags[j] for j in range(len(result)) if result[j]=='0'])
            n_graph_1.append([tags[j] for j in range(len(result)) if result[j]=='1'])        
            # print(f'\tLa division es: {n_graph_0[i-1]} | {n_graph_1[i-1]}')
            
            if n_graph_0[i-1] and n_graph_1[i-1]:
                ncuts.append(n_cut(minim,n_graph_0[i-1],n_graph_1[i-1],matrix))
                
        
    
    # Get the cuts created by the minimum ncut value
    index = np.argmin(ncuts)
    print(f'Se selecciona la separacion: {n_graph_0[index]} | {n_graph_1[index]}')
    
    node = TreeNode(tags)
    
    # Recursivity in the first graph
    if len(n_graph_0[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(qaoa_phylo_tree_qmio(matrix,tags=n_graph_0[index],timer=kwargs['timer']))
        else:
            node.children.append(qaoa_phylo_tree_qmio(matrix,tags=n_graph_0[index]))
    else:
        leaf = TreeNode(n_graph_0[index])
        if len(n_graph_0[index]) == 2:
            leaf.children.append(TreeNode([n_graph_0[index][0]]))
            leaf.children.append(TreeNode([n_graph_0[index][1]]))
        node.children.append(leaf)
        
    # Recursivity in the first graph
    if len(n_graph_1[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(qaoa_phylo_tree_qmio(matrix,tags=n_graph_1[index],timer=kwargs['timer']))
        else:
            node.children.append(qaoa_phylo_tree_qmio(matrix,tags=n_graph_1[index]))
    else:
        leaf = TreeNode(n_graph_1[index])
        if len(n_graph_1[index]) == 2:
            leaf.children.append(TreeNode([n_graph_1[index][0]]))
            leaf.children.append(TreeNode([n_graph_1[index][1]]))
        node.children.append(leaf)
    
    return node

# expr = "+2454.0Z0Z1+2489.5Z2Z1+2463.5Z2Z0+2475.5Z3Z1+2461.0Z3Z0+2482.5Z3Z2+2483.0Z4Z1+2454.0Z4Z0+2468.5Z4Z2+2485.5Z4Z3+2500.0Z1+2500.0Z0+2500.0Z2+2500.0Z3+2500.0Z4"
# hamiltonian = to_pauli_string(expr, 5)

# # Create the QAOA instance
# sampler = BackendSampler(backend=QmioBackend())
# qaoa = QAOA(sampler=sampler, optimizer=COBYLA(maxiter=30), reps=1)

# # Compute the minimum eigenvalue (i.e., approximate ground state)
# result = qaoa.compute_minimum_eigenvalue(hamiltonian)

# # Extract the measurement distribution
# best_result = result.best_measurement['bitstring']
# print("Best result:", best_result)

# # Convert bitstring to array of bits
# # bit_array = [int(b) for b in most_likely_bitstring]

# # Print the result
# print("Estimated ground state energy:", result.eigenvalue.real)


# Test with size 8 matrix
file_names = os.listdir('matrices')

# Remove the folder test
file_names = file_names[:-1]

for folder in file_names:
    if int(folder) == 9:
        files = os.listdir(f'matrices/{folder}')
        for file in files:            
            
            # Load the sequences and create simmilarity matrix
            distance_matrix = np.load(f'matrices/{folder}/{file}')
            print(f'matrices/{folder}/{file}')
            # Begin NMcutQA
            timer = Timer(0.0)
            tree_qa = qaoa_phylo_tree_qmio(distance_matrix,timer=timer)

            file_ext = re.search(r'_([0-9]+)\.',file).group(1)  
            new_file = f'qaoa_qmio_tree_{file_ext}.newick'
            tree_qa.create_newick_file(f'./trees/{folder}/{new_file}')
            print(f'Finished {file} in {timer.value}ms')

            # Erase memory
            del distance_matrix
            del tree_qa
