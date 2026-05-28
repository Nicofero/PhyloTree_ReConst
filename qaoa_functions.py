"""
    Script Name: qaoa_functions.py
    Author: Nicolás Fernández Otero
    Creation date: 09/12/2024
    Last update: 09/12/2024
    Description: This script contains all the functions for the creation of Phylogenetic trees using a QAOA approach for the optimization problem.
"""

from qiskit.quantum_info import SparsePauliOp,Statevector
from typing import Union
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import real_amplitudes
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from scipy.optimize import minimize
from qa_functions import TreeNode,min_cut_c,n_cut,Timer
from qiskit_algorithms import QAOA
from qiskit_algorithms.optimizers import COBYLA
from qiskit.quantum_info import SparsePauliOp
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit_optimization.translators import from_ising
from qiskit_optimization import QuadraticProgram
from qiskit_aer.primitives import EstimatorV2 as AerEstimator
from qiskit_aer.primitives import SamplerV2 as AerSampler
from qiskit_optimization.algorithms.qrao import (
    QuantumRandomAccessEncoding,
    QuantumRandomAccessOptimizer,
    MagicRounding,
    SemideterministicRounding,
)
from qiskit_optimization.minimum_eigensolvers import VQE
import re
import time

import warnings
import functools

def deprecated(func):
    """Decorator to mark functions as deprecated."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        warnings.warn(
            f"{func.__name__} is deprecated and may be removed in the future.",
            category=DeprecationWarning,
            stacklevel=2
        )
        return func(*args, **kwargs)
    return wrapper

##############################################
#                                            #
#          ___      _    ___    _            #
#         / _ \    / \  / _ \  / \           #
#        | | | |  / _ \| | | |/ _ \          #
#        | |_| | / ___ \ |_| / ___ \         #
#         \__\_\/_/   \_\___/_/   \_\        #
#                                            #
##############################################

#


def interaction_term(qc: QuantumCircuit, phi, control, target)->QuantumCircuit:
    r"""
    Builds the gate from the interaction term in the Ising Hamiltonian
    
    Args:
        `qc`: The quantum circuit to be added
        `phi`: The parameter $\phi$ from the optimization
        `control`: The position of the control qubit
        `target`: The position of the target qubit
    
    Returns:
        The modified QuantumCircuit    
    """
    inverse_target = qc.num_qubits - target - 1
    inverse_control = qc.num_qubits - control - 1
    qc.cx(inverse_control,inverse_target)
    qc.rz(phi*2,inverse_target)
    qc.cx(inverse_control,inverse_target)
    qc.barrier()
    
    return qc



def circuit_size(qc:QuantumCircuit)->int:
    r"""Utility function to get the number of registers of a QuantumCircuit"""
    size = qc.num_qubits
    return size



# Transforms a Pauli expression to a supported quadratic expression
def pauliop_to_exp(op:SparsePauliOp)->str:
    r"""
    Transforms a Pauli expression to a supported quadratic expression
    
    Args:
        `op`: Pauli operation
        
    Returns:
        The str expression    
    """
    
    exp = ''
    for i,pauli in enumerate(op.paulis):
        exp+=str(float(op.coeffs[i].real))
        for i,char in enumerate(pauli.to_label()):
            if char == 'Z':
                exp+='Z'+str(i)
        exp+='+'
    exp = exp[:-1]
    return exp


# Function to automatize layers of an ansatz
def create_ansatz_layer(qc:QuantumCircuit,expression:Union[str,SparsePauliOp],phi=1,beta=1)->QuantumCircuit:
    r"""
    Creates an ansatz layer. However you need to start your circuit with h gates to enter superposition.
    
    Args:
        `qc`: QuantumCircuit to be modified.
        `expression`: Hamiltonian to optimize. Either in str or pauli expression.
        `phi`: Parameter for the Z rotations.
        `beta`: Parameter for the X rotations.
    
    Returns:
        The QuantumCircuit with the added ansatz.
    
    """
    if type(expression) == SparsePauliOp:
        expression = pauliop_to_exp(expression)
    #Read expression
    expression = expression.replace(' ', '')
    
    # Standardize the polynomial string to handle positive terms properly
    # Esto creo que esta mal, habria que comprobarlo
    expression = re.sub(r'(?<=[0-9])(-)(?=[0-9])','+-',expression)
    if expression[0] == '+':
        expression = expression[1:]
    
    # Split the string into terms
    terms = re.split(r'(?<=[0-9])(?:\+)(?=[0-9]|-|Z)',expression)
    
    for term in terms:
        gate = term.count('Z')
        coefs = term.split('Z')
        if coefs[0] == '':
            coefs[0]=1
        if coefs[0] == '-':
            coefs[0]=-1
        if gate == 1:
            qc.rz(2*float(coefs[0])*phi,int(coefs[1]))
            qc.barrier()
        else:
            qc = interaction_term(qc,float(coefs[0])*phi,int(coefs[1]),int(coefs[2]))
            
    size = circuit_size(qc)
    for i in range(size):
        qc.rx(2*beta,i)
    qc.barrier()
    return qc


# Create full ansatz (kind of stupid)
def create_ansatz(expression:str,qubits:int,layers=1,phi=[1],beta=[1]):
    r"""
    Creates an QAOA ansatz.
    
    Args:
        `expression`: Hamiltonian to optimize. Either in str or pauli expression.
        `qubits`: number of qubits of the circuit.
        `layers`:Number of layers to be added.
        `phi`: Parameter for the Z rotations. Needs to be of size == layers.
        `beta`: Parameter for the X rotations. Needs to be of size == layers.
    
    Returns:
        The QuantumCircuit with the added ansatz.
    
    Raises:
        ValueError: In case the size of the parameters is wrong.
    
    """
    qc = QuantumCircuit(qubits)
    for i in range(qubits):
        qc.h(i)
        
    if layers != len(phi) | layers!= len(beta):
        raise ValueError(f'For {layers} layers you need to input size {layers} parameters')
        
    for i in range(layers):
        qc = create_ansatz_layer(qc,expression,phi[i],beta[i])
    return qc

# Eval energy from the expression and value
@deprecated
def eval_energy(expression:Union[str,SparsePauliOp],factor:str):
    r"""
    Evaluates the energy of the factor given. This is the energy expectation from the state.
    
    Args:
        `expression`: Hamiltonian to optimize. Either in str or pauli expression.
        `factor`: Binary string to calculate the expectation value for the Hamiltonian.
    
    Returns:
    The energy value
    
    """
    energy = 0
    
    if type(expression) == SparsePauliOp:
        expression = pauliop_to_exp(expression)
        
    expression = expression.replace(' ', '')
    
    # Standardize the polynomial string to handle positive terms properly
    expression = re.sub(r'(?<=[0-9])(-)(?=[0-9])','+-',expression)
    if expression[0] == '+':
        expression = expression[1:]
    
    # Split the string into terms
    terms = re.split(r'(?<=[0-9])(?:\+)(?=[0-9]|-|Z)',expression)
    
    
    for term in terms:
        gate = term.count('Z')
        coefs = term.split('Z')
        if coefs[0] == '':
            coefs[0]=1
        if coefs[0] == '-':
            coefs[0]=-1
        if gate == 1:
            energy += (-2*int(factor[int(coefs[1])])+1)*float(coefs[0])
        else:
            energy += (2*np.abs(int(factor[int(coefs[1])])+int(factor[int(coefs[2])])-1)-1)*float(coefs[0])
    
    return energy

# Get energy full
@deprecated
def get_energy(qc:QuantumCircuit,expression,shots=1024):
    r"""
    Returns the energy from an execution of a QuantumCircuit
    
    Args:
        `qc`: QuantumCircuit to run.
        `shots`: Number of shots for the circuit.
        
    Returns:
        The energy value.
    """
    
    # Define the simulator, in a future version, this would be a parameter
    sim = AerSimulator()
    
    # Transpile the circuit for the simulator or real QPU
    qc.measure_all()
    qc = transpile(qc,sim)

    # Run the circuit and collect results
    sampler = Sampler(mode = AerSimulator())
    job = sampler.run([qc],shots=shots)
    job_result = job.result()
    counts=job_result[0].data.meas.get_counts()

    # Using the formula from [1]
    energy = 0
    for key in counts:
        energy+= (counts[key]/shots)*eval_energy(expression,key)
    
    return energy

# Theoretical perfection
@deprecated
def get_energy_statevector(qc:QuantumCircuit,expression):
    r"""
    Returns the energy from an execution of a QuantumCircuit
    
    Args:
        `qc`: QuantumCircuit to run.
        `expression`: Expression of the Hamiltonian.
        `size`: Size of the circuit.
        
    Returns:
        The energy value.
    """
    
    # Define the simulator, in a future version, this would be a parameter
    st = Statevector(qc)
    size = circuit_size(qc)
    # Using the formula from [1]
    energy = 0
    for i in range(2**size):
        energy+= (st[i]**2).real*eval_energy(expression,f"{i:0{size}b}")
    
    return energy


# Now the important part, the minimization. We will use a class

class MyQAOA:    
    
    def __init__(self,exp:Union[str,SparsePauliOp],size:int,layers=1,method='COBYLA',shots=1024,x0=None,exact=False,backend=None):
        if type(exp) == str:
            self.exp = exp
        else:
            self.exp = pauliop_to_exp(exp)
        self.size = size
        self.layers = layers
        self.method = method
        self.param = [0.]*2*layers
        self.min = np.inf
        self.shots = shots
        self.counts = None
        self.exact = exact
        self.qc = None
        if backend:
            self.backend = backend
        else:
            self.backend = AerSimulator()
            
        if x0 is None:
            self.x0 = np.random.random(layers*2) 
        elif layers == len(x0)/2:
            self.x0 = x0
        else:
            self.x0 = np.random.random(layers*2)
        self.eval = self.eval_energy_all()
        self.pauli = to_pauli_string(self.exp,self.size)
    
    def eval_energy_all(self):
        r"""
        Evaluates all the energies. This is the energy expectation from the state.
        
        Returns:
        A dictionary containing all the energy expectations
        
        """
        
        factors = [format(i, f'0{self.size}b') for i in range(2**self.size)]
        energy_exp = {}
        for factor in factors:
            energy = 0
                            
            self.exp = self.exp.replace(' ', '')
            
            # Standardize the polynomial string to handle positive terms properly
            self.exp = re.sub(r'(?<=[0-9])(-)(?=[0-9])','+-',self.exp)
            if self.exp[0] == '+':
                self.exp = self.exp[1:]
            
            # Split the string into terms
            terms = re.split(r'(?<=[0-9])(?:\+)(?=[0-9]|-|Z)',self.exp)
            
            
            for term in terms:
                gate = term.count('Z')
                coefs = term.split('Z')
                if coefs[0] == '':
                    coefs[0]=1
                if coefs[0] == '-':
                    coefs[0]=-1
                if gate == 1:
                    energy += (-2*int(factor[int(coefs[1])])+1)*float(coefs[0])
                else:
                    energy += (2*np.abs(int(factor[int(coefs[1])])+int(factor[int(coefs[2])])-1)-1)*float(coefs[0])
                    
            energy_exp[factor] = energy
        return energy_exp
    
    def get_energy(self,qc:QuantumCircuit):
        r"""
        Returns the energy from an execution of a QuantumCircuit
        
        Args:
            `qc`: QuantumCircuit to run.
            
        Returns:
            The energy value.
        """
               
        # Transpile the circuit for the simulator or real QPU
        qc.measure_all()
        qc = transpile(qc,self.backend)

        # Run the circuit and collect results
        sampler = Sampler(self.backend)
        job = sampler.run([qc],shots=self.shots)
        job_result = job.result()
        counts=job_result[0].data.meas.get_counts()

        # Using the formula from [1]
        # energy = 0
        # for key in counts:
        #     energy+= (counts[key]/self.shots)*self.eval[key]
        state = np.array([0]*(2**self.size),dtype=float)
        for key in counts:
            state[int(key,2)] = np.sqrt(counts[key]/self.shots)
        energy = state @ self.pauli.to_matrix().real @ state

        # Save circuit and counts
        if energy < self.min:
            self.qc = qc
            self.counts = counts
            self.min = energy
        
        return energy
    
    def get_energy_statevector(self,qc:QuantumCircuit):
        r"""
        Returns the energy from an execution of a QuantumCircuit
        
        Args:
            `qc`: QuantumCircuit to run.
            `expression`: Expression of the Hamiltonian.
            `size`: Size of the circuit.
            
        Returns:
            The energy value.
        """
        
        # Define the simulator, in a future version, this would be a parameter
        st = Statevector(qc)
        size = circuit_size(qc)
        # Using the formula from [1]
        energy = 0
        for i in range(2**size):
            energy+= (np.abs(st[i])**2)*self.eval[f"{i:0{size}b}"]
            
        if energy < self.min:
            self.qc = qc
            self.counts = {f"{i:0{size}b}":np.abs(st[i])**2 for i in range(2**size)}
            self.min = energy
        return energy
    
    def objective (self,x):
    
        # Setting up gamma and beta
        gamma = x[:self.layers]
        beta = x[self.layers:]
        
        qc = create_ansatz(self.exp,self.size,self.layers,gamma,beta)
        
        if not self.exact:
            energy = self.get_energy(qc)
        else:
            energy = self.get_energy_statevector(qc)
        return energy
    
    # @dask.delayed
    def get_min(self):
        # Run minimization 
        bds = tuple((0,2*np.pi) for _ in range(self.layers)) + tuple((0,np.pi) for _ in range(self.layers))
        result = minimize(self.objective,self.x0,method=self.method,bounds=bds)
        self.param = result.x
        self.min = result.fun
        return self.min,self.param
        
    def get_opt_circ(self):

        self.get_min()
        self.qc = create_ansatz(self.exp,self.size,self.layers,self.param[:self.layers],self.param[self.layers:])
        
        return self.qc
    
# Prepare expression
def prepare_exp(matrix:np.ndarray,c=0):
    r"""
    Prepares the Hamiltonian needed for the QAOA circuit
    
    Args:
        `matrix` (np.ndarray): Numpy array containing the symmetric (or lower triangular matrix) that defines the problem.
        `c` (int): Number of non-zero variables.
        `tags` (list): Defines the name of the variables. **MUST BE AN INT LIST**
        
    Returns:
        The expression as a string.        
    """
    
    
    alpha = matrix.shape[0]*100
    h, J, _ = min_cut_c(matrix,c=c,alpha=alpha).to_ising()
    # print(f'Linar coeffs: {h}')
    # print(f'Quadratic terms: {J}')
    exp = ""

    for term in J:
        exp+="+"+str(J[term])+"Z"+str(term[0])+"Z"+str(term[1])

    for term in h:
        if h[term]>0:
            exp+="+"+str(h[term])+"Z"+str(term)
        elif h[term]<0:
            exp+=str(h[term])+"Z"+str(term)
            
    return exp 

# Combination of solutions
def combine_inverse_keys(original_dict):
    """
    Combines the reciprocal solutions as only one solution
    """
    new_dict = {}
    
    for key, value in original_dict.items():
        # Compute the inverse key
        inverse_key = ''.join('1' if bit == '0' else '0' for bit in key)
        
        # Sort the pair to ensure consistent representation
        combined_key = min(key, inverse_key)
        
        # Add value to the combined key
        if combined_key in new_dict:
            new_dict[combined_key] += value
        else:
            new_dict[combined_key] = value
    
    return new_dict

# Best solution
def get_bes_sol(qc:QuantumCircuit)->dict:
    """   Runs the quantic circuit simulation and returns the value with more counts   """
    sim = AerSimulator()

    # Transpile the circuit for the simulator or real QPU
    qc.measure_all()
    qc = transpile(qc,sim)

    # Run the circuit and collect results
    sampler = Sampler(mode = AerSimulator())
    job = sampler.run([qc],shots=2048)
    job_result = job.result()
    counts=job_result[0].data.meas.get_counts()
    n_counts = combine_inverse_keys(counts)
    best = max(n_counts,key=n_counts.get)
    return best

def qaoa_phylo_tree(matrix:np.ndarray,tags=[],client=None,backend=AerSimulator(),**kwargs):
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
            if client:  
                problem = prepare_exp(sub_mat,c=i)
                num_problems = len(client.scheduler_info()['workers'])
                instances = [MyQAOA(problem,rows,layers=3) for _ in range(num_problems)]
                futures = [client.submit(instance.get_min,pure=False) for instance in instances]
                # Gather results as they complete
                results = client.gather(futures)
                
                best_tot = min(results)
                best = best_tot[1]
                minim = best_tot[0]
                qc = create_ansatz(problem,rows,phi=[best[0]],beta=[best[1]])
                result = get_bes_sol(qc)
            else:
                problem = prepare_exp(sub_mat,c=i)
                qaoa = MyQAOA(problem,rows,layers=1,backend=backend)
                qaoa.get_min()
                minim = qaoa.min
                counts = combine_inverse_keys(qaoa.counts)
                result = max(counts,key=counts.get)

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
    # print(f'Se selecciona la separacion: {n_graph_0[index]} | {n_graph_1[index]}')
    
    node = TreeNode(tags)
    
    # Recursivity in the first graph
    if len(n_graph_0[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(qaoa_phylo_tree(matrix,tags=n_graph_0[index],client=client,backend=backend,timer=kwargs['timer']))
        else:
            node.children.append(qaoa_phylo_tree(matrix,tags=n_graph_0[index],client=client,backend=backend))
    else:
        leaf = TreeNode(n_graph_0[index])
        if len(n_graph_0[index]) == 2:
            leaf.children.append(TreeNode([n_graph_0[index][0]]))
            leaf.children.append(TreeNode([n_graph_0[index][1]]))
        node.children.append(leaf)
        
    # Recursivity in the first graph
    if len(n_graph_1[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(qaoa_phylo_tree(matrix,tags=n_graph_1[index],client=client,backend=backend,timer=kwargs['timer']))
        else:
            node.children.append(qaoa_phylo_tree(matrix,tags=n_graph_1[index],client=client,backend=backend))
    else:
        leaf = TreeNode(n_graph_1[index])
        if len(n_graph_1[index]) == 2:
            leaf.children.append(TreeNode([n_graph_1[index][0]]))
            leaf.children.append(TreeNode([n_graph_1[index][1]]))
        node.children.append(leaf)
    
    return node


def bqm_to_quadratic_program(bqm) -> QuadraticProgram:
    # Ensure BINARY vartype (not SPIN)
    bqm_binary = bqm.change_vartype("BINARY", inplace=False)

    qp = QuadraticProgram()

    # Add binary variables
    for var in bqm_binary.variables:
        qp.binary_var(str(var))

    # Build linear and quadratic dicts
    linear    = {str(v): bias for v, bias in bqm_binary.linear.items()}
    quadratic = {(str(u), str(v)): bias for (u, v), bias in bqm_binary.quadratic.items()}

    # Set objective (BQMs are minimization problems)
    qp.minimize(
        linear=linear,
        quadratic=quadratic,
        constant=bqm_binary.offset
    )

    return qp

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

def qaoa_phylo_tree_qiskit(matrix:np.ndarray,tags=[],backend=AerSimulator(),layers = 1,**kwargs):
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
    
    alpha = rows *300
    
    var = int(np.floor(rows/2.0))+1
    
    sampler = Sampler(mode=backend)    
    
    pm = generate_preset_pass_manager(optimization_level=3, backend=backend)
    
    qaoa = QAOA(sampler=sampler, optimizer=COBYLA(maxiter=300), reps=layers, transpiler=pm)
                
    eigen_optimizer = MinimumEigenOptimizer(min_eigen_solver = qaoa)
    
    while not ncuts:
        
        n_graph_0 = []
        n_graph_1 = []
        # Run min_cut for each configuration
        for i in range(1,var):
            # print(f'Corte con {i}')
            if 'timer' in kwargs:
                start = time.time_ns()/1000000
            # Prepare the expression and run the QAOA    
            problem = min_cut_c(sub_mat,c=i,alpha=alpha)            
            
            qp = bqm_to_quadratic_program(problem)
                
            # Compute the minimum eigenvalue (i.e., approximate ground state)
            rest = eigen_optimizer.solve(qp)

            # Extract the measurement distribution
            result = [str(int(x)) for x in rest.x]
            minim = rest.fval
                    
            # Time measurement
            if 'timer' in kwargs:
                end = time.time_ns()/1000000
                kwargs['timer'].update(end-start)
                
            n_graph_0.append([tags[j] for j in range(len(result)) if result[j]=='0'])
            n_graph_1.append([tags[j] for j in range(len(result)) if result[j]=='1'])        
            # print(f'\tLa division es: {n_graph_0[i-1]} | {n_graph_1[i-1]}')
            
            # print(n_cut(minim,n_graph_0[i-1],n_graph_1[i-1],matrix))
            
            if n_graph_0[i-1] and n_graph_1[i-1]:
                ncuts.append(n_cut(minim,n_graph_0[i-1],n_graph_1[i-1],matrix))
                
    
    # Get the cuts created by the minimum ncut value
    index = np.argmin(ncuts)
    # print(f'Se selecciona la separacion: {n_graph_0[index]} | {n_graph_1[index]}')
    
    node = TreeNode(tags)
    
    # Recursivity in the first graph
    if len(n_graph_0[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(qaoa_phylo_tree_qiskit(matrix,tags=n_graph_0[index],backend=backend,layers=layers,timer=kwargs['timer']))
        else:
            node.children.append(qaoa_phylo_tree_qiskit(matrix,tags=n_graph_0[index],backend=backend,layers=layers))
    else:
        leaf = TreeNode(n_graph_0[index])
        if len(n_graph_0[index]) == 2:
            leaf.children.append(TreeNode([n_graph_0[index][0]]))
            leaf.children.append(TreeNode([n_graph_0[index][1]]))
        node.children.append(leaf)
        
    # Recursivity in the first graph
    if len(n_graph_1[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(qaoa_phylo_tree_qiskit(matrix,tags=n_graph_1[index],backend=backend,layers=layers,timer=kwargs['timer']))
        else:
            node.children.append(qaoa_phylo_tree_qiskit(matrix,tags=n_graph_1[index],backend=backend,layers=layers))
    else:
        leaf = TreeNode(n_graph_1[index])
        if len(n_graph_1[index]) == 2:
            leaf.children.append(TreeNode([n_graph_1[index][0]]))
            leaf.children.append(TreeNode([n_graph_1[index][1]]))
        node.children.append(leaf)
    
    return node

##################################################      
#         ____  _____            ____            #
#        / __ \|  __ \     /\   / __ \           #
#       | |  | | |__) |   /  \ | |  | |          #
#       | |  | |  _  /   / /\ \| |  | |          #
#       | |__| | | \ \  / ____ \ |__| |          #
#        \___\_\_|  \_\/_/    \_\____/           #
#                                                #
##################################################

def min_cut_qp(matrix: np.ndarray, tags=[]) -> QuadraticProgram:
    r"""
    Creates a QuadraticProgram from a numpy matrix using the Min-cut formulation.
    Both symmetrical matrices and matrices with 0 above the main diagonal work.

    Objective (minimise):
        sum_{i>j} w_ij * (x_i - x_j)^2
        = sum_{i>j} w_ij * (x_i + x_j - 2*x_i*x_j)    [since x_i^2 = x_i for binary]

    Args:
        matrix: NxN weight matrix defining the problem.
        tags:   Optional variable names. Defaults to "0", "1", ..., "N-1".

    Returns:
        QuadraticProgram encoding the min-cut objective.
    """
    rows = matrix.shape[0]
    qp = QuadraticProgram()

    var_names = [str(t) for t in tags] if tags else [str(i) for i in range(rows)]

    for name in var_names:
        qp.binary_var(name)

    linear = {}
    quadratic = {}

    for i in range(rows):
        for j in range(i):          # lower triangle only — handles both symmetric
            w = matrix[i, j]        # and upper-triangle-zero matrices correctly
            if w == 0:
                continue

            # Linear part: w*(x_i + x_j)
            linear[var_names[i]] = linear.get(var_names[i], 0) + w
            linear[var_names[j]] = linear.get(var_names[j], 0) + w

            # Quadratic part: -2w * x_i * x_j
            key = (var_names[i], var_names[j])
            quadratic[key] = quadratic.get(key, 0) - 2 * w

    qp.minimize(linear=linear, quadratic=quadratic)
    return qp

import numpy as np
from scipy import stats

def threshold_similarity_matrix(W: np.ndarray, method: str = "otsu") -> tuple[np.ndarray, float]:
    """
    Sparsify a similarity matrix by zeroing out edges below a data-driven threshold.
    
    Args:
        W:      Symmetric similarity matrix (NxN), diagonal ignored.
        method: One of "otsu", "percentile", "zscore", "knee"
    
    Returns:
        W_sparse: Thresholded matrix
        threshold: The threshold value used
    """
    # Extract upper-triangle values (avoid diagonal and duplicate edges)
    idx = np.triu_indices_from(W, k=1)
    weights = W[idx]

    if method == "otsu":
        threshold = _otsu_threshold(weights)

    elif method == "percentile":
        # Keep the top 50% of edges by default — adjust as needed
        threshold = np.percentile(weights, 50)

    elif method == "zscore":
        # Zero out anything more than 1 std below the mean
        threshold = weights.mean() - weights.std()

    elif method == "knee":
        threshold = _knee_threshold(weights)

    else:
        raise ValueError(f"Unknown method: {method}")

    # Apply threshold symmetrically
    W_sparse = W.copy()
    W_sparse[W_sparse < threshold] = 0
    np.fill_diagonal(W_sparse, 0)  # Nodes don't connect to themselves in QRAO
    return W_sparse, threshold


def _otsu_threshold(values: np.ndarray) -> float:
    """
    Otsu's method: finds the threshold that minimises intra-class variance.
    Works well when similarity values are bimodal (a cluster of weak edges
    and a cluster of strong edges).
    """
    bins = 256
    counts, bin_edges = np.histogram(values, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    total = counts.sum()

    best_thresh, best_var = 0, -1
    w0_sum = 0
    mu0_sum = 0

    for i in range(1, bins):
        w0 = counts[:i].sum() / total
        w1 = 1 - w0
        if w0 == 0 or w1 == 0:
            continue
        mu0 = (counts[:i] * bin_centers[:i]).sum() / (counts[:i].sum() + 1e-12)
        mu1 = (counts[i:] * bin_centers[i:]).sum() / (counts[i:].sum() + 1e-12)
        between_var = w0 * w1 * (mu0 - mu1) ** 2
        if between_var > best_var:
            best_var = between_var
            best_thresh = bin_centers[i]

    return best_thresh


def _knee_threshold(values: np.ndarray) -> float:
    """
    Knee/elbow detection on the sorted weight curve.
    Finds where the curve bends — a natural 'drop-off' point.
    """
    sorted_vals = np.sort(values)[::-1]  # Descending
    n = len(sorted_vals)
    x = np.linspace(0, 1, n)

    # Normalise
    y = (sorted_vals - sorted_vals.min()) / (sorted_vals.max() - sorted_vals.min() + 1e-12)

    # Distance from each point to the line connecting first and last
    line_vec = np.array([x[-1] - x[0], y[-1] - y[0]])
    line_vec /= np.linalg.norm(line_vec)
    point_vecs = np.stack([x - x[0], y - y[0]], axis=1)
    distances = np.abs(np.cross(line_vec, point_vecs))
    knee_idx = np.argmax(distances)

    return sorted_vals[knee_idx]

def _eval(qp: QuadraticProgram, x: np.ndarray) -> float:
    """Evaluate the QP objective at a binary vector."""
    return qp.objective.evaluate(x)


def _project_to_cardinality(
    x: np.ndarray, c: int, qp: QuadraticProgram
) -> tuple[np.ndarray, float]:
    """
    Greedily flip bits one at a time until exactly c are 1,
    choosing each flip to keep the objective as low as possible.
    """
    x = x.copy().astype(float)

    while int(x.sum()) != c:
        excess = int(x.sum()) - c          # positive → too many 1s
        candidates = np.where(x == (1 if excess > 0 else 0))[0]
        target_value = 0.0 if excess > 0 else 1.0

        best_idx = min(
            candidates,
            key=lambda idx: _eval(qp, np.where(
                np.arange(len(x)) == idx, target_value, x
            )),
        )
        x[best_idx] = target_value

    return x, _eval(qp, x)


def _best_with_cardinality(
    samples, c: int, qp: QuadraticProgram
) -> tuple[np.ndarray, float]:
    """
    Return the best solution with exactly c ones.

    Strategy:
      1. Filter samples that already satisfy the constraint → pick best fval.
      2. If none exist, project every sample to cardinality c and pick best.
    """
    feasible = [s for s in samples if int(s.x.sum()) == c]

    if feasible:
        best = min(feasible, key=lambda s: s.fval)   # minimisation → lowest fval
        return best.x.copy().astype(float), float(best.fval)

    # No sample satisfies the constraint: project all and keep the best
    best_x, best_fval = None, float("inf")
    for sample in samples:
        x_proj, fval = _project_to_cardinality(sample.x, c, qp)
        if fval < best_fval:
            best_fval = fval
            best_x = x_proj

    return best_x, best_fval


# ── main function ──────────────────────────────────────────────────────────────

def run_qrao_min_cut(
    matrix: np.ndarray,
    c: int,
    sparsity_method = "zscore", # "zscore" | "otsu" | "knee" | "percentile"
    tags: list = [],
    max_vars_per_qubit: int = 3,
    rounding: str = "magic",   # "magic" | "semideterministic"
    shots: int = 10_000,
    seed: int = 42,
) -> dict:
    """
    Solve a min-cut problem with QRAO, post-processing the result
    to enforce exactly c variables equal to 1.

    Args:
        matrix:             NxN weight matrix (symmetric or lower-triangle).
        c:                  Required number of selected assets (ones in solution).
        tags:               Optional variable names; defaults to "0".."N-1".
        max_vars_per_qubit: QRAC compression level (1, 2, or 3).
        rounding:           "magic" (multiple samples) or "semideterministic".
        shots:              Number of magic rounding shots.
        seed:               RNG seed for reproducibility.

    Returns:
        dict with keys:
            x                — binary solution array with exactly c ones
            fval             — objective value at x
            raw_result       — full QuantumRandomAccessOptimizationResult
            num_qubits       — qubits used after compression
            compression_ratio
            feasible_found   — whether any sample satisfied the constraint directly
    """
    # ── 1. Sparsify the matrix and build the QUBO ──────────────────────────────────────────────────────
    W_sparse,thresh = threshold_similarity_matrix(matrix, method=sparsity_method)
    
    qp = min_cut_qp(W_sparse,tags)

    # ── 2. Encode ──────────────────────────────────────────────────────────────
    encoding = QuantumRandomAccessEncoding(max_vars_per_qubit=max_vars_per_qubit)
    encoding.encode(qp)
    # print(
    #     f"Encoding: {encoding.num_vars} variables → {encoding.num_qubits} qubits "
    #     f"(compression {encoding.compression_ratio:.2f}x)"
    # )

    # ── 3. Solver (VQE) ────────────────────────────────────────────────────────
    pm = generate_preset_pass_manager(optimization_level=3, backend=AerSimulator())
    estimator = AerEstimator(
                        options={
                            "backend_options": {
                                "method": "statevector",   # exact simulation — same as StatevectorEstimator
                                "device": "CPU",           # swap to "GPU" if you have a CUDA-capable card
                                "max_parallel_threads": 0, # 0 = use all available cores
                                "max_parallel_experiments": 0,
                            }
                        }
)
    ansatz = real_amplitudes(encoding.num_qubits)
    vqe = VQE(ansatz=ansatz, optimizer=COBYLA(maxiter=300), estimator=estimator, pass_manager= pm)

    # ── 4. Rounding scheme ─────────────────────────────────────────────────────
    if rounding == "magic":
        sampler = AerSampler(
                    options={
                        "backend_options": {
                            "method": "statevector",
                            "max_parallel_threads": 0, # 0 = use all available cores
                            "max_parallel_experiments": 0,
                        }
                    }
                )
        rounding_scheme = MagicRounding(sampler=sampler, pass_manager=pm)
    elif rounding == "semideterministic":
        rounding_scheme = SemideterministicRounding()
    else:
        raise ValueError(f"Unknown rounding '{rounding}'. Use 'magic' or 'semideterministic'.")

    # ── 5. Solve ───────────────────────────────────────────────────────────────
    qrao = QuantumRandomAccessOptimizer(
        min_eigen_solver=vqe,
        rounding_scheme=rounding_scheme,
    )
    raw_result = qrao.solve(qp)

    # ── 6. Post-process: enforce exactly c ones ────────────────────────────────
    feasible_found = any(int(s.x.sum()) == c for s in raw_result.samples)
    x, fval = _best_with_cardinality(raw_result.samples, c, qp)

    if not feasible_found:
        print(
            f"No sample had exactly {c} ones — "
            f"projected best sample to cardinality {c}."
        )

    return {
        "x":                  x,
        "fval":               fval,
        "raw_result":         raw_result,
        "num_qubits":         encoding.num_qubits,
        "compression_ratio":  encoding.compression_ratio,
        "feasible_found":     feasible_found,
    }



def qrao_phylo_tree_qiskit(matrix:np.ndarray,tags=[],backend=AerSimulator(),**kwargs):
    r"""
    Recursive function that uses QRAO to create the Phylogenetic tree using Ncut
    
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
    
    while not ncuts:
        
        n_graph_0 = []
        n_graph_1 = []
        # Run min_cut for each configuration
        for i in range(1,var):
            # print(f'Corte con {i}')
            if 'timer' in kwargs:
                start = time.time_ns()/1000000
            # Prepare the expression and run the QRAO    
            res = run_qrao_min_cut(sub_mat,c = i)
            
            result = [str(int(x)) for x in res['x']]
            minim = res['fval']
                    
            # Time measurement
            if 'timer' in kwargs:
                end = time.time_ns()/1000000
                kwargs['timer'].update(end-start)
                
            n_graph_0.append([tags[j] for j in range(len(result)) if result[j]=='0'])
            n_graph_1.append([tags[j] for j in range(len(result)) if result[j]=='1'])        
            # print(f'\tLa division es: {n_graph_0[i-1]} | {n_graph_1[i-1]}')
            
            # print(n_cut(minim,n_graph_0[i-1],n_graph_1[i-1],matrix))
            
            if n_graph_0[i-1] and n_graph_1[i-1]:
                ncuts.append(n_cut(minim,n_graph_0[i-1],n_graph_1[i-1],matrix))
                
    
    # Get the cuts created by the minimum ncut value
    index = np.argmin(ncuts)
    # print(f'Se selecciona la separacion: {n_graph_0[index]} | {n_graph_1[index]}')
    
    node = TreeNode(tags)
    
    # Recursivity in the first graph
    if len(n_graph_0[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(qrao_phylo_tree_qiskit(matrix,tags=n_graph_0[index],backend=backend,timer=kwargs['timer']))
        else:
            node.children.append(qrao_phylo_tree_qiskit(matrix,tags=n_graph_0[index],backend=backend,))
    else:
        leaf = TreeNode(n_graph_0[index])
        if len(n_graph_0[index]) == 2:
            leaf.children.append(TreeNode([n_graph_0[index][0]]))
            leaf.children.append(TreeNode([n_graph_0[index][1]]))
        node.children.append(leaf)
        
    # Recursivity in the first graph
    if len(n_graph_1[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(qrao_phylo_tree_qiskit(matrix,tags=n_graph_1[index],backend=backend,timer=kwargs['timer']))
        else:
            node.children.append(qrao_phylo_tree_qiskit(matrix,tags=n_graph_1[index],backend=backend,))
    else:
        leaf = TreeNode(n_graph_1[index])
        if len(n_graph_1[index]) == 2:
            leaf.children.append(TreeNode([n_graph_1[index][0]]))
            leaf.children.append(TreeNode([n_graph_1[index][1]]))
        node.children.append(leaf)
    
    return node