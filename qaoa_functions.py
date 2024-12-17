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
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import SamplerV2
from scipy.optimize import minimize
from qa_functions import TreeNode,min_cut_c,n_cut,Timer
import dask

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
    
    qc.cx(control,target)
    qc.rz(phi*2,target)
    qc.cx(control,target)
    qc.barrier()
    
    return qc



def circuit_size(qc:QuantumCircuit)->int:
    r"""Utility function to get the number of registers of a QuantumCircuit"""
    size = 0
    for reg in qc.qregs:
        size+=reg.size
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
    expression = expression.replace('-', '+-')
    if expression[0] == '+':
        expression = expression[1:]
    
    # Split the string into terms
    terms = expression.split('+')
    
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
    Creates an ansatz. However you need to start your circuit with h gates to enter superposition.
    
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
    expression = expression.replace('-', '+-')
    if expression[0] == '+':
        expression = expression[1:]
    
    # Split the string into terms
    terms = expression.split('+')
    
    
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
    sampler = SamplerV2()
    job = sampler.run([qc],shots=shots)
    job_result = job.result()
    counts=job_result[0].data.meas.get_counts()

    # Using the formula from [1]
    energy = 0
    for key in counts:
        energy+= (counts[key]/shots)*eval_energy(expression,key)
    
    return energy

# Theoretical perfection
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

class QAOA:    
    
    def __init__(self,exp:Union[str,SparsePauliOp],size:int,layers=1,method='COBYLA',shots=1024,x0=[0.,0.]):
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
        if layers == len(x0)/2:
            self.x0 = x0
        else:
            self.x0 = np.random.random(layers*2)
        
        
    def objective (self,x):
    
        # Setting up gamma and beta
        gamma = x[:self.layers]
        beta = x[self.layers:]
        
        qc = create_ansatz(self.exp,self.size,self.layers,gamma,beta)
        
        energy = get_energy(qc,self.exp,shots=self.shots)
        
        return energy
    
    @dask.delayed
    def get_min(self):
        # Run minimization 
        result = minimize(self.objective,self.x0,method=self.method)
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
        `matrix` (np.ndarray): Numpy array containing the simmetric (or lower triangular matrix) that defines the problem.
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
    sampler = SamplerV2()
    job = sampler.run([qc],shots=2048)
    job_result = job.result()
    counts=job_result[0].data.meas.get_counts()
    n_counts = combine_inverse_keys(counts)
    best = max(n_counts,key=n_counts.get)
    return best