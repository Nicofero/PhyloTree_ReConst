import numpy as np
import dimod
from dimod import BinaryQuadraticModel, BINARY
from typing import Optional, Union
from dwave.system import DWaveSampler, EmbeddingComposite
import time
import matplotlib.pyplot as plt
import sys
import re
from Bio import Align
from Bio.Align import substitution_matrices
from Bio import SeqIO
import pandas as pd
from ete3 import Tree

############################################################################################################################

###########################################################################
#                                                                         #
#             ____  _       ____        _   _                             #
#            | __ )(_) ___ |  _ \ _   _| |_| |__   ___  _ __              #
#            |  _ \| |/ _ \| |_) | | | | __| '_ \ / _ \| '_ \             #
#            | |_) | | (_) |  __/| |_| | |_| | | | (_) | | | |            #
#            |____/|_|\___/|_|    \__, |\__|_| |_|\___/|_| |_|            #
#                                 |___/                                   #
#                                                                         #
###########################################################################

# Load sequences from FASTA file
def load_sequences(fasta_file,format='fasta'):
    sequences = []
    labels = []
    for record in SeqIO.parse(fasta_file, format):
        sequences.append(str(record.seq))
        labels.append(str(record.id))
    return sequences,labels

# Compute pairwise alignment scores using BLOSUM62 and normalizing as the paper do
def compute_distance_matrix(sequences):
    aligner = Align.PairwiseAligner()
    blosum62 = substitution_matrices.load("BLOSUM62")
    aligner.mode='global'
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = -10  # Gap opening penalty
    aligner.extend_gap_score = -0.5  # Gap extension penalty
    
    n = len(sequences)
    distance_matrix = np.zeros((n, n))
    bit = [0]*n

    for i in range(n):
        bit[i] = aligner.score(sequences[i], sequences[i])

    for i in range(n):
        for j in range(i):
            sequences[i] = sequences[i].replace('-',"")
            sequences[j] = sequences[j].replace('-',"")
            score = aligner.score(sequences[i], sequences[j])
            # Convert score to a distance (example: max_score - score)
            distance_matrix[i, j] = 100*score/np.mean([bit[i],bit[j]])  
            distance_matrix[j, i] = distance_matrix[i, j] # Symetrical matrix

    return distance_matrix

# Save the distance matrix as a CSV or visualize it
def save_distance_matrix(matrix, labels, output_file):
    df = pd.DataFrame(matrix, index=labels, columns=labels)
    df.to_csv(output_file)
    print(f"Distance matrix saved to {output_file}")
    
def transform_distance(matrix):
    rows = matrix.shape[0]
    dmatrix = []
    names = []
    for i in range(rows):
        names.append(str(i))
        aux = []
        for j in range(i+1):
            aux.append(100 - matrix[i,j])
        dmatrix.append(aux)
    
    return dmatrix,names


############################################################################################################################

##############################################################
#                                                            #
#             ____     __        __                          #
#            |  _ \    \ \      / /_ ___   _____             #
#            | | | |____\ \ /\ / / _` \ \ / / _ \            #
#            | |_| |_____\ V  V / (_| |\ V /  __/            #
#            |____/       \_/\_/ \__,_| \_/ \___|            #
#                                                            #
##############################################################


# Class to create a tree
class TreeNode:
    def __init__(self, value):
        self.value = value
        self.children = []  # A list of child selfs
    
    def print_tree(self, prefix="", is_last=True):
        r"""
        Prints the tree structure with a visually appealing layout.

        Args:
            `node` (TreeNode): The current node to print.
            `prefix` (str): The prefix for the current node (e.g., spaces and vertical bars).
            `is_last` (bool): True if the node is the last child of its parent.
        """
        if self is None:
            return

        # Prettier connectors
        connector = "└── " if is_last else "├── "
        if len(self.value)!=1:
            print(prefix + connector + '|')
        else:
            print(prefix + connector + str(self.value))

        # Update the prefix for child nodes
        if self.children:
            for i, child in enumerate(self.children):
                # Use "│   " for intermediate children, "    " for the last child
                next_prefix = prefix + ("    " if is_last else "│   ")
                child.print_tree(prefix=next_prefix, is_last=(i == len(self.children) - 1))
    
    def _newick_tree(self,exp=[""]):
        r"""
        Prints the tree structure with a visually appealing layout.

        Args:
            `node` (TreeNode): The current node to print.
            `is_last` (bool): True if the node is the last child of its parent.
        """
        if self is None:
            return
        
        # Update the prefix for child nodes
        if self.children:
            exp[0]+='('
            for child in self.children:             
                child._newick_tree(exp)  
            exp[0]+='),'
        else:
            exp[0]+=str(self.value)[1:-1]+','
            
    def to_newick(self):
        exp = [""]
        self._newick_tree(exp)
        exp = re.sub(r",\s*\)", ")",exp[0])
        exp = re.sub(r",$", ";",exp)
        return exp
            
    def create_newick_file(self,file_name:str='tree'):
        r"""
        Creates a newick file from the chosen tree
        
        Args:
            `tree` (TreeNode): The full tree
            `file_name` (str): The name of the file to be written without extension
        """
        exp = [""]
        self._newick_tree(exp)
        exp = re.sub(r",\s*\)", ")",exp[0])
        exp = re.sub(r",$", ";",exp)
        with open(file_name,'w+') as file:
            file.write(exp)
            
    # More visual way
    # Recursive function to calculate positions for nodes
    def _calculate_positions(self, x=0, y=0, dx=1,dy=0.5, positions=None, parent_positions=None, parent=None):
        if positions is None:
            positions = {}
        if parent_positions is None:
            parent_positions = []

        # Assign a position to the current node
        node_label = str(self.value)  # Convert list to string for display
        positions[node_label] = (x, y)
        if parent is not None:
            parent_positions.append((parent, node_label))

        # Layout children
        num_children = len(self.children)
        if num_children > 0:
            child_dx = dx / num_children
            child_x = x - dx / 2 + child_dx / 2
            for child in self.children:
                child._calculate_positions(child_x, y - dy, child_dx, dy, positions, parent_positions, node_label)
                child_x += child_dx

        return positions, parent_positions

    def _display(self,positions,parent_positions,tree_size=1,square=True,names=[]):
        plt.figure(figsize=(10*tree_size, 5*tree_size))
        
        if not square:
            for parent, child in parent_positions:
                x1, y1 = positions[parent]
                x2, y2 = positions[child]
                plt.plot([x1, x2], [y1, y2], "k-")  # Draw edges
        
        if square:
        # Draw edges with straight lines
            for parent, child in parent_positions:
                x1, y1 = positions[parent]
                x2, y2 = positions[child]

                # Horizontal line to connect parent to the child's vertical position
                plt.plot([x1, x2], [y1, y1], "k-")
                # Vertical line to connect child's position to the horizontal line
                plt.plot([x2, x2], [y1, y2], "k-")
        
        positions = {key: value for key, value in positions.items() if len(eval(key)) == 1}
        # Display nodes with single-element values
        for node, (x, y) in positions.items():
            
            if len(eval(node)) == 1:  # Check the list length
                plt.scatter(x, y, s=500, color="lightblue",edgecolors='black')
                if len(names)!= len(positions):
                    if len(node) == 4:
                        plt.text(x, y, node[1:-1], ha="center", weight='bold', va="center", fontsize=10, bbox=dict(facecolor="lightblue", edgecolor="none",boxstyle='circle,pad=0.3'))
                    else:
                        plt.text(x, y, node[1:-1], ha="center", weight='bold', va="center", fontsize=10, bbox=dict(facecolor="lightblue", edgecolor="none",boxstyle='circle,pad=0.5'))
                else:
                    if len(node) == 4:
                        plt.text(x, y, names[int(node[1:-1])], ha="center", weight='bold', va="center", fontsize=10, bbox=dict(facecolor="lightblue", edgecolor="none",boxstyle='circle,pad=0.3'))
                    else:
                        plt.text(x, y, names[int(node[1:-1])], ha="center", weight='bold', va="center", fontsize=10, bbox=dict(facecolor="lightblue", edgecolor="none",boxstyle='circle,pad=0.5'))
                    

        plt.axis("off")
        plt.title("Phylogenetic tree")
        plt.show()  
        
    def display_tree(self,tree_size=1,square=True,names=[]):
        r""""Displays the tree as a graph. Encapsulates 2 functions to easily handle the fancy display.
        
        Args:
            `tree_size` (int): Factor of the renderization size (not aspect ratio)
            `square` (bool): Horizontal and vertical lines (True) or diagonal lines (False)
            `names` (list of str): Names of the nodes       
        """
        pos, par_pos = self._calculate_positions()
        self._display(pos, par_pos,tree_size=tree_size,square=square,names=names)

# Timer class
class Timer:
    def __init__(self,value=0.0) -> None:
        self.value = value
    
    def __str__(self) -> str:
        return f"{self.value}"
    
    def update(self,value):
        self.value += value

# Function to create BinaryQuadraticModel from a numpy matrix
def create_cqm_problem (matrix:np.ndarray,tags=[],c=0,alpha=0)->dimod.ConstrainedQuadraticModel:
    r"""
    Creates a ConstrainedQuadraticModel from a numpy matrix using the Min-cut formulation. Both simmetrical matrices and matrices with 0 above the main diagonal work.
    
    Args:
        `matrix`: Matrix that defines the problem.
        `tags`: Name for the variables. Used in the Philogenetic Tree to mantain tags for each matrix.
        `c`: Number of non zero results wanted.
        `alpha`: Factor to amplify. Higher alpha (>1) achieve the number of `c` of non zero results easily. If 0, no defined `c` is set, thus 11...11 and 00...00 become the best solution.
    
    Returns:
        The ConstrainedQuadraticModel from dimod that defines the problem.
    """
    
    rows = matrix.shape[0]
    var = []
    
    if not tags:
        for i in range(rows):
            var.append(dimod.Binary(i))
    else:
        var = [dimod.Binary(i) for i in tags]
        
    obj = BinaryQuadraticModel({},{},0.0,BINARY)
    
    for i in range(rows):
        for j in range(i):
            obj+=matrix[i,j]*(var[i]-var[j])**2
    
    suma = np.sum(var)      
    # Add restriction term
    obj+=alpha*(suma-c)**2
    
    problem = dimod.ConstrainedQuadraticModel()
    
    problem.set_objective(obj)
    
    return problem

# Creates a BinaryQuadraticModel from a numpy matrix using the Min-cut formulation. Both simmetrical matrices and matrices with 0 above the main diagonal work.
def min_cut_c (matrix:np.ndarray,tags=[],c=0,alpha=0)->BinaryQuadraticModel:
    r"""
    Creates a BinaryQuadraticModel from a numpy matrix using the Min-cut formulation. Both simmetrical matrices and matrices with 0 above the main diagonal work.
    
    Args:
        `matrix`: Matrix that defines the problem.
        `c`: Number of non zero results desired.
        `alpha`: Factor to amplify. Higher alpha (>1) achieve the number of `c` of non zero results easily. If 0, no defined `c` is set, thus 11...11 and 00...00 become the best solution.
    
    Returns:
        The BinaryQuadraticModel from dimod that defines the problem.
    """
    
    rows = matrix.shape[0]
    var = []
    
    if not tags:
        for i in range(rows):
            var.append(dimod.Binary(i))
    else:
        var = [dimod.Binary(i) for i in tags]
    obj = BinaryQuadraticModel({},{},0.0,BINARY)
    
    for i in range(rows):
        for j in range(i):
            obj+=matrix[i,j]*(var[i]-var[j])**2
    
    suma = np.sum(var)      
    # Add restriction term
    obj+=alpha*(suma-c)**2
    
    return obj

# Returns the Ncut from a cut
def n_cut(score,ng0,ng1,og):
    r"""
    Returns the Ncut from a cut.
    
    Args:
        `score`: The energy level obtained from the Mincut method.
        `ng0`: The nodes defining one of the new graphs.
        `ng1`: The nodes defining the other new graph.
        `og`: The original matrix.
    Returns:
        Ncut energy level.
    """
    asoc0 = np.sum(og[ng0, :])
    asoc1 = np.sum(og[ng1, :])   
    
    # print(f'asoc0= {asoc0}, asoc1 = {asoc1}') 
    
    if asoc0 <10e-8 or asoc1 <10e-8:
        ncut = np.inf
    else:
        ncut = (score/asoc0) + (score/asoc1)
    
    return ncut

# Recursive function that uses D-Wave QA to create the Phylogenetic tree using Ncut
def phylo_tree(matrix:np.ndarray,tags=[],**kwargs):
    r"""
    Recursive function that uses D-Wave QA to create the Phylogenetic tree using Ncut
    
    Args:
        `matrix`: The matrix defining the graph.
        `tags`: Tags defining the names of the nodes, used for recursivity. **MUST BE AN INT LIST**
    Returns:
        The `TreeNode` containing the full tree. 
    """
    
    
    ncuts = []
    n_graph_0 = []
    n_graph_1 = []
    sampler = EmbeddingComposite(DWaveSampler())
    
    if not tags:
        sub_mat = matrix
    else:
        sub_mat = matrix[np.ix_(tags, tags)]
        
    rows = sub_mat.shape[0]
    
    alpha = rows*100
    
    var = int(np.floor(rows/2.0))+1

    # Run min_cut for each configuration
    for i in range(1,var):
        # print(f'Corte con {i}')
        if not tags:
            problem = min_cut_c(sub_mat,c=i,alpha=alpha)
        else:
            problem = min_cut_c(sub_mat,tags=tags,c=i,alpha=alpha)
        result = sampler.sample(problem, num_reads=32)
        
        # Time measurement
        if 'timer' in kwargs:
            kwargs['timer'].update(result.info['timing']['qpu_access_time']/1000)
            
        n_graph_0.append([j for j in result.first.sample if result.first.sample[j]==0])
        n_graph_1.append([j for j in result.first.sample if result.first.sample[j]==1])        
        # print(f'\tLa division es: {n_graph_0[i-1]} | {n_graph_1[i-1]}')
        
        if not n_graph_0[i-1] or not n_graph_1[i-1]:
            n_graph_0.pop()
            n_graph_1.pop()
        else:
            ncuts.append(n_cut(result.first.energy,n_graph_0[i-1],n_graph_1[i-1],matrix))
        
    
    # Get the cuts created by the minimum ncut value
    index = np.argmin(ncuts)
    # print(f'Se selecciona la separacion: {n_graph_0[index]} | {n_graph_1[index]}')
    
    node = TreeNode(tags)
    
    # Recursivity in the first graph
    if len(n_graph_0[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(phylo_tree(matrix,n_graph_0[index],timer=kwargs['timer']))
        else:
            node.children.append(phylo_tree(matrix,n_graph_0[index]))
    else:
        leaf = TreeNode(n_graph_0[index])
        if len(n_graph_0[index]) == 2:
            leaf.children.append(TreeNode([n_graph_0[index][0]]))
            leaf.children.append(TreeNode([n_graph_0[index][1]]))
        node.children.append(leaf)
        
    # Recursivity in the first graph
    if len(n_graph_1[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(phylo_tree(matrix,n_graph_1[index],timer=kwargs['timer']))
        else:
            node.children.append(phylo_tree(matrix,n_graph_1[index]))
    else:
        leaf = TreeNode(n_graph_1[index])
        if len(n_graph_1[index]) == 2:
            leaf.children.append(TreeNode([n_graph_1[index][0]]))
            leaf.children.append(TreeNode([n_graph_1[index][1]]))
        node.children.append(leaf)
    
    return node


def bf_tree(matrix,tags=[],**kwargs):
    r"""
    Recursive function that uses D-Wave brute force solver to create the Phylogenetic tree using Ncut.
    
    Args:
        `matrix`: The matrix defining the graph.
        `tags`: Tags defining the names of the nodes, used for recursivity. **MUST BE AN INT LIST**
        
    Returns:
        The `TreeNode` containing the full tree. 
    """
    
    ncuts = []
    n_graph_0 = []
    n_graph_1 = []
    solver = dimod.ExactCQMSolver()
    
    if not tags:
        sub_mat = matrix
    else:
        sub_mat = matrix[np.ix_(tags, tags)]
        
    rows = sub_mat.shape[0]
    
    var = int(np.floor(rows/2.0))+1
    
    alpha = rows*100

    # Run min_cut for each configuration
    for i in range(1,var):
        # print(f'Corte con {i}')
        if 'timer' in kwargs:
            start = time.time_ns()/1000000
            
        if not tags:
            problem = create_cqm_problem(sub_mat,c=i,alpha=alpha)
        else:
            problem = create_cqm_problem(sub_mat,tags=tags,c=i,alpha=alpha)

        sol = solver.sample_cqm(problem)

        # We want the best feasible solution. We can filter by its feasibility and take the first element
        feas_sol = sol.filter(lambda s: s.is_feasible)
        
        if 'timer' in kwargs:
            end = time.time_ns()/1000000
            kwargs['timer'].update(end-start)
            
        n_graph_0.append([j for j in feas_sol.first.sample if feas_sol.first.sample[j]==0])
        n_graph_1.append([j for j in feas_sol.first.sample if feas_sol.first.sample[j]==1])   
        # print(f'\tLa division es: {n_graph_0[i-1]} | {n_graph_1[i-1]}')
        
        if not n_graph_0[i-1] or not n_graph_1[i-1]:
            n_graph_0.pop()
            n_graph_1.pop()
        else:
            ncuts.append(n_cut(feas_sol.first.energy,n_graph_0[i-1],n_graph_1[i-1],matrix))
    
    index = np.argmin(ncuts)
    # print(f'Se selecciona la separacion: {n_graph_0[index]} | {n_graph_1[index]}')
    
    node = TreeNode(tags)
    
    # Recursivity in the first graph
    if len(n_graph_0[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(bf_tree(matrix,n_graph_0[index],timer=kwargs['timer']))
        else:
            node.children.append(bf_tree(matrix,n_graph_0[index]))
    else:
        leaf = TreeNode(n_graph_0[index])
        if len(n_graph_0[index]) == 2:
            leaf.children.append(TreeNode([n_graph_0[index][0]]))
            leaf.children.append(TreeNode([n_graph_0[index][1]]))
        node.children.append(leaf)
        
    # Recursivity in the first graph
    if len(n_graph_1[index]) > 2:
        if 'timer' in kwargs:
            node.children.append(bf_tree(matrix,n_graph_1[index],timer=kwargs['timer']))
        else:
            node.children.append(bf_tree(matrix,n_graph_1[index]))
    else:
        leaf = TreeNode(n_graph_1[index])
        if len(n_graph_1[index]) == 2:
            leaf.children.append(TreeNode([n_graph_1[index][0]]))
            leaf.children.append(TreeNode([n_graph_1[index][1]]))
        node.children.append(leaf)
    
    return node

# Compare trees
def treecmp(tree1:Optional[Union[str , Tree] ],tree2:Optional[Union[str , Tree]]):
    r"""
    Compare two input trees and returns the Robinson-Foulds distance between them
    
    Args:
		`tree1` (Tree or Newick string): First tree.
    	`tree2` (Tree or Newick string): Second tree.
    
    Returns:
		The Robinson-Foulds distance between trees.
  
	Raises:
		Exception.    
    """
    try:
        if type(tree1) == str:
            tree1 = Tree(tree1)
        if type(tree2) == str:
            tree2 = Tree(tree2)
            
        rf, max_parts = tree1.robinson_foulds(tree2, unrooted_trees=True)[:2]
        accuracy = round((1 - (rf / max_parts)) * 100, 2)
        return accuracy
        
    except Exception as e:
        print(str(e), file=sys.stderr)
        
        
####################################################################################
#                                                                                  #
#             __  __ _              _                                              #
#            |  \/  (_)___  ___ ___| | __ _ _ __   ___  ___  _   _ ___             #
#            | |\/| | / __|/ __/ _ \ |/ _` | '_ \ / _ \/ _ \| | | / __|            #
#            | |  | | \__ \ (_|  __/ | (_| | | | |  __/ (_) | |_| \__ \            #
#            |_|  |_|_|___/\___\___|_|\__,_|_| |_|\___|\___/ \__,_|___/            #
#                                                                                  #
####################################################################################

def slice_alignment(file,number,new_name=""):
    r"""Slices alignment of size N to size number, and creates a new phylip file to save it.
    
    Args:
        `file` (str): Name of the file where the alignment is located.
        `number` (int): Number of sequences desired.
        `new_name` (str): Name of the new sliced file. If not given, the same name will be inputed.
    """
        
    with open(file,'r+') as f:
        
        # Get sizes
        fl = f.readline().split()
        size = int(fl[0])
        length = int(fl[1])
        
        
        i = 0
        lines = [""]*size
        for line in f:
            if line == '\n':
                i=0
            else:
                lines[i]+=line.strip()+' '
                i+=1
    
    # New file name
    if new_name != "":
        new_path = "/".join(file.split("/")[:-1] + [new_name])
    else:
        new_path = ".".join(file.split(".")[:-1] + ['_sliced_'+str(number)+'.'+file.split(".")[-1]])
        
    with open(new_path,'w+') as f:
        f.write(' '+str(number)+' '+str(length)+'\n')
        for i in range(number):
            f.write(lines[i]+'\n')