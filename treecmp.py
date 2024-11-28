import sys
from typing import Optional
try:
	from ete3 import Tree
except ImportError:
	print("Module ete3 required. Please install with pip install ete3", file=sys.stderr)
	exit(-1)


def treecmp(tree1:Optional[str | Tree ],tree2:Optional[str | Tree ]):
    """
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
