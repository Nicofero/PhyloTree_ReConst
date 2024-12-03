# Reconstruction of Phylogenetic Trees
 This repository contains all the code from my Thesis. The thesis tries to improve the reconstruction of Phylogenetic Trees by using Quantum Computing. The reconstruction of the trees will be made using both Quantum Annealers and QAOA.

## How to use
- Create a virtual or conda environment (Optional but highly recomended)
  ```
  python -m venv phylo_tree
  ```
  - Run the enviroment
    ```
    phylo_tree\Scripts\activate
    ```

- Install the requirements:
  ```
  pip install -r requirements.txt
  ```

- Use the notebooks

## Description of the files

- **qaoa_test.ipynb:** This notebook contains the different tests for creating a gate level QAOA algorithm.
- **qa_pyhlo_tree.ipynb:** This notebook cointains the creation of phylogenetic trees using D-Wave Quantum Annealer.  
- **qa_functions.py:** Python library that contains the functions for the QA.
- **test.ipynb:** Random tests.

### Description of the directories

- **/alignments:** Contains the aligned sequences of proteins or nucleotides. The files are in `phylip` or `fasta` format.
- **/trees:** Contains the trees created by the algorithm in newick format.

## Main references

- [1] Combarro, E. F., & Gonzalez-Castillo, S. (2023). A practical guide to quantum machine learning and quantum optimisation: Hands-On Approach to Modern Quantum Algorithms. Packt Publishing.
- [2] Onodera, W., Hara, N., Aoki, S., Asahi, T., & Sawamura, N. (2022). Phylogenetic tree reconstruction via graph cut presented using a quantum-inspired computer. Molecular Phylogenetics and Evolution, 178, 107636. https://doi.org/10.1016/j.ympev.2022.107636
