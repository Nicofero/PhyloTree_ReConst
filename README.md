# **Reconstruction of Phylogenetic Trees using Quantum Computing**

This repository contains the complete codebase developed for my Thesis.  
The goal of the thesis is to improve **phylogenetic tree reconstruction** by leveraging **quantum computing techniques**, specifically:

- **Quantum Annealing (QA)** using D-Wave systems  
- **Quantum Approximate Optimization Algorithm (QAOA)**

Both approaches are applied to the reconstruction of phylogenetic trees from biological sequence alignments and are benchmarked using multiple quantitative metrics.

---

## **Repository Structure Overview**

The project is organized into Jupyter notebooks, Python modules, and structured directories to facilitate experimentation, benchmarking, and reproducibility.

---

## **Installation and Usage**

### 1. (Optional but Recommended) Create a Virtual Environment

```bash
python -m venv phylo_tree
```

#### Activate the environment

**Windows**
```bash
phylo_tree\Scripts\activate
```

**Linux / macOS**
```bash
source phylo_tree/bin/activate
```

---

### 2. Install Dependencies

```bash
pip install -r requirements.txt
```

---

### 3. Run the Notebooks

All experiments, reconstructions, and analyses are primarily conducted through the provided Jupyter notebooks.

---

## **File Descriptions**

### Notebooks

- **qa_phylo_tree.ipynb**  
  Constructs phylogenetic trees using a **D-Wave Quantum Annealer**.

- **qaoa_phylo_tree.ipynb**  
  Constructs phylogenetic trees using the **Quantum Approximate Optimization Algorithm (QAOA)**.

- **metrics_study.ipynb**  
  Exploratory analysis aimed at identifying, extracting, and evaluating metrics from the generated trees.

---

### Python Modules

- **qa_functions.py**  
  Core utility functions for Quantum Annealing–based phylogenetic reconstruction.

- **qaoa_functions.py**  
  Core utility functions for QAOA-based phylogenetic reconstruction.

- **qa_big_trees.py**  
  Script to execute the **NMcutQA algorithm** on new or large-scale alignments.

- **qa_hybrid_trees.py**  
  Variant of the above script implementing a **hybrid quantum–classical approach**.

- **requirements.txt**  
  List of Python dependencies required to reproduce the experimental environment.

---

## **Directory Descriptions**

- **/alignments**  
  Contains aligned protein or nucleotide sequences.  
  - Files are ordered by number of taxa.  
  - Supported formats: `phylip`, `fasta`.

- **/matrices**  
  Stores similarity matrices derived from alignments to avoid recomputation.  
  The `distance_matrices` subdirectory contains distance matrices obtained by inverting similarity values.

- **/trees**  
  Output phylogenetic trees in **Newick format**, ordered by number of taxa.

- **/metrics**  
  Stores extracted metrics such as execution times, distances, and reconstruction quality indicators.

- **/benchmarking_results**  
  Notebooks and results related to benchmarking Quantum Annealing and hybrid approaches.

- **/miscellaneous**  
  Experimental files, prototypes, and tests not currently integrated into the main workflow.

---

## **Reproducibility Notes**

- Some experiments require access to **D-Wave quantum hardware or hybrid solvers**.
- Results may vary depending on solver configuration, backend availability, and random seeds.

---

## **Main References**

1. Combarro, E. F., & Gonzalez-Castillo, S. (2023). *A Practical Guide to Quantum Machine Learning and Quantum Optimisation*. Packt Publishing.

2. Onodera, W., Hara, N., Aoki, S., Asahi, T., & Sawamura, N. (2022). *Phylogenetic tree reconstruction via graph cut presented using a quantum-inspired computer*. Molecular Phylogenetics and Evolution, 178, 107636. https://doi.org/10.1016/j.ympev.2022.107636

3. Farhi, E., Goldstone, J., & Gutmann, S. (2014). *A Quantum Approximate Optimization Algorithm*. arXiv:1411.4028

---

## **License**

This project is intended for academic and research purposes.  
Please consult the accompanying license file for usage and redistribution terms.
