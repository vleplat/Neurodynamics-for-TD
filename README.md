# Nonnegative Tensor Decomposition Via Collaborative Neurodynamic Optimization
This project is dedicated to the development of Neurodymanics and ODE-based methods for computing Canonical Polyadic Decompositions (CPD) of tensors under constraints. A specific emphasis is put on the nonnegative constraints, but methods can be easily extended to other (simple) constraints such as box, unit ball and unit simplex constraints.

This MATLAB software reproduces the results from the following paper:

```
@misc{SalmanAhmadi-Asl2023,
      title={Nonnegative Tensor Decomposition Via Collaborative Neurodynamic Optimization}, 
      author={Salman Ahmadi-Asl and Valentin Leplat and Anh Huy Phan and Jun Wang and Andrzej Cichocki},
      year={2023},
      eprint={},
      archivePrefix={arXiv},
      primaryClass={cs.LG}
}
```

## Acknowledgements

The baseline algorithms used in the manuscript are courtesy of their respective authors.

## Content
 
 - /libraries : contains helpful libaries.
 
 - /Datasets : contains test data sets.

 - /Utils : contains helpful files, the MatLab implementations of the Algorithms developped in the paper, and MatLab routines to run the demos.
   
 - test files detailed in next Section
   
## Test files
 
 Test files are available. To proceed, open and start one of the following files:
 
- demo_real_data_HSI.m : demo file for standalone HSI tests
- benchmark_HSI.m : run benchmark for Hyperspectral Images data set, see section 7-Example 4 of the paper.
- demo_experiment_1.m : the code associated with Example 1 of Section 7 of the paper.
- demo_experiment_2.m : the code associated with Example 2 of Section 7 of the paper.

