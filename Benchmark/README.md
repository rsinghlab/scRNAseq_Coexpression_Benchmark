# Benchmark

This directory contains all the scripts of benchmarks. All the experiment results (i.e., model metrics) will be saved to the ``res`` sub-directory of the project root path.

<!-- You can download experiment results at [here](https://drive.google.com/file/d/1WNwrXj1JnoS4H55R5rpmWboKQyAxvrXK/view?usp=sharing). It should be unzipped and put in the project root directory. -->

## Test Methods on Simulated and Experimental Data

- [1_Benchmark_100hvg.R](./1_Benchmark_100hvg.R): Test methods on NORTA and SERGIO simulations with 100 genes.
- [1_Benchmark_HighDim.R](./1_Benchmark_HighDim.R): Test methods on NORTA (500 genes) and SERGIO (400 genes) simulations.
- [1_Benchmark_OldSim.R](./1_Benchmark_OldSim.R): Test methods on ZI-Gaussian and ZI-Poisson simulations.

## Test Methods on Simulations with Different Settings

- [2_Benchmark_DiffCellNum.R](./2_Benchmark_DiffCellNum.R): Test on simulations with different number of cells.
- [2_Benchmark_DiffGraph.R](./2_Benchmark_DiffGraph.R): Test on simulations with different network topologies.
- [2_Benchmark_DiffSparsity.R](./2_Benchmark_DiffSparsity.R): Test on simulations with different sparsity levels.

## Test Methods on Pre-processed Simulations

- [3_Benchmark_Normalization.R](./3_Benchmark_Normalization.R): Test methods on normalized simulations.
- [3_Benchmark_Pseudobulk.R](./3_Benchmark_Pseudobulk.R): Test methods on pseudo-bulks of simulations.
- [3_Benchmark_Imputation.R](./3_Benchmark_Imputation.R): Test methods on imputed simulations.