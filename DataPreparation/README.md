# Data Preparation

This directory contains all the scripts used for preparing data in our study. All the data will be saved to the ``data``directory.

There is no need to re-run all the pre-processing scripts, you can download all the data from [here](https://drive.google.com/file/d/1Tmw9mKN20jTcsVcDVxIgbM3xV_7ozDGw/view?usp=sharing). It should be unzipped and put in the project root directory.


## Experimental Data Preparation

We use eight publicly available scRNA-seq datasets in this study. 
We extract data of four different sequencing protocols from [mouse cortex](https://singlecell.broadinstitute.org/single_cell/study/SCP425/single-cell-comparison-cortex-data\#study-summary) dataset and [human peripheral blood mononuclear cells (PBMC)](https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data\#study-summary) dataset downloaded from [Single Cell Protal](https://singlecell.broadinstitute.org/single_cell). 

- [1_ExperimentalPreprocessing.py](./1_ExperimentalPreprocessing.py): Pre-process the experimental datasets.
- [1_ExtractRefNet.py](./1_ExtractRefNet.py): Obtain known gene expression networks.


## Data Simulation

- [2_Old_Simulation.R](./2_Old_Simulation.R): ZI-Gaussian and ZI-Poisson simulations.
- [2_NORTA_Simulation.R](./2_NORTA_Simulation.R): NORTA simulation.
- [2_NORTA_Simulation_DiffSetting.R](./2_NORTA_Simulation_DiffSetting.R): NORTA simulation with different settings.
- [2_SERGIO_Simulation.py](./2_SERGIO_Simulation.py): SERGIO simulation.


## Simulation Pre-processing

- [3_PseudoBulkData.py](./3_PseudoBulkData.py): compute pseudo-bulks of simulations.
- [3_ALRA_Imputation.R](./3_ALRA_Imputation.R): ALRA imputation.
- [3_DCA_Imputation.R](./3_DCA_Imputation.py): DCA imputation.
- [3_MAGIC_Imputation.R](./3_MAGIC_Imputation.py): MAGIC imputation.
- [3_SAVER_Imputation.R](./3_SAVER_Imputation.R): SAVER imputation.