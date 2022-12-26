# Description: SAVER imputation.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

library(Matrix)
library(SAVER) # Codes are available at https://github.com/mohuangx/SAVER
source("./DataPreparation/SimulationUtils.R")

# ======================================

SAVERImpute <- function(data){
  # data: cell x gene
  imputed_data <- saver(t(data), estimates.only = TRUE)
  return (t(imputed_data))
}

# ======================================

# Imputation for NORTA simulations (mouse cortex)
if (TRUE) {
  exp_types <- c("Cortex1", "Cortex2")
  protocol_types <- c("10xChromium", "Smart_seq2")
  num_genes <- c("100hvg")
  sim_name <- "NORTA"
  data_path <- "./data/simulation/NORTA/"
  save_path <- "./data/pre-processing/imputation/NORTA"
  data_name_list <- expand.grid(exp_types, protocol_types, num_genes)
  # -----
  for (i in 1:dim(data_name_list)[1]){
    data_name <- unlist(data_name_list[i, ])
    data_name <- sprintf("%s-%s-%s", data_name["Var1"], data_name["Var2"], data_name["Var3"])
    print(data_name)
    data_mat <- loadData(sprintf("%s/%s-%s-data_mat.csv", data_path, data_name, sim_name))
     # -----
     saver_imputation <- SAVERImpute(data_mat)
     write.csv(saver_imputation, sprintf("%s/%s-%s-SAVER-data_mat.csv", save_path, data_name, sim_name))
  }
}


# Imputation for NORTA simulations (PBMC)
if (TRUE) {
  exp_types <- c("pbmc1", "pbmc2")
  protocol_types <- c("Drop", "inDrops")
  num_genes <- c("100hvg")
  sim_name <- "NORTA"
  data_path <- "./data/simulation/NORTA/"
  save_path <- "./data/pre-processing/imputation/NORTA"
  data_name_list <- expand.grid(exp_types, protocol_types, num_genes)
  # -----
  for (i in 1:dim(data_name_list)[1]){
    data_name <- unlist(data_name_list[i, ])
    data_name <- sprintf("%s-%s-%s", data_name["Var1"], data_name["Var2"], data_name["Var3"])
    print(data_name)
    data_mat <- loadData(sprintf("%s/%s-%s-data_mat.csv", data_path, data_name, sim_name))
    # -----
    saver_imputation <- SAVERImpute(data_mat)
    write.csv(saver_imputation, sprintf("%s/%s-%s-SAVER-data_mat.csv", save_path, data_name, sim_name))
  }
}


# Imputation for SERGIO simulations
if (TRUE) {
  sparsity_list <- c(1, 5, 10, 15, 20)
  data_path <- "./data/simulation/SERGIO/"
  save_path <- "./data/pre-processing/imputation/SERGIO"
  d <- "100gene-9groups"
  # -----
  for (s in sparsity_list){
    data_name <- sprintf("%s-%s", d, s)
    print(data_name)
    data_mat <- loadData(sprintf("%s/%ssparsity.csv", data_path, data_name))
    # -----
    saver_imputation <- SAVERImpute(data_mat)
    write.csv(saver_imputation, sprintf("%s/%s-SAVER-data_mat.csv", save_path, data_name))
  }
}