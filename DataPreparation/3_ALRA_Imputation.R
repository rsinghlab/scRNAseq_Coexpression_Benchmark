# Description: ALRA imputation.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

library(Matrix)
source("./DataPreparation/ALRA_src/alra.R") # Source codes of ALRA are available at: https://github.com/KlugerLab/ALRA
source("./DataPreparation/SimulationUtils.R")

# ======================================

ALRAImpute <- function(data){
  data <- as.matrix(data)
  data_norm <- normalize_data(data)
  k_choice <- choose_k(data_norm, K=dim(data_norm)[2], noise_start=5)
  k_choice <- ifelse(k_choice$k > 1, k_choice$k, 2)
  data_norm_completed <- alra(data_norm,k=k_choice)[[3]]
  return (data_norm_completed)
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
    alra_imputation <- ALRAImpute(data_mat)
    write.csv(alra_imputation, sprintf("%s/%s-%s-ALRA-data_mat.csv", save_path, data_name, sim_name))
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
    alra_imputation <- ALRAImpute(data_mat)
    write.csv(alra_imputation, sprintf("%s/%s-%s-ALRA-data_mat.csv", save_path, data_name, sim_name))
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
    alra_imputation <- ALRAImpute(data_mat)
    write.csv(alra_imputation, sprintf("%s/%s-ALRA-data_mat.csv", save_path, data_name))
  }
}