# Description: Benchmark of exisiting models on simulations with normalization.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(Matrix)
  library(argparse) # argument parsing

  library(sctransform) # for sctransform
  library(PsiNorm) # for psiNorm

  source("./Benchmark/BenchmarkUtils.R")
})


# =====================================
#           NORMALIZATION
# =====================================

normData <- function(data, type){
  if (type == "lognorm"){
    norm_data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data)) # this is equivalent to lognorm essentially
    norm_data[which(is.na(norm_data))] <- 0.0
  } else if (type == "sctransform"){
    gene_by_cell <- t(data)
    cell_sparsity <- colSums(gene_by_cell)
    non_zero_idx <- which(cell_sparsity != 0)
    norm_mat <- matrix(0.0, dim(gene_by_cell)[1], dim(gene_by_cell)[2])
    row.names(norm_mat) <- row.names(gene_by_cell)
    colnames(norm_mat) <- colnames(gene_by_cell)
    tmp_data <- gene_by_cell[,non_zero_idx]
    res <- sctransform::vst(tmp_data)$y
    res_rows <- rownames(res)
    res_cols <- colnames(res)
    norm_mat[res_rows,res_cols] <- res
    norm_data <- t(norm_mat)
  } else if (type == "psiNorm"){
    res <- PsiNorm(t(data)) # requires gene x cell
    res[which(is.na(res))] <- 0.0
    norm_data <- t(res)
  } else {
    stop(sprintf("Unknown normalziation type: %s", type))
  }
  return (norm_data)
}

# ----------------------------------------

# Test nine exisiting methods
# Note: ZILGM has a high time cost.
all_models_list <- c("pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet", "ZILGM")


# ----------------------------------------
# Three normalization types: lognorm, sctransform, and psiNorm


# Model running on eight NORTA simulation datasets (100 genes)
# Simulation dataset names:
#   Cortex1-10xChromium-100hvg, Cortex2-10xChromium-100hvg
#   Cortex2-Smart_seq2-100hvg, Cortex2-Smart_seq2-100hvg
#   pbmc1-Drop-100hvg, pbmc2-Drop-100hvg
#   pbmc1-inDrops-100hvg, pbmc2-inDrops-100hvg
if (TRUE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = 'Cortex1-10xChromium-100hvg', help = "scRNA-seq data name.")
  parser$add_argument("-n", "--num_steps", default = '100', help = "Number of steps.")
  parser$add_argument("-p", "--norm_type", default = "lognorm", help = "Type of normalization (lognorm, sctransform, or psiNorm).")
  args <- parser$parse_args()
  data_name <- args$data_name
  simulation_name <- "NORTA"
  num_steps <- as.integer(args$num_steps)
  norm_type <- args$norm_type
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data...", data_name, simulation_name))
  data_filename <- sprintf("./data/simulation/NORTA/%s-%s-data_mat.csv", data_name, simulation_name)
  net_filename <- sprintf("./data/simulation/NORTA/%s-net_mat.csv", data_name)
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  print(sprintf("Normalizing (%s)...", norm_type))
  data <- normData(data, norm_type)
  # Model running
  model_metrics <- list()
  for (i in 1:range(length(all_models_list))) {
    model_name <- all_models_list[i]
    print(paste(rep("=", 70), collapse = ""))
    print(sprintf("[ %s ] Total num of steps = %d", model_name, num_steps))
    est_res <- runModel(name = model_name, data = data, num_steps = num_steps)
    # Evaluation
    metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
    metric_table <- metric_summary$metric
    metric_table["Pars"] <- unlist(est_res$par_list)
    metric_table["Model"] <- model_name
    metric_table["Simulation"] <- simulation_name
    # -----
    model_metrics[[i]] <- metric_table
  }
  all_model_metrics <- do.call("rbind", model_metrics)
  options(warn = defaultW)
  # Save data
  metric_save_path <- "./res/pre-processing/normalization/NORTA/"
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-%s-metrics.csv", metric_save_path, data_name, norm_type))
}


# Model running on SERGIO simualted data (100 genes)
# Simulation dataset names:
#   100gene-9groups-1, 100gene-9groups-5
#   100gene-9groups-10, 100gene-9groups-15
#   100gene-9groups-20
if (TRUE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = '100gene-9groups-1', help = "scRNA-seq data name.")
  parser$add_argument("-n", "--num_steps", default = '100', help = "Number of steps.")
  parser$add_argument("-p", "--norm_type", default = "lognorm", help = "Type of normalization (lognorm, sctransform, or psiNorm).")
  args <- parser$parse_args()
  data_name <- args$data_name
  simulation_name <- "NORTA"
  num_steps <- as.integer(args$num_steps)
  norm_type <- args$norm_type
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data...", data_name, simulation_name))
  data_filename <- sprintf("./data/simulation/SERGIO/%ssparsity.csv", data_name)
  net_filename <- sprintf("./data/simulation/SERGIO/true_network.csv")
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  print(sprintf("Normalizing (%s)...", norm_type))
  data <- normData(data, norm_type)
  # Model running
  model_metrics <- list()
  for (i in 1:range(length(all_models_list))) {
    model_name <- all_models_list[i]
    print(paste(rep("=", 70), collapse = ""))
    print(sprintf("[ %s ] Total num of steps = %d", model_name, num_steps))
    est_res <- runModel(name = model_name, data = data, num_steps = num_steps)
    # Evaluation
    metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
    metric_table <- metric_summary$metric
    metric_table["Pars"] <- unlist(est_res$par_list)
    metric_table["Model"] <- model_name
    metric_table["Simulation"] <- simulation_name
    # -----
    model_metrics[[i]] <- metric_table
  }
  all_model_metrics <- do.call("rbind", model_metrics)
  options(warn = defaultW)
  # Save data
  metric_save_path <- "./res/pre-processing/normalization/SERGIO/"
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-%s-metrics.csv", metric_save_path, data_name, norm_type))
}

