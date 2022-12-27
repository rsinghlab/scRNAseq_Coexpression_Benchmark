# Description: Benchmark of exisiting models on NORTA simulations with diff sparsity levels.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(Matrix)
  library(argparse) # argument parsing
  source("./Benchmark/BenchmarkUtils.R")
})


# Test nine exisiting methods
# Note: ZILGM has a high time cost.
all_models_list <- c("pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet", "ZILGM")


# Model running on NORTA simulations with different number of cells
# Sparsity level: 0.1, 0.2, ..., 0.9
if (TRUE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = '100gene-hub', help = "scRNA-seq data name.")
  parser$add_argument("-s", "--sparsity", default = '0.9', help = "Sparsity level.")
  parser$add_argument("-n", "--num_steps", default = '100', help = "Number of steps.")
  args <- parser$parse_args()
  data_name <- args$data_name
  sparsity_level <- args$sparsity
  num_steps <- as.integer(args$num_steps)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data...", data_name, sparsity_level))
  data_filename <- sprintf("./data/diff_settings/diff_sparsity/%s-%ssparsity-NORTA-data_mat.csv", data_name, sparsity_level)
  net_filename <- sprintf("./data/diff_settings/diff_sparsity/100gene-hub-net_mat.csv")
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  print(sprintf("Data shape : #cells=%d, #genes=%d", dim(data)[1], dim(data)[2]))
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
    # -----
    model_metrics[[i]] <- metric_table
  }
  all_model_metrics <- reduce(model_metrics, rbind)
  options(warn = defaultW)
  # Save data
  metric_save_path <- "./res/diff_settings/diff_sparsity/"
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-%ssparsity-metrics.csv", metric_save_path, data_name, sparsity_level))
}