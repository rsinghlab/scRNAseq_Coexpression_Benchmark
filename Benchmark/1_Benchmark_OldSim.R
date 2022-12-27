# Description: Benchmark of exisiting models on ZI-Gaussian and ZI-Poisson simulations.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(Matrix)
  library(argparse) # argument parsing
  source("./Benchmark/BenchmarkUtils.R")
})


# Test nine exisiting methods
# Note: ZILGM has a high time cost.
all_models_list <- c("pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet", "ZILGM")


# Model running on ZI-Gaussian simulations
# rho: 0.07, 0.10, 0.13, 0.16
if (TRUE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-r", "--rho", default = '0.07', help = "rho parameter value")
  parser$add_argument("-n", "--num_steps", default = '100', help = "Number of steps.")
  args <- parser$parse_args()
  data_name <- sprintf("100hvg-scLink_%srho", args$rho)
  num_steps <- as.integer(args$num_steps)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s ] Loading data...", data_name))
  data_filename <- sprintf("./data/simulation/ZI_Gaussian/%s-data_mat.csv", data_name)
  net_filename <- sprintf("./data/simulation/ZI_Gaussian/100hvg-net_mat.csv")
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
  metric_save_path <- "./res/simulation/ZI_Gaussian/"
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-metrics.csv", metric_save_path, data_name))
}


# Model running on ZI-Poisson data
# pi: 0.8, 0.9, 1.0
if (TRUE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-i", "--pi", default = '0.8', help = "pi parameter value.")
  parser$add_argument("-n", "--num_steps", default = '100', help = "Number of steps.")
  args <- parser$parse_args()
  data_name <- sprintf("100hvg-ZILGM_%spi", args$pi)
  num_steps <- as.integer(args$num_steps)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s ] Loading data...", data_name))
  data_filename <- sprintf("./data/simulation/ZI_Poisson/%s-data_mat.csv", data_name)
  net_filename <- sprintf("./data/simulation/ZI_Poisson/100hvg-net_mat.csv")
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
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
  metric_save_path <- "./res/simulation/ZI_Poisson/"
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-metrics.csv", metric_save_path, data_name))
}