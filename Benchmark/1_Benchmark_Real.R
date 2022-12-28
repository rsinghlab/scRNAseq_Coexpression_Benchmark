# Description: Benchmark of exisiting models on experimental datasets (with 1000HVGs).
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(Matrix)
  library(argparse) # argument parsing
  library(stringr)
  source("./Benchmark/BenchmarkUtils.R")
})

loadRealData <- function(data_filename, net_filename) {
  data <- read.csv(data_filename, header = 1, row.names = 1)
  gene_names <- colnames(data)
  gene_names <- str_split(gene_names, "_")
  gene_names <- unlist(lapply(gene_names, function(x){return (x[[2]])}))
  colnames(data) <- gene_names
  network <- read.csv(net_filename, header = 1, row.names = 1)
  return(list(data = as.matrix(data), network = as.matrix(network)))
}

# Test nine exisiting methods
# all_models_list <- c("PIDC", "pearson", "spearman", "glasso", "scLink", "GENIE3", "scDesign2", "minet")
all_models_list <- c("pearson")


# Model running on real mouse cortex data
# Experimental dataset names:
#   Cortex1-10xChromium, Cortex2-10xChromium
#   Cortex2-Smart_seq2, Cortex2-Smart_seq2
if (TRUE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = 'Cortex1-10xChromium', help = "scRNA-seq data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  args <- parser$parse_args()
  species_type <- "mouse"
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data...", data_name, species_type))
  data_dir_path <- "./data/experimental/mouse_cortex/processed/expr/"
  net_dir_path <- "./data/experimental/reference_net/"
  data_filename <- sprintf("%s/%s-1000hvg.csv", data_dir_path, data_name)
  net_filename <- sprintf("%s/%s-mouse_TF_PCC-sub_net_mat.csv", net_dir_path, data_name)
  data_obj <- loadRealData(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  TF_list <- rownames(true_network)
  print(sprintf("Data shape : #cells=%d, #genes=%d", dim(data)[1], dim(data)[2]))
  # Model running
  model_metrics <- list()
  for (i in 1:range(length(all_models_list))) {
    model_name <- all_models_list[i]
    print(paste(rep("=", 70), collapse = ""))
    print(sprintf("[ %s ] Total num of steps = %d", model_name, num_steps))
    est_res <- runModel(name = model_name, data = data, num_steps = num_steps)
    # Evaluation
    est_sub_network_list <- lapply(est_res$est_network_list, function (net){
      return (net[TF_list, TF_list])
    })
    metric_summary <- seqEstEvaluation(true_network, est_sub_network_list)
    metric_table <- metric_summary$metric
    metric_table["Pars"] <- unlist(est_res$par_list)
    metric_table["Model"] <- model_name
    # -----
    model_metrics[[i]] <- metric_table
  }
  all_model_metrics <- reduce(model_metrics, rbind)
  options(warn = defaultW)
  # Save data
  metric_save_path <- "./res/experimental/mouse_cortex/"
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-metrics.csv", metric_save_path, data_name))
}


# Model running on real PBMC data
# Experimental dataset names:
#   pbmc1-Drop, pbmc2-Drop
#   pbmc1-inDrops, pbmc2-inDrops
if (TRUE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = 'pbmc1-Drop', help = "scRNA-seq data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  args <- parser$parse_args()
  species_type <- "human"
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data...", data_name, species_type))
  data_dir_path <- "./data/experimental/PBMC/processed/expr/"
  net_dir_path <- "./data/experimental/reference_net/"
  data_filename <- sprintf("%s/%s-1000hvg.csv", data_dir_path, data_name)
  net_filename <- sprintf("%s/%s-human_TF_similarity-sub_net_mat.csv", net_dir_path, data_name)
  data_obj <- loadRealData(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  TF_list <- rownames(true_network)
  print(sprintf("Data shape : #cells=%d, #genes=%d", dim(data)[1], dim(data)[2]))
  # Model running
  model_metrics <- list()
  for (i in 1:range(length(all_models_list))) {
    model_name <- all_models_list[i]
    print(paste(rep("=", 70), collapse = ""))
    print(sprintf("[ %s ] Total num of steps = %d", model_name, num_steps))
    est_res <- runModel(name = model_name, data = data, num_steps = num_steps)
    # Evaluation
    idx_list <- intersect(TF_list, colnames(data))
    est_sub_network_list <- lapply(est_res$est_network_list, function (net){
      return (net[idx_list, idx_list])
    })
    true_sub_network <- true_network[idx_list, idx_list]
    metric_summary <- seqEstEvaluation(true_sub_network, est_sub_network_list)
    metric_table <- metric_summary$metric
    metric_table["Pars"] <- unlist(est_res$par_list)
    metric_table["Model"] <- model_name
    # -----
    model_metrics[[i]] <- metric_table
  }
  all_model_metrics <- reduce(model_metrics, rbind)
  options(warn = defaultW)
  # Save data
  metric_save_path <- "./res/experimental/PBMC/"
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-metrics.csv", metric_save_path, data_name))
}