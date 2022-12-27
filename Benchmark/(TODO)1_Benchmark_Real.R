# Description: Benchmark of exisiting models on experimental datasets (with 1000HVGs).
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(Matrix)
  library(argparse) # argument parsing
  source("./Benchmark/BenchmarkUtils.R")
})

#TODO: network evaluation


# Test nine exisiting methods
# all_models_list <- c("PIDC", "pearson", "spearman", "glasso", "scLink", "GENIE3", "scDesign2", "minet")
all_models_list <- c("pearson", "spearman")


# Model running on experimental data
# Experimental dataset names:
#   Cortex1-10xChromium, Cortex2-10xChromium
#   Cortex2-Smart_seq2, Cortex2-Smart_seq2
#   pbmc1-Drop, pbmc2-Drop
#   pbmc1-inDrops, pbmc2-inDrops
if (TRUE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = 'pbmc1-inDrops', help = "scRNA-seq data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  args <- parser$parse_args()
  species_type <- ifelse(grepl("Cortex", data_name, fixed = TRUE), "mouse", "human")
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data %s...", data_name, species_type))
  data_dir_path <- ifelse(
    species_type == "mouse",
    "./data/experimental/mouse_cortex/processed/expr/",
    "./data/experimental/PBMC/processed/expr/"
  )
  net_dir_path <- "./data/experimental/reference_net/"

  data_filename <- sprintf("%s/%s-1000hvg.csv", data_dir_path, data_name)
  net_filename <- sprintf("%s/%s-%s.mtx", net_dir_path, data_name, true_net_name)
  data_obj <- loadData(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  print(sprintf("Data shape : #cells=%d, #genes=%d", dim(data)[1], dim(data)[2]))
  print(sprintf("#edges=%d", netEdgeNum(true_network)))
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  model_metrics <- list()
  # model_conf_mat <- list()
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
    metric_table["Network"] <- true_net_name
    metric_table["Num_Edge"] <- netEdgeNum(true_network)
    # conf_mat <- metric_summary$conf_mat
    # -----
    model_metrics[[i]] <- metric_table
    model_predictions[[model_name]] <- est_res
    # model_conf_mat[[model_name]] <- conf_mat
  }
  all_model_metrics <- reduce(model_metrics, rbind)
  options(warn = defaultW)
  # Save data
  # metric_save_path <- "./res/logs/evaluation/real/"
  # pred_save_path <- "./res/prediction/real/"
  metric_save_path <- "/gpfs/data/rsingh47/jzhan322/Coexpression_Benchmark/res/logs/evaluation/real/"
  pred_save_path <- "/gpfs/data/rsingh47/jzhan322/Coexpression_Benchmark/res/prediction/real/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-%s%s-metrics.csv", metric_save_path, data_name, true_net_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s-%s%s-predictions.rds", pred_save_path, data_name, true_net_name, need_preprocess_str))
  # saveRDS(model_conf_mat, sprintf("%s/%s-%s%s-sign_conf_mats.rds", pred_save_path, data_name, true_net_name, need_preprocess_str))
}


# Model running on Tabula Muris data
if (FALSE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = 'T_cell', help = "scRNA-seq data name.")
  parser$add_argument("-r", "--true_net", default = 'TRRUST', help = "Reference network name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  true_net_name <- args$true_net
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data %s...", data_name, true_net_name, ifelse(need_preprocess, "(pre-process)", "")))

  data_dir_path <- "/gpfs/data/rsingh47/jzhan322/Coexpression_Benchmark/data/Tabula_Muris/500hvg/"
  net_dir_path <- "/gpfs/data/rsingh47/jzhan322/Coexpression_Benchmark/data/Tabula_Muris/500hvg/"
  # data_dir_path <- "./data/Tabula_Muris/500hvg/"
  # net_dir_path <- "./data/Tabula_Muris/500hvg/"

  data_filename <- sprintf("%s/%s.csv", data_dir_path, data_name)
  net_filename <- sprintf("%s/%s-%s_net.mtx", net_dir_path, data_name, true_net_name)
  data_obj <- loadData(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  print(sprintf("Data shape : #cells=%d, #genes=%d", dim(data)[1], dim(data)[2]))
  print(sprintf("#edges=%d", netEdgeNum(true_network)))
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  # model_predictions <- list()
  model_metrics <- list()
  # model_conf_mat <- list()
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
    metric_table["Network"] <- true_net_name
    metric_table["Num_Edge"] <- netEdgeNum(true_network)
    # conf_mat <- metric_summary$conf_mat
    # -----
    model_metrics[[i]] <- metric_table
    # model_predictions[[model_name]] <- est_res
    # model_conf_mat[[model_name]] <- conf_mat
  }
  all_model_metrics <- reduce(model_metrics, rbind)
  options(warn = defaultW)
  # Save data
  metric_save_path <- "./res/logs/evaluation/real/"
  pred_save_path <- "./res/prediction/real/"
  # metric_save_path <- "/gpfs/data/rsingh47/jzhan322/Coexpression_Benchmark/res/logs/evaluation/real/"
  # pred_save_path <- "/gpfs/data/rsingh47/jzhan322/Coexpression_Benchmark/res/prediction/real/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-%s%s-metrics.csv", metric_save_path, data_name, true_net_name, need_preprocess_str))
  # saveRDS(model_predictions, sprintf("%s/%s-%s%s-predictions.rds", pred_save_path, data_name, true_net_name, need_preprocess_str))
  # saveRDS(model_conf_mat, sprintf("%s/%s-%s%s-sign_conf_mats.rds", pred_save_path, data_name, true_net_name, need_preprocess_str))
}
