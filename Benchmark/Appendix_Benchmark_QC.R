# Description: Benchmark of exisiting models on simulations after removing lowly expressed cells.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(Matrix)
  library(argparse) # argument parsing
  source("./Benchmark/BenchmarkUtils.R")
})

# Test nine exisiting methods
# Note: ZILGM has a high time cost.
all_models_list <- c("pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet", "ZILGM")


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
  parser$add_argument("-d", "--data_name", default = 'pbmc1-Drop-100hvg', help = "scRNA-seq data name.")
  parser$add_argument("-r", "--qc", default = '1', help = "Quality control parameter. 1, 5, 10, 15, 20")
  parser$add_argument("-n", "--num_steps", default = '100', help = "Number of steps.")
  args <- parser$parse_args()
  qc_percent <- args$qc
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  sim_name <- "NORTA"
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s percent ] Loading data...", data_name, qc_percent))
  data_filename <- sprintf("./data/appendix/quality_control/%s-%s-%spercentile-data_mat.csv", data_name, sim_name, qc_percent)
  net_filename <- sprintf("./data/simulation/NORTA/%s-net_mat.csv", data_name)
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
    metric_table["QC_Percent"] <- qc_percent
    # -----
    model_metrics[[i]] <- metric_table
  }
  all_model_metrics <- reduce(model_metrics, rbind)
  options(warn = defaultW)
  # Save data
  metric_save_path <- "./res/appendix/quality_control/"
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-%s-%spercentile-metrics.csv", metric_save_path, data_name, sim_name, qc_percent))
}



