# Description: Benchmark for SOTA graphical models on simulated data.
# Author: Jiaqi Zhang

suppressPackageStartupMessages({
  library(Matrix)
  library(MASS) # for covariance estimation
  library(igraph) # network plotting
  library(cvms) # for evaluation metrics
  library(prg) # precision-recall-gain curve
  library(zoo) # rolling means
  library(argparse) # argument parsing
  library(caret) # confusion matrix
  library(purrr)
  library(ZILGM)
  library(scLink)
})


# =====================================
#           LOAD DATA
# =====================================

loadSimulation <- function(data_filename, net_filename) {
  data <- read.csv(data_filename, header = 1, row.names = 1)
  network <- read.csv(net_filename, header = 1, row.names = 1)
  return(list(data = as.matrix(data), network = as.matrix(network)))
}

# =====================================
#           GRAPHICAL MODEL
# =====================================

runZILGM <- function(data, model_args) {
  est_network <- zilgm(
    X = data, lambda = c(model_args), nlambda = 1,
    family = "NBII", update_type = "IRLS",
    # do_boot = TRUE, boot_num = 5,
    do_boot = FALSE,
    sym = "OR", verbose=2
  )
  est_network <- as.matrix(est_network$network[[1]])
  diag(est_network) <- 1.0
  est_network[which(is.na(est_network))] <- 0.0
  return(est_network)
}

# =====================================
#             EVALUATION
# =====================================

preprocessNet <- function(net) {
  net_triu <- net[upper.tri(net, diag = FALSE)]
  adj_net_triu <- 1 * (net_triu != 0)
  sign_net_triu <- sign(net_triu)
  return(list(original = net_triu, adj = adj_net_triu, sign = sign_net_triu))
}


areaUnderPRC <- function(precision, recall) {
  precision <- unlist(precision)
  recall <- unlist(recall)
  order_ind <- order(precision)
  precision <- precision[order_ind]
  recall <- recall[order_ind]
  precision <- c(precision[1], precision, 1.0)
  recall <- c(1.0, recall, 0.0)
  auprc <- sum(diff(precision) * rollmean(recall, 2))
  return(auprc)
}


areaUnderROC <- function(fpr, tpr) {
  fpr <- unlist(fpr)
  tpr <- unlist(tpr)
  order_ind <- order(fpr)
  fpr <- fpr[order_ind]
  tpr <- tpr[order_ind]
  fpr <- c(0.0, fpr, 1.0)
  tpr <- c(0.0, tpr, 1.0)
  auroc <- sum(diff(fpr) * rollmean(tpr, 2))
  return(auroc)
}


metricWithSigns <- function(est_network, true_network) {
  conf_mat <- confusionMatrix(
    factor(est_network, levels = c("-1", "0", "1")),
    factor(true_network, levels = c("-1", "0", "1")),
    mode = "everything"
  )
  kappa_val <- as.numeric(conf_mat$overall["Kappa"])
  # One vs all metrics
  cm <- conf_mat$table
  n <- sum(cm) # number of instances
  row_sums <- rowSums(cm) # number of instances per class
  col_sums <- colSums(cm) # number of predictions per class
  one_vs_all <- lapply(1:3, function(i) {
    v <- c(cm[i, i], row_sums[i] - cm[i, i], col_sums[i] - cm[i, i], n - row_sums[i] - col_sums[i] + cm[i, i])
    return(matrix(v, nrow = 2, byrow = T))
  })
  s <- Reduce("+", one_vs_all)
  TP <- s[1, 1]
  FN <- s[1, 2]
  FP <- s[2, 1]
  TN <- s[2, 2]
  avg_precision <- TP / (TP + FP)
  avg_recall <- TP / (TP + FN)
  avg_accuracy <- sum(diag(s)) / sum(s)
  avg_F1 <- (2 * avg_precision * avg_recall) / (avg_precision + avg_recall)
  return(list(
    sign_TP = TP, sign_FN = FN, sign_FP = FP, sign_TN = TN,
    sign_Precision = avg_precision, sign_Recall = avg_recall,
    sign_Accuracy = avg_accuracy, sign_F1 = avg_F1, sign_Kappa = kappa_val,
    conf_mat = cm
  ))
}


networkMetric <- function(est_network, true_network) {
  # Pre-processing
  metric_ind <- c("Balanced Accuracy", "F1", "Sensitivity", "Specificity", "Kappa", "MCC")
  true_object <- preprocessNet(true_network)
  est_object <- preprocessNet(est_network)
  # Metrics of adjacent matrix
  data_table <- as.data.frame(list(true = true_object$adj, prediction = est_object$adj))
  metrics <- evaluate(data_table, target_col = "true", prediction_cols = "prediction", type = "binomial")
  confusion_mat <- metrics$"Confusion Matrix"[[1]]
  TN <- confusion_mat$N[1]
  FP <- confusion_mat$N[2]
  FN <- confusion_mat$N[3]
  TP <- confusion_mat$N[4]
  adj_metrics <- cbind(list(TN = TN, FP = FP, FN = FN, TP = TP), metrics[metric_ind])
  # # Metrics of weight matrix
  # abs_vec_est <- abs(est_object$original)
  # normalized_vec_est <- abs_vec_est / max(abs_vec_est)
  # normalized_vec_est[is.na(normalized_vec_est)] <- 0
  # curve_metric <- t(
  #   sapply(seq(from = 0.0, to = 1.0, by = 0.01), function(thresh) {
  #     tmp_est_network <- (normalized_vec_est > thresh) * 1
  #     tmp_data_table <- as.data.frame(list(true = true_object$adj, prediction = tmp_est_network))
  #     tmp_confusion_mat <- evaluate(tmp_data_table, target_col = "true", prediction_cols = "prediction", type = "binomial")$"Confusion Matrix"[[1]]
  #     TN <- tmp_confusion_mat$N[1]
  #     FP <- tmp_confusion_mat$N[2]
  #     FN <- tmp_confusion_mat$N[3]
  #     TP <- tmp_confusion_mat$N[4]
  #     TPR <- TP / (TP + FN)
  #     FPR <- FP / (FP + TN)
  #     precision <- TP / (TP + FP)
  #     recall <- TPR
  #     return(as.data.frame(list(TPR = TPR, FPR = FPR, precision = precision, recall = recall)))
  #   })
  # )
  # curve_metric <- curve_metric[1:100,]
  # auprc <- areaUnderPRC(curve_metric[, c("precision")], curve_metric[, c("recall")])
  # auroc <- areaUnderROC(curve_metric[, c("FPR")], curve_metric[, c("TPR")])
  # prg_curve <- create_prg_curve(true_object$adj, normalized_vec_est) # precision-recall-gain-curve
  # auprg <- calc_auprg(prg_curve)
  # weight_metrics <- as.data.frame(list(AUPRG = auprg, AUPRC = auprc, AUROC = auroc))
  # # Metric considering edge signs
  # sign_metric <- metricWithSigns(est_object$sign, true_object$sign)
  # conf_mat <- sign_metric$conf_mat
  # sign_metric <- as.data.frame(list(
  #   sign_F1 = as.numeric(sign_metric$sign_F1), sign_Kappa = as.numeric(sign_metric$sign_Kappa),
  #   sign_TP = as.numeric(sign_metric$sign_TP), sign_FP = as.numeric(sign_metric$sign_FP),
  #   sign_TN = as.numeric(sign_metric$sign_TN), sign_FN = as.numeric(sign_metric$sign_FN)
  # ))
  # Summary
  # all_metrics <- rbind(c(adj_metrics, weight_metrics, sign_metric))
  all_metrics <- adj_metrics
  return(list(
    metric = all_metrics,
    # conf_mat = conf_mat,
    sparsity = mean(est_object$sign)
    # MSE = netMSE(true_object$original, est_object$original)
  ))
}


netSparsity <- function(net) {
  adj_net <- sign(abs(net))
  diag(adj_net) <- NA
  net_vec <- adj_net[which(!is.na(adj_net))]
  return(mean(net_vec))
}


netMSE <- function(true_net, est_net) {
  mse <- mean((true_net - est_net)^2) / length(true_net)
  return(mse)
}

# =====================================
#           MODEL RUNNING
# =====================================

space.sampling <- function(upper, num_steps) {
  lower <- 1e-4 * upper
  upper <- upper
  lambda_seq <- seq(lower, upper, length.out = num_steps)
  # lambda_seq <- exp(lambda_seq)
  # lambda_seq <- ifelse(lambda_seq < 10, lambda_seq, 10)
  return(lambda_seq)
}


ZILGMRunning <- function(data, num_steps) {
  lambda_max <- find_lammax(data)
  # -----
  par_space <- space.sampling(lambda_max, num_steps)
  est_network_list <- list()
  for (t in 1:length(par_space)) {
    if (t %% 5 == 0) { print(sprintf("%d iteration done", t)) }
    step_pars <- par_space[t]
    est_network_list[[t]] <- runZILGM(data, model_args = step_pars)
  }
  res <- list(est_network_list = est_network_list, par_list = par_space)
  return(res)
}


seqEstEvaluation <- function(network, est_list) {
  metric_list <- list()
  conf_mat_list <- list()
  sparsity_list <- list()
  mse_list <- list()
  for (i in 1:range(length(est_list))) {
    iter_summary <- networkMetric(est_list[[i]], network)
    metric_list[[i]] <- iter_summary$metric
    conf_mat_list[[i]] <- iter_summary$conf_mat
    sparsity_list[[i]] <- iter_summary$sparsity
    mse_list[[i]] <- iter_summary$MSE
  }
  # -----
  metric_table <- as.data.frame(reduce(metric_list, rbind))
  metric_table["Sparsity"] <- sparsity_list
  metric_table["MSE"] <- mse_list
  return(list(metric = metric_table, conf_mat = conf_mat_list))
}


# Model running on NORTA simualted data
if (TRUE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-m", "--simulation_model", default = 'NORTA', help = "Simulation method,NORTA or netSmooth.")
  # parser$add_argument("-d", "--data_name", default = 'Cortex1-10xChromium-50hvg', help = "Simulation data name.")
  parser$add_argument("-d", "--data_name", default = 'pbmc1-inDrops-100hvg', help = "Simulation data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  simulation_name <- args$simulation_model
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data %s...", data_name, simulation_name, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/simulated/new/%s-%s-data_mat.csv", data_name, simulation_name)
  net_filename <- sprintf("./data/simulated/new/%s-net_mat.csv", data_name)
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ ZILGM ] Total num of steps = %d", num_steps))
  est_res <- ZILGMRunning(data = data, num_steps = num_steps)
  # Evaluation
  metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- "ZILGM"
  metric_table["Simulation"] <- simulation_name
  model_predictions[["ZILGM"]] <- est_res
  all_model_metrics <- metric_table
  options(warn = defaultW)
  # Save data
  # a <- apply(all_model_metrics, 2, as.character)
  metric_save_path <- "./res/logs/evaluation/ZILGM_new_sim/"
  pred_save_path <- "./res/prediction/ZILGM_new_sim/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-%s%s-metrics.csv", metric_save_path, data_name, simulation_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s-%s%s-predictions.rds", pred_save_path, data_name, simulation_name, need_preprocess_str))
}


# Model running on pseudo-bulk data
if (FALSE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-m", "--simulation_model", default = 'NORTA', help = "Simulation method,NORTA or netSmooth.")
  # parser$add_argument("-d", "--data_name", default = 'Cortex1-10xChromium-50hvg', help = "Simulation data name.")
  parser$add_argument("-d", "--data_name", default = 'pbmc1-inDrops-100hvg', help = "Simulation data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  simulation_name <- args$simulation_model
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data %s...", data_name, simulation_name, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/pseudo_bulk/%s-%s-data_mat.csv", data_name, simulation_name)
  net_filename <- sprintf("./data/simulated/new/%s-net_mat.csv", data_name)
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ ZILGM ] Total num of steps = %d", num_steps))
  est_res <- ZILGMRunning(data = data, num_steps = num_steps)
  # Evaluation
  metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- "ZILGM"
  metric_table["Simulation"] <- simulation_name
  model_predictions[["ZILGM"]] <- est_res
  all_model_metrics <- metric_table
  options(warn = defaultW)
  # Save data
  # a <- apply(all_model_metrics, 2, as.character)
  metric_save_path <- "./res/logs/evaluation/ZILGM_model_pseudobulk/"
  pred_save_path <- "./res/prediction/ZILGM_model_pseudobulk/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-%s%s-metrics.csv", metric_save_path, data_name, simulation_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s-%s%s-predictions.rds", pred_save_path, data_name, simulation_name, need_preprocess_str))
}


# Model running on SERGIO simualted data
if (FALSE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = '100gene-9groups-1sparsity', help = "Simulation data name.")
  parser$add_argument("-n", "--num_steps", default = '1', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s ] Loading data %s...", data_name, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/SERGIO_simulation_all/%s.csv", data_name)
  net_filename <- sprintf("./data/SERGIO_simulation_all/true_network.csv")
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ ZILGM ] Total num of steps = %d", num_steps))
  est_res <- ZILGMRunning(data = data, num_steps = num_steps)
  # Evaluation
  metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- "ZILGM"
  model_predictions[["ZILGM"]] <- est_res
  all_model_metrics <- metric_table
  options(warn = defaultW)
  # Save data
  # a <- apply(all_model_metrics, 2, as.character)
  metric_save_path <- "./res/logs/evaluation/ZILGM_model_SERGIO/"
  pred_save_path <- "./res/prediction/ZILGM_model_SERGIO/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s%s-metrics.csv", metric_save_path, data_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s%s-predictions.rds", pred_save_path, data_name, need_preprocess_str))
}


# Model running on ZILGM simualted data
if (FALSE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = '100hvg-ZILGM_0.8pi', help = "Simulation data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s ] Loading data %s...", data_name, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/simulated/old/%s-data_mat.csv", data_name)
  net_filename <- sprintf("./data/simulated/old/100hvg-net_mat.csv")
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ ZILGM ] Total num of steps = %d", num_steps))
  est_res <- ZILGMRunning(data = data, num_steps = num_steps)
  # Evaluation
  metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- "ZILGM"
  model_predictions[["ZILGM"]] <- est_res
  all_model_metrics <- metric_table
  options(warn = defaultW)
  # Save data
  # a <- apply(all_model_metrics, 2, as.character)
  metric_save_path <- "./res/logs/evaluation/ZILGM_simulation/"
  pred_save_path <- "./res/prediction/ZILGM_simulation/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s%s-ZILGM-metrics.csv", metric_save_path, data_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s%s-ZILGM-predictions.rds", pred_save_path, data_name, need_preprocess_str))
}


# Model running on scLink simualted data
if (FALSE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = '50hvg-scLink_0.07rho', help = "Simulation data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s ] Loading data %s...", data_name, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/simulated/old/%s-data_mat.csv", data_name)
  net_filename <- sprintf("./data/simulated/old/50hvg-net_mat.csv")
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ ZILGM ] Total num of steps = %d", num_steps))
  est_res <- ZILGMRunning(data = data, num_steps = num_steps)
  # Evaluation
  metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- "ZILGM"
  model_predictions[["ZILGM"]] <- est_res
  all_model_metrics <- metric_table
  options(warn = defaultW)
  # Save data
  # a <- apply(all_model_metrics, 2, as.character)
  metric_save_path <- "./res/logs/evaluation/scLink_simulation/"
  pred_save_path <- "./res/prediction/scLink_simulation/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s%s-ZILGM-metrics.csv", metric_save_path, data_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s%s-ZILGM-predictions.rds", pred_save_path, data_name, need_preprocess_str))
}


# Model running on SERGIO pseudo-bulk data
if (FALSE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-m", "--simulation_model", default = 'SERGIO', help = "Simulation method,NORTA or netSmooth.")
  parser$add_argument("-d", "--data_name", default = '100gene-9groups-20', help = "Simulation data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  simulation_name <- args$simulation_model
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data %s...", data_name, simulation_name, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/SERGIO_pseudobulk/%s.csv", data_name)
  net_filename <- sprintf("./data/SERGIO_simulation_all/true_network.csv")
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ ZILGM ] Total num of steps = %d", num_steps))
  est_res <- ZILGMRunning(data = data, num_steps = num_steps)
  # Evaluation
  metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- "ZILGM"
  metric_table["Simulation"] <- simulation_name
  model_predictions[["ZILGM"]] <- est_res
  all_model_metrics <- metric_table
  options(warn = defaultW)
  # Save data
  # a <- apply(all_model_metrics, 2, as.character)
  metric_save_path <- "./res/logs/evaluation/ZILGM_model_SERGIO_pseudobulk/"
  pred_save_path <- "./res/prediction/ZILGM_model_SERGIO_pseudobulk/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s%s-metrics.csv", metric_save_path, data_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s%s-predictions.rds", pred_save_path, data_name, need_preprocess_str))
}

# --------------------------------
# Model running on different cells
if (FALSE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-m", "--simulation_model", default = 'NORTA', help = "Simulation method,NORTA or netSmooth.")
  parser$add_argument("-d", "--data_name", default = 'Cortex1-10xChromium-100hvg-2.0cell', help = "Simulation data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  simulation_name <- args$simulation_model
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data %s...", data_name, simulation_name, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/diff_cell/%s-%s-data_mat.csv", data_name, simulation_name)
  net_filename <- sprintf("./data/diff_cell/Cortex1-10xChromium-100hvg-net_mat.csv")
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ ZILGM ] Total num of steps = %d", num_steps))
  est_res <- ZILGMRunning(data = data, num_steps = num_steps)
  # Evaluation
  metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- "ZILGM"
  metric_table["Simulation"] <- simulation_name
  model_predictions[["ZILGM"]] <- est_res
  all_model_metrics <- metric_table
  options(warn = defaultW)
  # Save data
  # a <- apply(all_model_metrics, 2, as.character)
  metric_save_path <- "./res/logs/evaluation/diff_cell/"
  pred_save_path <- "./res/prediction/diff_cell/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s%s-ZILGM-metrics.csv", metric_save_path, data_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s%s-ZILGM-predictions.rds", pred_save_path, data_name, need_preprocess_str))
}


# Model running on different graph
if (FALSE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-m", "--simulation_model", default = 'NORTA', help = "Simulation method,NORTA or netSmooth.")
  parser$add_argument("-d", "--data_name", default = 'Cortex1-10xChromium-100hvg-scale', help = "Simulation data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  simulation_name <- args$simulation_model
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data %s...", data_name, simulation_name, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/diff_graph/%s-%s-data_mat.csv", data_name, simulation_name)
  net_filename <- sprintf("./data/diff_graph/%s-net_mat.csv", data_name)
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ ZILGM ] Total num of steps = %d", num_steps))
  est_res <- ZILGMRunning(data = data, num_steps = num_steps)
  # Evaluation
  metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- "ZILGM"
  metric_table["Simulation"] <- simulation_name
  model_predictions[["ZILGM"]] <- est_res
  all_model_metrics <- metric_table
  options(warn = defaultW)
  # Save data
  # a <- apply(all_model_metrics, 2, as.character)
  metric_save_path <- "./res/logs/evaluation/diff_graph/"
  pred_save_path <- "./res/prediction/diff_graph/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s%s-ZILGM-metrics.csv", metric_save_path, data_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s%s-ZILGM-predictions.rds", pred_save_path, data_name, need_preprocess_str))
}


# Model running on different sparsity
if (FALSE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-m", "--simulation_model", default = 'NORTA', help = "Simulation method,NORTA or netSmooth.")
  parser$add_argument("-d", "--data_name", default = '100gene-hub-full', help = "Simulation data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  simulation_name <- args$simulation_model
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data %s...", data_name, simulation_name, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/diff_sparsity/%s-%s-data_mat.csv", data_name, simulation_name)
  net_filename <- sprintf("./data/diff_sparsity/100gene-hub-net_mat.csv")
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ ZILGM ] Total num of steps = %d", num_steps))
  est_res <- ZILGMRunning(data = data, num_steps = num_steps)
  # Evaluation
  metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- "ZILGM"
  metric_table["Simulation"] <- simulation_name
  model_predictions[["ZILGM"]] <- est_res
  all_model_metrics <- metric_table
  options(warn = defaultW)
  # Save data
  # a <- apply(all_model_metrics, 2, as.character)
  metric_save_path <- "./res/logs/evaluation/diff_sparsity/"
  pred_save_path <- "./res/prediction/diff_sparsity/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s%s-ZILGM-metrics.csv", metric_save_path, data_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s%s-ZILGM-predictions.rds", pred_save_path, data_name, need_preprocess_str))
}


# Model running on quality control
if (FALSE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-m", "--simulation_model", default = 'NORTA-1percentile', help = "Simulation method,NORTA or netSmooth.")
  parser$add_argument("-d", "--data_name", default = 'Cortex1-Smart_seq2-100hvg', help = "Simulation data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  simulation_name <- args$simulation_model
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data %s...", data_name, simulation_name, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/quality_control/%s-%s-data_mat.csv", data_name, simulation_name)
  net_filename <- sprintf("./data/simulated/new/%s-net_mat.csv", data_name)
  data_obj <- loadSimulation(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  # Pre-processing
  if (need_preprocess) {
    data <- sclink_norm(data, scale.factor = 1e4, gene.names = colnames(data))
    data[which(is.na(data))] <- 0.0
  }
  # Model running
  model_predictions <- list()
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ ZILGM ] Total num of steps = %d", num_steps))
  est_res <- ZILGMRunning(data = data, num_steps = num_steps)
  # Evaluation
  metric_summary <- seqEstEvaluation(true_network, est_res$est_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- "ZILGM"
  metric_table["Simulation"] <- simulation_name
  model_predictions[["ZILGM"]] <- est_res
  all_model_metrics <- metric_table
  options(warn = defaultW)
  # Save data
  # a <- apply(all_model_metrics, 2, as.character)
  metric_save_path <- "./res/logs/evaluation/ZILGM_model_qc/"
  pred_save_path <- "./res/prediction/ZILGM_model_qc/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-%s%s-metrics.csv", metric_save_path, data_name, simulation_name, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s-%s%s-predictions.rds", pred_save_path, data_name, simulation_name, need_preprocess_str))
}