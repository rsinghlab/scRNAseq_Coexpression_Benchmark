# Description: Benchmark for SOTA graphical models on filtered (ie. removing cells with low counts) simulation data.
# Author: Jiaqi Zhang

suppressPackageStartupMessages({
  library(Matrix)
  library(MASS) # for covariance estimation
  library(QUIC) # QUIC model
  library(igraph) # network plotting
  library(cvms) # for evaluation metrics
  library(prg) # precision-recall-gain curve
  library(zoo) # rolling means
  library(argparse) # argument parsing
  library(scLink) # for normalization and covaraince estimation
  library(caret) # confusion matrix
  library(GENIE3)
  # library(WGCNA)
  library(PIDC)
  library(minet)
  files.sources <- list.files("./util/scDesign2/R/", full.names = TRUE)
  sapply(files.sources, source)
})


# =====================================
#           LOAD DATA
# =====================================

loadData <- function(data_filename, net_filename) {
  data <- read.csv(data_filename, header = 1, row.names = 1)
  network <- read.csv(net_filename, header = 1, row.names = 1)
  return(list(data = as.matrix(data), network = as.matrix(network)))
}

# =====================================
#        COVARIANCE COMPUTATION
# =====================================

computeCov <- function(data, type) {
  if (type == "pearson") {
    cov_mat <- cor(data, method = "pearson")
  } else if (type == "spearman") {
    cov_mat <- cor(data, method = "spearman")
  } else if (type == "covariance") {
    data <- apply(data, MARGIN = 2, FUN = function(x) { return(x - mean(x)) })
    cov_mat <- cov(data)
  } else if (type == "scLink") {
    # ========================================
    cov_mat <- sclink_cor(data, ncores = 1, nthre = 10, dthre = 0.9)
  } else if (type == "scDesign2") {
    cov_mat <- fit_Gaussian_copula(t(data), jitter = FALSE, zp_cutoff = 1.1)$cov_mat # requires gene x cell
  } else if (type == "GENIE3") {
    cov_mat <- GENIE3(t(data)) # requires gene x cell
  } else if (type == "PIDC") {
    cov_mat <- PIDC(t(data), verbose = FALSE) # requires gene x cell
  } else if (type == "minet") {
    cov_mat <- minet(data) # requires cell x gene; default parameters method="mrnet", estimator="spearman"
  } else {
    stop(sprintf("Unkown covariance type %s!", type))
  }
  cov_mat[which(is.na(cov_mat))] <- 0.0
  return(cov_mat)
}

# =====================================
#           GRAPHICAL MODEL
# =====================================

thresholdingMat <- function(mat, thr) {
  # mat_dim <- dim(mat)
  # for (i in 1:(mat_dim[1] - 1)) {
  #   for (j in (i + 1):mat_dim[2]) {
  #     if (abs(mat[i, j]) < thr) {
  #       mat[i, j] <- 0.0
  #       mat[j, i] <- 0.0
  #     }
  #   }
  # }
  # return(mat)
  diag_val <- diag(mat)
  mat <- as.matrix(forceSymmetric(mat))
  mat[abs(mat) < thr] <- 0.0
  diag(mat) <- diag_val
  return(mat)
}


estimateNet <- function(cov_mat, model_args, type) {
  if (type == "thresholding") {
    model_args$mat <- cov_mat
    est_network <- do.call("thresholdingMat", model_args)
  } else if (type == "glasso") {
    model_args$S <- cov_mat
    model_args$msg <- 0
    est_network <- do.call("QUIC", model_args)$X
  } else {
    stop(sprintf("Unknown graphical model %s!", type))
  }
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
  # Metric considering edge signs
  # sign_metric <- metricWithSigns(est_object$sign, true_object$sign)
  # conf_mat <- sign_metric$conf_mat
  # sign_metric <- as.data.frame(list(
  #   sign_F1 = as.numeric(sign_metric$sign_F1), sign_Kappa = as.numeric(sign_metric$sign_Kappa),
  #   sign_TP = as.numeric(sign_metric$sign_TP), sign_FP = as.numeric(sign_metric$sign_FP),
  #   sign_TN = as.numeric(sign_metric$sign_TN), sign_FN = as.numeric(sign_metric$sign_FN)
  # ))
  # Summary
  # all_metrics <- rbind(c(adj_metrics, weight_metrics, sign_metric))
  # all_metrics <- rbind(c(adj_metrics, sign_metric))
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


netEdgeNum <- function (network){
  network[which(is.na(network))] <- 0.0
  network <- sign(abs(network))
  triu_net <- network[upper.tri(network, diag = FALSE)]
  edge_num <- sum(triu_net)
  return (edge_num)
}

# =====================================
#           MODEL RUNNING
# =====================================

space.sampling <- function(model_name, upper, num_steps) {
  if (model_name == "thresholding") {
    lower <- 1e-3
    upper <- upper - 1e-3
    step <- (upper - lower) / num_steps
    thr_seq <- seq(lower, upper, by = step)
    random_thr_seq <- sample(num_steps)
    model_args <- lapply(1:num_steps, function(i) { return(list(thr = thr_seq[random_thr_seq[i]])) })
  } else if (model_name == "glasso" | model_name == "quic") {
    lower <- 0.001
    upper <- 1.0
    step <- (upper - lower) / num_steps
    rho_seq <- seq(lower, upper, by = step)
    random_rho_seq <- sample(num_steps)
    model_args <- lapply(1:num_steps, function(i) { return(list(rho = rho_seq[random_rho_seq[i]])) })
  } else {
    stop(sprintf("Unknown model name %s!", model_name))
  }
  return(model_args)
}


runModel <- function(name, data, num_steps) {
  # pearson, spearman, glasso, scLink, GENIE3, PIDC, scDesign2, minet
  if (name == "pearson") {
    cov_type <- "pearson"
    estimate_type <- "thresholding"
  } else if (name == "spearman") {
    cov_type <- "spearman"
    estimate_type <- "thresholding"
  } else if (name == "glasso") {
    cov_type <- "covariance"
    estimate_type <- "glasso"
  } else if (name == "scLink") {
    cov_type <- "scLink"
    estimate_type <- "glasso"
  } else if (name == "GENIE3") {
    cov_type <- "GENIE3"
    estimate_type <- "thresholding"
  } else if (name == "PIDC") {
    cov_type <- "PIDC"
    estimate_type <- "thresholding"
  } else if (name == "scDesign2") {
    cov_type <- "scDesign2"
    estimate_type <- "glasso"
  } else if (name == "minet") {
    cov_type <- "minet"
    estimate_type <- "thresholding"
  }else {
    stop(sprintf("Unknown model name %s!", name))
  }
  # -----
  cov_mat <- computeCov(data, type = cov_type)
  par_space <- space.sampling(estimate_type, max(abs(cov_mat)), num_steps)
  est_network_list <- list()
  for (t in 1:length(par_space)) {
    if (t %% 5 == 0) { print(sprintf("%d iteration done", t)) }
    step_pars <- par_space[[t]]
    est_network_list[[t]] <- estimateNet(cov_mat, model_args = step_pars, type = estimate_type)
  }
  res <- list(name = name, cov_type = cov_type, estimate_type = estimate_type, est_network_list = est_network_list, par_list = par_space)
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

all_models_list <- c("pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet")
# all_models_list <- c("pearson", "spearman")

# Model running on filtered simulated data
if (TRUE) {
  defaultW <- getOption("warn")
  options(warn = -1)
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = 'pbmc1-Drop-50hvg', help = "scRNA-seq data name.")
  parser$add_argument("-r", "--qc", default = '5', help = "Imputaion type.")
  parser$add_argument("-n", "--num_steps", default = '2', help = "Number of steps.")
  parser$add_argument("-p", "--need_preprocessing", default = "0", help = "Need preprocessing.")
  args <- parser$parse_args()
  qc_percent <- args$qc
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  sim_name <- "NORTA"
  need_preprocess <- ifelse(args$need_preprocessing == "1", TRUE, FALSE)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s percent ] Loading data %s...", data_name, qc_percent, ifelse(need_preprocess, "(pre-process)", "")))
  data_filename <- sprintf("./data/quality_control/%s-%s-%spercentile-data_mat.csv", data_name, sim_name, qc_percent)
  net_filename <- sprintf("./data/simulated/new/%s-net_mat.csv", data_name)
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
    metric_table["Num_Edge"] <- netEdgeNum(true_network)
    # -----
    model_metrics[[i]] <- metric_table
    model_predictions[[model_name]] <- est_res
  }
  all_model_metrics <- reduce(model_metrics, rbind)
  options(warn = defaultW)
  # Save data
  metric_save_path <- "./res/logs/evaluation/quality_control/"
  pred_save_path <- "./res/prediction/quality_control/"
  need_preprocess_str <- ifelse(need_preprocess, "-process", "")
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-%s-%spercentile%s-metrics.csv", metric_save_path, data_name, sim_name, qc_percent, need_preprocess_str))
  saveRDS(model_predictions, sprintf("%s/%s-%s-%spercentile%s-predictions.rds", pred_save_path, data_name, sim_name, qc_percent, need_preprocess_str))
}



