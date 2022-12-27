suppressPackageStartupMessages({
  library(Matrix)
  library(igraph)
})

# =====================================
#           LOAD DATA
# =====================================

loadData <- function(metrics_filename, pred_filename, conf_mat_filename) {
  metric <- read.csv(metrics_filename, header = 1, row.names = 1)
  predictions <- readRDS(pred_filename)
  conf_mats <- readRDS(conf_mat_filename)
  return(list(metric = metric, predictions = predictions, conf_mats=conf_mats))
}

loadNet <- function (net_filename){
  net <- read.csv(net_filename, header = 1, row.names = 1)
  return (net)
}
# =====================================
#           VISUALIZATION
# =====================================

plotPredTaskNet <- function(true_network, est_network, name) {
  diag(true_network) <- NA
  diag(est_network) <- NA
  pdf(sprintf("./res/fig/example_network/%s-pred-network.pdf", name), width=10, height = 6)
  par(mfrow = c(1, 2), cex = 1.2)
  true_t <- which(true_network != 0.0 & lower.tri(true_network), arr.ind = TRUE)
  true_t <- cbind(true_t, true_network[which(true_network != 0.0 & lower.tri(true_network), arr.ind = TRUE)])
  append_t <- as.data.frame(do.call(rbind, lapply(1:dim(true_network)[1], function(i) { return(list(i, i, NA)) })))
  colnames(append_t) <- colnames(true_t)
  colnames(append_t)[3] <- "V3"
  true_t <- rbind(append_t, true_t) # to specify node orders in a manual way
  true_graph <- graph.data.frame(true_t, directed = F)
  E(true_graph)$color <- ifelse(E(true_graph)$V3 > 0.0, "red4", "dodgerblue4") # specify edge colors
  l <- layout.fruchterman.reingold(true_graph) # specify node layouts
  plot(true_graph, vertex.label = NA, vertex.size = 7, edge.width = 4, layout = l)
  # legend(
  #   "topleft", legend = c("positive", "negative"), col = c("red4", "dodgerblue4"),
  #   bty = "n", lty = 1, lwd = 5
  # )

  est_t <- which(est_network != 0.0 & lower.tri(est_network), arr.ind = TRUE)
  est_t <- cbind(est_t, est_network[which(est_network != 0.0 & lower.tri(est_network), arr.ind = TRUE)])
  append_t <- as.data.frame(do.call(rbind, lapply(1:dim(est_network)[1], function(i) { return(list(i, i, NA)) })))
  colnames(append_t) <- colnames(est_t)
  colnames(append_t)[3] <- "V3"
  est_t <- rbind(append_t, est_t)
  est_graph <- graph.data.frame(est_t, directed = F)
  E(est_graph)$color <- apply(est_t, MARGIN = 1, FUN = function(each) {
    if (is.na(each$V3)) {
      return("white")
    }
    if (true_network[each$row, each$col] > 0.0 & est_network[each$row, each$col] > 0.0) {
      return("red4")
    } else if (true_network[each$row, each$col] < 0.0 & est_network[each$row, each$col] < 0.0) {
      return("dodgerblue4")
    } else {
      # return("snow2")
      return("grey60")
    }
  })
  plot(est_graph, vertex.label = NA, vertex.size = 7, edge.width = 4, layout = l)
  # legend(
  #   "topright", legend = c("true positive", "true negative", "others"), col = c("red4", "dodgerblue4", "snow2"),
  #   bty = "n", lty = 1, lwd = 5)
  legend(
    "bottom", legend = c("true positive", "false positive"), col = c("red4", "grey60"),
    bty = "n", lty = 1, lwd = 5, horiz=TRUE, xjust=1.5)
  dev.off()
}

plotNetwork <- function(true_network, name) {
  diag(true_network) <- NA
  pdf(sprintf("./res/fig/example_network/%s-network.pdf", name), width=8, height = 6)
  true_t <- which(true_network != 0.0 & lower.tri(true_network), arr.ind = TRUE)
  true_t <- cbind(true_t, true_network[which(true_network != 0.0 & lower.tri(true_network), arr.ind = TRUE)])
  append_t <- as.data.frame(do.call(rbind, lapply(1:dim(true_network)[1], function(i) { return(list(i, i, NA)) })))
  colnames(append_t) <- colnames(true_t)
  colnames(append_t)[3] <- "V3"
  true_t <- rbind(append_t, true_t) # to specify node orders in a manual way
  true_graph <- graph.data.frame(true_t, directed = F)
  E(true_graph)$color <- ifelse(E(true_graph)$V3 > 0.0, "red4", "dodgerblue4") # specify edge colors
  l <- layout.fruchterman.reingold(true_graph) # specify node layouts
  plot(true_graph, vertex.label = NA, vertex.size = 7, edge.width = 4, layout = l)
  dev.off()
}

plotRefNetwork <- function(true_network, name) {
  diag(true_network) <- NA
  pdf(sprintf("./res/fig/example_network/%s-network.pdf", name), width=8, height = 6)
  true_t <- which(true_network != 0.0 & lower.tri(true_network), arr.ind = TRUE)
  true_t <- cbind(true_t, true_network[which(true_network != 0.0 & lower.tri(true_network), arr.ind = TRUE)])
  append_t <- as.data.frame(do.call(rbind, lapply(1:dim(true_network)[1], function(i) { return(list(i, i, NA)) })))
  colnames(append_t) <- colnames(true_t)
  colnames(append_t)[3] <- "V3"
  true_t <- rbind(append_t, true_t) # to specify node orders in a manual way
  true_graph <- graph.data.frame(true_t, directed = F)
  E(true_graph)$color <- ifelse(E(true_graph)$V3 > 0.0, "red4", "dodgerblue4") # specify edge colors
  # l <- layout.fruchterman.reingold(true_graph) # specify node layouts
  l <- layout.circle(true_graph)
  plot(true_graph, vertex.label = NA, vertex.size = 7, edge.width = 4, layout = l)
  dev.off()
}


plotGRNNetwork <- function(true_network, reg_names, name) {
  diag(true_network) <- NA
  pdf(sprintf("./res/fig/example_network/%s-network.pdf", name), width=8, height = 6)
  true_t <- which(true_network != 0.0 & lower.tri(true_network), arr.ind = TRUE)
  true_t <- cbind(true_t, true_network[which(true_network != 0.0 & lower.tri(true_network), arr.ind = TRUE)])
  append_t <- as.data.frame(do.call(rbind, lapply(1:dim(true_network)[1], function(i) { return(list(i, i, NA)) })))
  colnames(append_t) <- colnames(true_t)
  colnames(append_t)[3] <- "V3"
  true_t <- rbind(append_t, true_t) # to specify node orders in a manual way

  true_graph <- graph.data.frame(true_t, directed = F)
  E(true_graph)$color <- ifelse(E(true_graph)$V3 > 0.0, "red4", "dodgerblue4") # specify edge colors
  V(true_graph)$color <- "goldenrod" # specify node colors
  V(true_graph)$color[reg_names] <- "dodgerblue4"
  l <- layout.fruchterman.reingold(true_graph) # specify node layouts
  # l <- layout.circle(true_graph)
  plt <- plot(true_graph, vertex.label = NA, vertex.size = 7, edge.width = 4, layout = l)
  dev.off()
}

# =====================================

all_models_list <- c("pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet")

if(FALSE){
  data_name <- 'Cortex1-Smart_seq2-50hvg'
  simulation_name <- "NORTA"
  net_filename <- sprintf("./data/simulated/new/%s-net_mat.csv", data_name)
  true_net <- loadNet(net_filename)

  # metric_filename <- sprintf("./res/logs/evaluation/simulation/%s-%s-metrics.csv", data_name, simulation_name)
  # pred_filename <- sprintf("./res/prediction/simulation/%s-%s-predictions.rds", data_name, simulation_name)
  # conf_mat_filename <- sprintf("./res/prediction/simulation/%s-%s-sign_conf_mats.rds", data_name, simulation_name)
  metric_filename <- sprintf("./res/logs/evaluation/simulation/%s-%s-process-metrics.csv", data_name, simulation_name)
  pred_filename <- sprintf("./res/prediction/simulation/%s-%s-process-predictions.rds", data_name, simulation_name)
  conf_mat_filename <- sprintf("./res/prediction/simulation/%s-%s-process-sign_conf_mats.rds", data_name, simulation_name)
  data_obj <- loadData(metric_filename, pred_filename, conf_mat_filename)
  # -----
  model_name <- "spearman"
  model_metric <- data_obj$metric
  model_metric <- model_metric[model_metric["Model"] == model_name,]
  model_predictions <- data_obj$predictions[model_name]
  model_conf_mats <- data_obj$conf_mats[model_name]
  best_arg_idx <- as.integer(which.max(unlist(model_metric["F1"])))
  best_F1 <- model_metric[best_arg_idx, "F1"]
  print(sprintf("Best F1 = %f", best_F1))

  best_TP <- model_metric[best_arg_idx, "TP"]
  best_FP <- model_metric[best_arg_idx, "FP"]
  print(sprintf("Best TP = %f, FP=%f", best_TP, best_FP))
  # best_arg_idx <- as.integer(which.max(unlist(model_metric["Sensitivity"])))
  # best_TPR <- model_metric[best_arg_idx, "Sensitivity"]
  # print(sprintf("Best TPR = %f", best_TPR))
  # best_arg_idx <- as.integer(which.max(unlist(model_metric["Specificity"])))
  # best_FPR <- 1 - model_metric[best_arg_idx, "Specificity"]
  # print(sprintf("Best FPR = %f", best_FPR))
  best_pred <- model_predictions[[model_name]]$est_network_list[[best_arg_idx]]
  best_conf_mat <- model_conf_mats[[model_name]][[best_arg_idx]]

  plotPredTaskNet(abs(as.matrix(true_net)), abs(as.matrix(best_pred)), name=model_name)
}


# Plot simulated network
if(FALSE){
  net <- loadNet("./data/diff_graph/Cortex1-10xChromium-100hvg-GRN-net_mat.csv")
  plotNetwork(abs(net), "GRN")
  net <- loadNet("./data/diff_graph/Cortex1-10xChromium-100hvg-power-net_mat.csv")
  plotNetwork(abs(net), "hub")
  net <- loadNet("./data/diff_graph/Cortex1-10xChromium-100hvg-scale-net_mat.csv")
  plotNetwork(abs(net), "scale_free")
}


# Plot reference network
if(FALSE){
  data_name <- "Cortex1-10xChromium"
  net <- loadNet(sprintf("./data/experimental/mouse_cortex/processed/network/%s-100hvg-STRING_net.mtx", data_name))
  plotRefNetwork(abs(net), "STRING")
  net <- loadNet(sprintf("./data/experimental/mouse_cortex/processed/network/%s-100hvg-TRRUST_net.mtx", data_name))
  plotRefNetwork(abs(net), "TRRUSTv2")
}


if(FALSE){
  net <- loadNet("./data/SERGIO_simulation_all/true_network.csv")
  reg_names <- read.csv("./util/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1/Regs_cID_4.txt", header = FALSE)
  reg_names <- reg_names$V1
  plotGRNNetwork(abs(net), reg_names, "GRN")
}

