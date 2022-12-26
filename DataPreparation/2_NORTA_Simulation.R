# Description: NORTA simulation processes.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(SimCorMultRes) # NORTA algorithm
  library(fitdistrplus) # distribution fitting
  library(igraph) # network visualization
  library(Matrix) # symmetric matrix
  library(Metrics)
  source("./DataPreparation/SimulationUtils.R")
})

# ======================================
#    Marginal Distribution Parameters
# ======================================

fitMarginal <- function(data, dist_name, visualize = FALSE) {
  fit_res <- fitdist(data, 'nbinom', start = list(mu = 1, size = 0.1), method = "mme")
  if (visualize) {
    par(mar = c(1, 1, 1, 1))
    plot(fit_res)
  }
  return(fit_res)
}


evaluateMarginalFitting <- function(marg_dist_list) {
  all_stats <- lapply(marg_dist_list, function(marg_dist) {
    return(gofstat(marg_dist)) # goodness-of-fit stats
  })
  return(all_stats)
}

# ======================================
#         Data Simulation
# ======================================

NORTASimulation <- function(num_cells, cov_mat, dist_list, pars_list, noise = NA, seed = 1) {
  set.seed(seed)
  simulated_data <- rnorta(
    R = num_cells,
    cor.matrix = cov_mat,
    distr = dist_list,
    qparameters = pars_list
  )
  simulated_data[which(is.na(simulated_data))] <- 0.0
  if (!is.na(noise)) {
    noise_term <- rpois(dim(simulated_data), noise)
    simulated_data <- simulated_data + noise_term
  }
  return(simulated_data)
}

# ======================================
# NOTE: The network generated from "netSimulation" function might be non-positive-definite,
# suth that the simulation may raise errors. You can re-run the simulation to solve this problem.

# Simulate with distributions learned from mouse cortex data
if (TRUE) {
  # Settings
  dist_name <- "nbinom"
  num_cells_ratio <- 1.0
  noise_loc <- NA
  # -----
  # Load mouse cortex data
  dir_path <- "./data/experimental/mouse_cortex/processed/expr/"
  exp_type <- "Cortex1" # Cortex1, Cortex2
  protocal_type <- "10xChromium" # 10xChromium, Smart_seq2
  num_genes <- "100hvg" # 100hvg, 500hvg
  print(sprintf("[ %s-%s-%s ] Data simulation", exp_type, protocal_type, num_genes))
  if (num_genes == "100hvg") {
    gene_per_clust <- 10
    num_clust <- 10
    num_hubs <- 5
    hub_degree <- 2
    num_other_edges <- 2
  } else if (num_genes == "500hvg") {
    gene_per_clust <- 20
    num_clust <- 25
    num_hubs <- 10
    hub_degree <- 4
    num_other_edges <- 4
  }
  print("Start loading data...")
  sc_data <- loadExprMat(dir_path, exp_type = exp_type, protocal_type = protocal_type, num_gene = num_genes)
  num_cells <- round(length(sc_data[, 1]) * num_cells_ratio)
  # -----
  # Generate graph
  print("Start generating network...")
  net <- netSimulation(gene_per_clust, num_clust, num_hubs, hub_degree, num_other_edges)
  net_g <- graph_from_adjacency_matrix(sign(abs(net)), mode = "undirected", diag = FALSE, weighted = TRUE)
  plot(net_g, vertex.label = NA, vertex.size = 3, edge.width = 3, main = "True Graph")
  # -----
  # Fit marginal distributions
  print("Start fitting marginal distributions...")
  fitted_pars_list <- list()
  fitted_dist_list <- list()
  for (i in 1:dim(sc_data)[2]) {
    gene_data <- as.vector(sc_data[, i])
    gene_dist <- fitMarginal(gene_data, dist_name, visualize = FALSE)
    gene_pars <- as.list(gene_dist$estimate)
    fitted_pars_list[[i]] <- list(size = gene_pars$size, prob = gene_pars$size / (gene_pars$size + gene_pars$mu))
    fitted_dist_list[[i]] <- gene_dist
  }
  # Simulate data
  print("Start generating data (NORTA)...")
  norta_simulated_data <- NORTASimulation(
    num_cells, net, rep(sprintf("q%s", dist_name), dim(sc_data)[2]), fitted_pars_list,
    noise = noise_loc, seed = 1
  )
  # Evaluation
  evaluateSim(sc_data, norta_simulated_data)
  visDataHist(sc_data, norta_simulated_data)
  # -----
  # Save data
  colnames(norta_simulated_data) <- colnames(sc_data)
  rownames(net) <- colnames(net) <- colnames(sc_data)
  rownames(norta_simulated_data) <- paste0('cell', 1:(num_cells))
  print("Start saving data...")
  file_name <- sprintf("./data/simulation/NORTA/%s-%s-%s", exp_type, protocal_type, num_genes)
  write.csv(norta_simulated_data, sprintf("%s-NORTA-data_mat.csv", file_name))
  write.csv(net, sprintf("%s-net_mat.csv", file_name))
}


# Simulate with distributions learned from PBMC data
if (TRUE) {
  # Settings
  dist_name <- "nbinom"
  num_cells_ratio <- 1.0
  noise_loc <- NA
  # -----
  # Load mouse cortex data
  dir_path <- "./data/experimental/PBMC/processed/expr/"
  exp_type <- "pbmc1" # pbmc1, pbmc2
  protocal_type <- "inDrops" # Drop, inDrops
  num_genes <- "100hvg" # 100hvg, 500hvg
  print(sprintf("[ %s-%s-%s ] Data simulation", exp_type, protocal_type, num_genes))
  if (num_genes == "100hvg") {
    gene_per_clust <- 10
    num_clust <- 10
    num_hubs <- 5
    hub_degree <- 2
    num_other_edges <- 2
  } else if (num_genes == "500hvg") {
    gene_per_clust <- 20
    num_clust <- 25
    num_hubs <- 10
    hub_degree <- 4
    num_other_edges <- 4
  }
  print("Start loading data...")
  sc_data <- loadExprMat(dir_path, exp_type = exp_type, protocal_type = protocal_type, num_gene = num_genes)
  num_cells <- round(length(sc_data[, 1]) * num_cells_ratio)
  # -----
  # Generate graph
  print("Start generating network...")
  net <- netSimulation(gene_per_clust, num_clust, num_hubs, hub_degree, num_other_edges)
  net_g <- graph_from_adjacency_matrix(sign(abs(net)), mode = "undirected", diag = FALSE, weighted = TRUE)
  plot(net_g, vertex.label = NA, vertex.size = 3, edge.width = 3, main = "True Graph")
  # -----
  # Fit marginal distributions
  print("Start fitting marginal distributions...")
  fitted_pars_list <- list()
  fitted_dist_list <- list()
  for (i in 1:dim(sc_data)[2]) {
    # for (i in 1:5) {
    gene_data <- as.vector(sc_data[, i])
    gene_dist <- fitMarginal(gene_data, dist_name, visualize = FALSE)
    gene_pars <- as.list(gene_dist$estimate)
    fitted_pars_list[[i]] <- list(size = gene_pars$size, prob = gene_pars$size / (gene_pars$size + gene_pars$mu))
    fitted_dist_list[[i]] <- gene_dist
  }
  # Simulate data
  print("Start generating data (NORTA)...")
  norta_simulated_data <- NORTASimulation(
    num_cells, net, rep(sprintf("q%s", dist_name), dim(sc_data)[2]), fitted_pars_list,
    noise = noise_loc, seed = 1
  )
  # Evaluation
  evaluateSim(sc_data, norta_simulated_data)
  visDataHist(sc_data, norta_simulated_data)
  # -----
  # Save data
  colnames(norta_simulated_data) <- colnames(sc_data)
  rownames(net) <- colnames(net) <- colnames(sc_data)
  rownames(norta_simulated_data) <- paste0('cell', 1:(num_cells))
  print("Start saving data...")
  file_name <- sprintf("./data/simulation/NORTA/%s-%s-%s", exp_type, protocal_type, num_genes)
  write.csv(norta_simulated_data, sprintf("%s-NORTA-data_mat.csv", file_name))
  write.csv(net, sprintf("%s-net_mat.csv", file_name))
}


