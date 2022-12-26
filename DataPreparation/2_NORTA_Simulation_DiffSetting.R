# Description: NORTA simulations with different settings.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(SimCorMultRes) # NORTA algorithm
  library(fitdistrplus) # distribution fitting
  library(igraph) # network visualization
  library(Matrix) # symmetric matrix
  library(matrixcalc)
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


addDropouts <- function(cnt, s) {
  shape <- 1
  cnt_zero <- cnt
  log_cnt_zero <- log(cnt_zero + 1)
  log_mid_point <- as.numeric(quantile(log_cnt_zero, s))
  prob_b <- 1.0 / (1 + exp(-1 * shape * (log_cnt_zero - log_mid_point)))
  dropind <- rbinom(prod(dim(cnt_zero)), size = 1, prob = as.vector(prob_b))
  cnt_zero[dropind == 1] <- 0.0
  return(cnt_zero)
}

# ======================================

# Simulate with pars learned from mouse cortex data (different cell number)
if (TRUE) {
  # Settings
  dist_name <- "nbinom"
  noise_loc <- NA
  # -----
  # Load mouse cortex data
  dir_path <- "./data/experimental/mouse_cortex/processed/expr/"
  exp_type <- "Cortex1" # Cortex1, Cortex2
  protocal_type <- "10xChromium" # 10xChromium, Smart_seq2
  num_genes <- "100hvg"
  file_name <- sprintf("./data/diff_settings/diff_cell/%s-%s-%s", exp_type, protocal_type, num_genes)
  gene_per_clust <- 10
  num_clust <- 10
  num_hubs <- 5
  hub_degree <- 2
  num_other_edges <- 2
  print("Start loading data...")
  sc_data <- loadExprMat(dir_path, exp_type = exp_type, protocal_type = protocal_type, num_gene = num_genes)
  # -----
  # Generate graph
  print("Start generating network...")
  net <- netSimulation(gene_per_clust, num_clust, num_hubs, hub_degree, num_other_edges)
  net_g <- graph_from_adjacency_matrix(sign(abs(net)), mode = "undirected", diag = FALSE, weighted = TRUE)
  plot(net_g, vertex.label = NA, vertex.size = 3, edge.width = 3, main = "True Graph")
  write.csv(net, sprintf("%s-net_mat.csv", file_name))
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
  for (n in c(0.1, 0.5, 1.0, 1.5, 2.0)) {
    print(paste(rep("=", 70), collapse = ""))
    num_cells <- round(length(sc_data[1,]) * n)
    print(sprintf("[ Cell ratio = %f | # cells = %d] ...", n, num_cells))
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
    write.csv(norta_simulated_data, sprintf("%s-%.1fcell-NORTA-data_mat.csv", file_name, n))
  }
}


# Simulate with pars learned from mouse cortex data (different graph structure)
if (TRUE) {
  # Settings
  dist_name <- "nbinom"
  noise_loc <- NA
  # -----
  # Load mouse cortex data
  dir_path <- "./data/experimental/mouse_cortex/processed/expr/"
  exp_type <- "Cortex1" # Cortex1, Cortex2
  protocal_type <- "10xChromium" # 10xChromium, Smart_seq2
  num_genes <- "100hvg"
  file_name <- sprintf("./data/diff_settings/diff_graph/%s-%s-%s", exp_type, protocal_type, num_genes)
  print("Start loading data...")
  sc_data <- loadExprMat(dir_path, exp_type = exp_type, protocal_type = protocal_type, num_gene = num_genes)
  num_cells <- round(length(sc_data[1,]) * 1.0)
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
  # hub: hub-based network
  # power: scale-free network
  # GRN: gene regulatory network
  for (g_t in c("hub", "power", "GRN")) {
    print(paste(rep("=", 70), collapse = ""))
    print(sprintf("[ Graph = %s ] ...", g_t))
    # Generate graph
    print("Start generating network...")
    if (g_t == "hub") {
      gene_per_clust <- 10
      num_clust <- 10
      num_hubs <- 5
      hub_degree <- 2
      num_other_edges <- 2
      net <- netSimulation(gene_per_clust, num_clust, num_hubs, hub_degree, num_other_edges)
    } else if (g_t == "power") {
      gene_per_clust <- 10
      num_clust <- 10
      low.strength <- 0.1
      sup.strength <- 0.9
      net <- powerLawSimulation(gene_per_clust, num_clust, low.strength, sup.strength)
    } else if (g_t == "GRN") {
      net <- GRNSimulation()
      if (!is.positive.definite(net)) {
        print("Net non-PD!")
        net <- nearPD(net, corr = TRUE, keepDiag = TRUE, base.matrix = TRUE, doSym = TRUE)$mat
        net <- as.matrix(forceSymmetric(net))
      }
    }
    tmp_net <- net
    tmp_net[abs(tmp_net) < 1e-1] <- 0.0
    net_g <- graph_from_adjacency_matrix(sign(abs(tmp_net)), mode = "undirected", diag = FALSE, weighted = TRUE)
    plot(net_g, vertex.label = NA, vertex.size = 3, edge.width = 3, main = "True Graph")
    write.csv(net, sprintf("%s-%s-net_mat.csv", file_name, g_t))

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
    write.csv(norta_simulated_data, sprintf("%s-%s-NORTA-data_mat.csv", file_name, g_t))
  }

}


# Simulate with pre-defined pars (different sparsity)
if (TRUE) {
  # Settings
  dist_name <- "nbinom" # zinb, nbinom, zip
  noise_loc <- NA
  prob <- 0.005
  size <- 0.5
  sparsity_list <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
  num_cells <- 500
  file_name <- sprintf("./data/diff_settings/diff_sparsity/100gene-hub")
  # -----
  # Generate data
  gene_per_clust <- 10
  num_clust <- 10
  num_hubs <- 5
  hub_degree <- 2
  num_other_edges <- 2
  gene_list <- paste0('gene', 1:(num_clust * gene_per_clust))
  cell_list <- paste0('cell', 1:(num_cells))

  net <- netSimulation(gene_per_clust, num_clust, num_hubs, hub_degree, num_other_edges)
  net_g <- graph_from_adjacency_matrix(sign(abs(net)), mode = "undirected", diag = FALSE, weighted = TRUE)
  plot(net_g, vertex.label = NA, vertex.size = 3, edge.width = 3, main = "True Graph")
  write.csv(net, sprintf("%s-net_mat.csv", file_name))
  # -----
  # Simulate data
  fitted_pars_list <- lapply(1:(num_clust * gene_per_clust), function(t) {
    p_size <- size + rnorm(1, mean = 0.001, sd = 0.001)
    p_prob <- prob + rnorm(1, mean = 0.001, sd = 0.001)
    if (p_size < 0) p_size <- size
    if (p_prob < 0) p_prob <- prob
    return(list(size = p_size, prob = p_prob))
  })
  print("Start generating data (NORTA)...")
  norta_simulated_data <- NORTASimulation(
    num_cells, net, rep(sprintf("q%s", dist_name), length(fitted_pars_list)), fitted_pars_list,
    noise = noise_loc, seed = 1
  )
  colnames(norta_simulated_data) <- gene_list
  rownames(norta_simulated_data) <- cell_list
  write.csv(norta_simulated_data, sprintf("%s-full-NORTA-data_mat.csv", file_name, s))
  # Evaluation
  sparsity <- mean(sign(abs(norta_simulated_data)))
  print(sprintf("Sparsity = %f [size=%f, prob=%f]", sparsity, size, prob))
  par(mfrow = c(1, 2))
  hist(norta_simulated_data, main = sprintf("Raw [size=%f, prob=%f]", size, prob))
  hist(log10(norta_simulated_data + 1), main = sprintf("Log-Transformed [size=%f, prob=%f]", size, prob))
  for (s in sparsity_list) {
    drop_data <- addDropouts(norta_simulated_data, s = s)
    sparsity <- mean(sign(abs(drop_data)))
    print(sprintf("Sparsity = %f [s=%f]", sparsity, s))
    par(mfrow = c(1, 2))
    hist(drop_data, main = sprintf("Raw [sparsity=%f]", sparsity))
    hist(log10(drop_data + 1), main = sprintf("Log-Transformed [sparsity=%f]", sparsity))
    # Save data
    colnames(drop_data) <- gene_list
    rownames(drop_data) <- cell_list
    print("Start saving data...")
    write.csv(drop_data, sprintf("%s-%.1fsparsity-NORTA-data_mat.csv", file_name, s))
  }
}


# Simulate with pars learned from mouse cortex data (more genes)
if (TRUE) {
  # Settings
  dist_name <- "nbinom" # zinb, nbinom, zip
  num_cells_ratio <- 1.0
  noise_loc <- NA
  # -----
  # Load mouse cortex data
  dir_path <- "./data/experimental/mouse_cortex/processed/expr/"
  exp_type <- "Cortex1" # Cortex1, Cortex2
  protocal_type <- "Smart_seq2" # 10xChromium, Smart_seq2
  num_genes <- "500hvg"
  print("Start loading data...")
  sc_data <- loadExprMat(dir_path, exp_type = exp_type, protocal_type = protocal_type, num_gene = num_genes)
  num_cells <- round(length(sc_data[, 1]) * num_cells_ratio)
  # -----
  # Generate graph
  print("Start generating network...")
  ene_per_clust <- 10
  num_clust <- 10
  num_hubs <- 5
  hub_degree <- 2
  num_other_edges <- 2
  num_blk <- 5
  net <- netSimulationBlks(gene_per_clust, num_clust, num_hubs, hub_degree, num_other_edges, num_blk)
  if (!is.positive.definite(net)) {
    print("Net non-PD!")
    net <- nearPD(net, corr = TRUE, keepDiag = TRUE, base.matrix = TRUE, doSym = TRUE)$mat
    net <- as.matrix(forceSymmetric(net))
  }
  net_g <- graph_from_adjacency_matrix(sign(abs(net)), mode = "undirected", diag = FALSE, weighted = TRUE)
  plot(net_g, vertex.label = NA, vertex.size = 3, edge.width = 3, main = "True Graph")
  # Simulate data
  print("Start generating data (NORTA)...")
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
  norta_simulated_data <- NORTASimulation(
    num_cells, net, rep(sprintf("q%s", dist_name), dim(sc_data)[2]), fitted_pars_list,
    noise = noise_loc, seed = 1
  )
  evaluateSim(sc_data, norta_simulated_data)
  visDataHist(sc_data, norta_simulated_data)
  # -----
  # Save data
  colnames(norta_simulated_data) <- colnames(sc_data)
  rownames(net) <- colnames(net) <- colnames(sc_data)
  rownames(norta_simulated_data) <- paste0('cell', 1:(num_cells))
  print("Start saving data...")
  file_name <- sprintf("./data/diff_settings/high_dim/NORTA/%s-%s-%s", exp_type, protocal_type, num_genes)
  write.csv(norta_simulated_data, sprintf("%s-NORTA-data_mat.csv", file_name))
  write.csv(net, sprintf("%s-net_mat.csv", file_name))
}


# Simulate with pars learned from PBMC data (more genes)
if (TRUE) {
  # Settings
  dist_name <- "nbinom" # zinb, nbinom, zip
  num_cells_ratio <- 1.0
  noise_loc <- NA
  # -----
  # Load mouse cortex data
  dir_path <- "./data/experimental/PBMC/processed/expr/"
  exp_type <- "pbmc1" # pbmc1, pbmc2
  protocal_type <- "inDrops" # Drop, inDrops
  num_genes <- "500hvg"
  print("Start loading data...")
  sc_data <- loadExprMat(dir_path, exp_type = exp_type, protocal_type = protocal_type, num_gene = num_genes)
  num_cells <- round(length(sc_data[, 1]) * num_cells_ratio)
  # -----
  # Generate graph
  gene_per_clust <- 10
  num_clust <- 10
  num_hubs <- 5
  hub_degree <- 2
  num_other_edges <- 2
  num_blk <- 5
  print("Start generating network...")
  net <- netSimulationBlks(gene_per_clust, num_clust, num_hubs, hub_degree, num_other_edges, num_blk)
  if (!is.positive.definite(net)) {
    print("Net non-PD!")
    net <- nearPD(net, corr = TRUE, keepDiag = TRUE, base.matrix = TRUE, doSym = TRUE)$mat
    net <- as.matrix(forceSymmetric(net))
  }
  net_g <- graph_from_adjacency_matrix(sign(abs(net)), mode = "undirected", diag = FALSE, weighted = TRUE)
  plot(net_g, vertex.label = NA, vertex.size = 3, edge.width = 3, main = "True Graph")
  # Simulate data
  print("Start generating data (NORTA)...")
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
  norta_simulated_data <- NORTASimulation(
    num_cells, net, rep(sprintf("q%s", dist_name), dim(sc_data)[2]), fitted_pars_list,
    noise = noise_loc, seed = 1
  )
  evaluateSim(sc_data, norta_simulated_data)
  visDataHist(sc_data, norta_simulated_data)
  # -----
  # Save data
  colnames(norta_simulated_data) <- colnames(sc_data)
  rownames(net) <- colnames(net) <- colnames(sc_data)
  rownames(norta_simulated_data) <- paste0('cell', 1:(num_cells))
  print("Start saving data...")
  file_name <- sprintf("./data/diff_settings/high_dim/NORTA/%s-%s-%s", exp_type, protocal_type, num_genes)
  write.csv(norta_simulated_data, sprintf("%s-NORTA-data_mat.csv", file_name))
  write.csv(net, sprintf("%s-net_mat.csv", file_name))
}