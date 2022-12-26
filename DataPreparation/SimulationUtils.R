# Description: Util function for data simulation.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(Matrix)
  library(MASS) # multivariate Gaussian
  library(matrixcalc) # for SPD check
  library(mc2d)

  # Codes for network generation, taken from:
  #   https://github.com/Vivianstats/scLink/blob/master/inst/docs/scripts/simulation.R
  # It uses the source file of pcorSimulator package (https://rdrr.io/cran/ldstatsHD/man/pcorSimulator.html)
  files.sources <- list.files("./DataPreparation/ZIGaussian_sim_utils/pcorSimulator", full.names = TRUE)
  sapply(files.sources, source)
})

# ======================================
#             Data Loading
# ======================================

loadExprMat <- function (dir_path, exp_type, protocal_type, num_gene){
  file_name <-sprintf("%s/%s-%s-%s.csv", dir_path, exp_type, protocal_type, num_gene)
  expr_mat <- read.csv(file_name, sep=",", header=1, row.names=1)
  expr_mat <- as.matrix(expr_mat)
  return (expr_mat)
}

# ======================================
#         Data Simulation
# ======================================

netSimulation <- function(gene_per_clust, num_clust, num_hubs, hub_degree, num_other_edges) {
  sim_cor <- pcorSimulatorSimple(
    nobs = 200,
    nclusters = num_clust,
    nnodesxcluster = rep(gene_per_clust, num_clust),
    pattern = "hub", nhubs = num_hubs, degree.hubs = hub_degree, nOtherEdges = num_other_edges
  )
  net <- sim_cor$omega
  net[net > 0] <- runif(length(net[net > 0]), min = 0.1, max = 0.5)
  net[net < 0] <- runif(length(net[net < 0]), min = -0.5, max = -0.1)
  diag(net) <- 1.0
  net <- as.matrix(forceSymmetric(net))
  return(net)
}


powerLawSimulation <- function(gene_per_clust, num_clust, low.strength, sup.strength) {
  files.sources <- list.files("./util/simulation/pcorSimulator", full.names = TRUE)
  sapply(files.sources, source)
  sim_cor <- pcorSimulatorSimple(
    nobs = 200,
    nclusters = num_clust,
    nnodesxcluster = rep(gene_per_clust, num_clust),
    pattern = "powerLaw", low.strength = low.strength, sup.strength = sup.strength
  )
  net <- sim_cor$omega
  net[net > 0] <- runif(length(net[net > 0]), min = 0.1, max = 0.5)
  net[net < 0] <- runif(length(net[net < 0]), min = -0.5, max = -0.1)
  diag(net) <- 1.0
  net <- as.matrix(forceSymmetric(net))
  return(net)
}


GRNSimulation <- function(){
  net <- read.csv("./data/simulation/SERGIO/100gene-true_network.csv", header=1, row.names = 1)
  net <- as.matrix(net)
  net[net != 0] <- runif(length(net[net > 0]), min = 0.1, max = 0.5)
  diag(net) <- 1.0
  net <- as.matrix(forceSymmetric(net))
  return(net)
}


netSimulationBlks <- function(gene_per_clust, num_clust, num_hubs, hub_degree, num_other_edges, num_blks) {
  blk_sim <- list()
  for (i in 1:num_blks){
    sim_cor <- pcorSimulatorSimple(
    nobs = 200,
    nclusters = num_clust,
    nnodesxcluster = rep(gene_per_clust, num_clust),
    pattern = "hub", nhubs = num_hubs, degree.hubs = hub_degree, nOtherEdges = num_other_edges
    )
    net <- sim_cor$omega
    net[net > 0] <- runif(length(net[net > 0]), min = 0.1, max = 0.5)
    net[net < 0] <- runif(length(net[net < 0]), min = -0.5, max = -0.1)
    diag(net) <- 1.0
    net <- as.matrix(forceSymmetric(net))
    blk_sim[[i]] <- net
  }
  comb_net <- as.matrix(bdiag(blk_sim))
  return (comb_net)
}

# ======================================
#         Simulations Evaluation
# ======================================

evaluateSimulationCov <- function(ref_cov, simulated_cov) {
  # Metrics of all elements
  ref_cov_label <- sign(ref_cov)
  simulated_cov_label <- sign(simulated_cov)
  full_mse_val <- mse(ref_cov, simulated_cov)
  full_accuracy_val <- accuracy(ref_cov_label, simulated_cov_label)
  # Metrics of existing edges
  edge_ind <- which(ref_cov != 0, arr.ind = T)
  ref_edge_cov <- ref_cov[edge_ind]
  simulated_edge_cov <- simulated_cov[edge_ind]
  ref_edge_label <- sign(ref_edge_cov)
  simulated_edge_label <- sign(simulated_edge_cov)
  true_edge_mse_val <- mse(ref_edge_cov, simulated_edge_cov)
  true_edge_accuracy_val <- accuracy(ref_edge_label, simulated_edge_label)
  ###
  return(list(
    full_mse = full_mse_val, full_accuracy = full_accuracy_val,
    edge_mse = true_edge_mse_val, edge_accuracy = true_edge_accuracy_val
  ))
}

dataStats <- function(data) {
  vec_data <- as.vector(data)
  sparsity <- mean(sign(abs(data)))
  max_val <- max(data)
  dispersion <- var(vec_data) / mean(vec_data)
  return(as.data.frame(list(sparsity = sparsity, max_val = max_val, dispersion = dispersion)))
}

evaluateSim <- function (ref_data, sim_data){
  # cov_stats <- evaluateSimulationCov(net, cov(sim_data))
  # print("Recovery of covariance : ")
  # print(cov_stats)
  ref_stats <- dataStats(ref_data)
  simulated_stats <- dataStats(sim_data)
  print("Reference data stats : ")
  print(ref_stats)
  print("Simulated data stats : ")
  print(simulated_stats)
}

# ======================================
#            Visualization
# ======================================

visDataHist <- function (ref_data, sim_data){
  # Visualization
  par(mfrow = c(1, 2))
  sc_data_filtered <- ref_data[ref_data > 0]
  simulated_data_filtered <- sim_data[sim_data > 0]
  hist(sc_data_filtered)
  hist(simulated_data_filtered)
  par(mfrow = c(1, 2))
  hist(log10(sc_data_filtered + 1))
  hist(log10(simulated_data_filtered + 1))
}

# ======================================
#    Util Functions for ZI-Poisson
# ======================================

vec2tri <- function(k, p) {
  i = ceiling(0.5 * (2 * p - 1 - sqrt((2 * p - 1)^2 - 8 * k)))
  j = k - p * (i - 1) + i * (i - 1) / 2 + i
  return(as.matrix(cbind(i, j)))
}


addBernDropouts <- function(cnt, pi) {
  cnt_shape <- dim(cnt)
  num_cells <- cnt_shape[1]
  num_genes <- cnt_shape[2]
  dropind <- rbern(num_cells * num_genes, pi)
  cnt_zero <- cnt * dropind
  return(cnt_zero)
}


zilgm_sim <- function(A, n, p, zlvs, family = c("poisson", "negbin"), signal, theta = NULL, noise) {
  family = match.arg(family)
  is.symm = TRUE

  npair = p * (p - 1) / 2
  Y = E = matrix(0, nrow = n, ncol = p)
  Ypair = matrix(0, n, npair)

  if (family == "poisson") {
    for (j in 1:p)
    {
      Y[, j] = rpois(n = n, lambda = signal) # rZeroPoisson(p+pair, w=pi0, mu=mu_t)
      E[, j] = rpois(n = n, lambda = noise) #rZeroPoisson(m, w=pi0, mu=mu_n) # rpois(p, mu_n)
    }

    for (j in 1:npair)
    {
      Ypair[, j] = rpois(n = n, lambda = signal) # rZeroPoisson(p + pair, w = pi0, mu = mu_t)
    }
  }

  if (family == "negbin") {
    if (is.null(theta)) { stop("Need a theta") }

    for (j in 1:p)
    {
      Y[, j] = rnegbin(n = n, mu = signal, theta = theta) # rZeroPoisson(p+pair, w=pi0, mu=mu_t)
      E[, j] = rpois(n = n, lambda = noise) #rZeroPoisson(m, w=pi0, mu=mu_n) # rpois(p, mu_n)
    }

    for (j in 1:npair)
    {
      Ypair[, j] = rnegbin(n = n, mu = signal, theta = theta) # rZeroPoisson(p + pair, w = pi0, mu = mu_t)
    }
  }

  Y = cbind(Y, Ypair)

  # B is real signal
  tri_A = t(A)[lower.tri(t(A), diag = FALSE)]
  B = matrix(0, nrow = npair, ncol = p)
  ix = vec2tri(1:npair, p)
  adj = cbind(ix, tri_A)

  nadj = (1:npair)[adj[, 3] != 0]
  for (i in nadj)
  {
    B[i, adj[i, 2]] = 1

    if (is.symm) { B[i, adj[i, 1]] = 1 }
  }
  B = rbind(diag(1, p), B)
  X = Y %*% B + E

  for (j in 1:p)
  {
    X[runif(n) <= zlvs, j] = 0
  }
  return(list(A = A, X = X))
}

# ======================================
#    Util Functions for Imputation
# ======================================

loadData <- function(filename) {
  # filename <- sprintf("%s/%s-%s-data_mat.csv", data_path, data_name, sim_name)
  expr_mat <- read.csv(filename, header = 1, row.names = 1)
  return(expr_mat)
}