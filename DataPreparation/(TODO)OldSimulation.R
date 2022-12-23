suppressPackageStartupMessages({
  library(igraph)
  library(Matrix) # symmetric matrix
  library(Metrics)
  library(MASS) # multivariate Gaussian
  library(matrixcalc) # for SPD check
  library(mc2d)
  source("./util/SimulationUtils.R")
})

# ======================================

scLinkSimulate <- function(net, num_cells, rho) {
  theta = net
  num_genes = length(net[1,])
  # simulate data
  pars = readRDS("./util/simulation/mixpara.rds")
  pars = pars[complete.cases(pars),]
  pars = pars[, c("rate", "mu", "sigma")]
  pars = pars[order(pars[, "mu"], decreasing = TRUE)[1:1000],]
  # pars = pars[order(pars[,"mu"],decreasing = TRUE)[sample(dim(pars)[1], 1000, replace=TRUE)], ]
  pars_selected = pars[sample(1:1000, num_genes),]
  sigma = pars_selected[, "sigma"]

  cov = solve(theta)
  dem = sqrt(diag(cov))
  dem[which(is.na(dem))] <- 0.0 #TODO: note this!
  cov = sweep(cov, MARGIN = 2, dem, FUN = "/")
  cov = sweep(cov, MARGIN = 1, dem, FUN = "/")
  cov = sweep(cov, MARGIN = 2, sigma, FUN = "*")
  cov = sweep(cov, MARGIN = 1, sigma, FUN = "*")
  cov[which(is.na(cov))] <- 0.0 #TODO: note this!
  cov[which(is.infinite(cov))] <- 0.0 #TODO: note this!
  cov <- as.matrix(forceSymmetric(cov))
  # Convert covariance to PD matrix
  if (!is.positive.definite(cov)) {
    print("covariance matrix is not PD! Added identity matrix")
    cov = cov + eigen(cov)$values[1] * diag(dim(cov)[1])
  }
  theta = solve(cov)
  theta[abs(theta) < 1e-5] = 0
  theta <- as.matrix(forceSymmetric(theta))
  if (!is.positive.definite(theta)) stop("precision matrix is not PD!")
  ### simulate count matrices for all time points -----------------------
  mu = pars_selected[, "mu"]
  cnt = mvrnorm(n = num_cells, mu = mu, Sigma = cov)
  cnt[cnt < 0] = 0
  cnt_zero = cnt
  droprate = exp(-rho * cnt^2)
  dropind = rbinom(prod(dim(droprate)), size = 1, prob = as.vector(droprate))
  dropvals = sample(c(0, log10(2)), sum(dropind == 1), prob = c(0.7, 0.3), replace = TRUE)
  cnt_zero[dropind == 1] = dropvals
  return(cnt_zero)
}

gaussianSimulate <- function(net, num_cells) {
  num_genes <- length(num[1,])
  simulated_data <- mvrnorm(num_genes)
}

# ======================================
vec2tri = function(k, p) {
  i = ceiling(0.5 * (2 * p - 1 - sqrt((2 * p - 1)^2 - 8 * k)))
  j = k - p * (i - 1) + i * (i - 1) / 2 + i
  return(as.matrix(cbind(i, j)))
}


addDropouts <- function(cnt, rho) {
  cnt[cnt < 0] <- 0
  cnt_zero <- cnt
  droprate <- exp(-rho * cnt^2)
  dropind <- rbinom(prod(dim(droprate)), size = 1, prob = as.vector(droprate))
  cnt_zero[dropind == 1] <- 0.0
  return(cnt_zero)
}

addBernDropouts <- function(cnt, pi) {
  cnt_shape <- dim(cnt)
  num_cells <- cnt_shape[1]
  num_genes <- cnt_shape[2]
  dropind <- rbern(num_cells * num_genes, pi)
  cnt_zero <- cnt * dropind
  return(cnt_zero)
}

zilgm_sim = function(A, n, p, zlvs, family = c("poisson", "negbin"),
                     signal, theta = NULL, noise)
{
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


ZILGMSimulate <- function(net, num_cells, pi) {
  num_genes <- length(net[1,])
  zlvs <- 0.0
  sim_res <- zilgm_sim(net, num_cells, num_genes, zlvs, family = "poisson", signal = 1.5, noise = 0.5)
  # Add dropouts
  cnt_data <- sim_res$X
  if (!is.na(pi)) {
    sim_res_w_dropout <- addBernDropouts(cnt_data, pi)
    return(sim_res_w_dropout)
  } else {
    return(cnt_data)
  }
}

# ======================================

#DEPRECATED
if (FALSE) {
  # Settings
  dist_name <- "nbinom" # zinb, nbinom, zip
  num_cells <- 1000
  num_genes <- "100hvg" # 50hvg, 100hvg
  if (num_genes == "50hvg") {
    gene_per_clust <- 10
    num_clust <- 5
    num_hubs <- 5
    hub_degree <- 2
    num_other_edges <- 2
  } else if (num_genes == "100hvg") {
    gene_per_clust <- 10
    num_clust <- 10
    num_hubs <- 5
    hub_degree <- 2
    num_other_edges <- 2
  }
  rho <- 0.07
  noise_loc <- NA
  # -----
  net <- netSimulation(gene_per_clust, num_clust, num_hubs, hub_degree, num_other_edges)
  net_g <- graph_from_adjacency_matrix(sign(abs(net)), mode = "undirected", diag = FALSE, weighted = TRUE)
  plot(net_g, vertex.label = NA, vertex.size = 3, edge.width = 3, main = "True Graph")
  # -----
  # gaussian_simulated_data <- gaussianSimulate(net, num_cells = 500)
  scLink_simulated_data <- scLinkSimulate(net, num_cells = num_cells, rho = rho)
  zilgm_simulated_data <- ZILGMSimulate(net, num_cells = num_cells, rho = rho)
  par(mfrow = c(1, 2))
  hist(scLink_simulated_data)
  hist(log10(scLink_simulated_data + 1))
  par(mfrow = c(1, 2))
  hist(zilgm_simulated_data)
  hist(log10(zilgm_simulated_data + 1))
  # Save data
  colnames(scLink_simulated_data) <- colnames(zilgm_simulated_data) <- paste0('gene', 1:(gene_per_clust*num_clust))
  rownames(net) <- colnames(net) <- paste0('gene', 1:(gene_per_clust*num_clust))
  rownames(scLink_simulated_data) <- rownames(zilgm_simulated_data) <- paste0('cell', 1:(num_cells))
  print("Start saving data...")
  file_name <- sprintf("./data/simulated/old/%s", num_genes)
  # write.csv(norta_simulated_data, sprintf("%s-Gaussian-data_mat.csv", gaussian_simulated_data))
  write.csv(scLink_simulated_data, sprintf("%s-scLink_rho%f-data_mat.csv", file_name, rho))
  write.csv(zilgm_simulated_data, sprintf("%s-ZILGM_rho%f-data_mat.csv", file_name, rho))
  write.csv(net, sprintf("%s-net_mat.csv", file_name))
}


if (FALSE) {
  num_genes <- 100
  num_cells <- 2 * num_genes
  file_name <- sprintf("./data/simulated/old/%dhvg", num_genes)
  # -----
  net <- read.csv(sprintf("./data/simulated/old/%dhvg-net_mat.csv", num_genes), header = 1, row.names = 1)
  net <- as.matrix(net)
  net_g <- graph_from_adjacency_matrix(sign(abs(net)), mode = "undirected", diag = FALSE, weighted = TRUE)
  plot(net_g, vertex.label = NA, vertex.size = 3, edge.width = 3, main = "True Graph")
  rownames(net) <- colnames(net) <- paste0('gene', 1:(num_genes))
  # -----
  for (pi in c(1.0, 0.9, 0.8)) {
    zilgm_simulated_data <- ZILGMSimulate(net, num_cells = num_cells, pi = pi)
    print(sprintf("[pi = %f] Sparsity = %f", pi, mean(sign(zilgm_simulated_data))))
    par(mfrow = c(1, 2))
    hist(zilgm_simulated_data)
    hist(log10(zilgm_simulated_data + 1))
    # Save data
    colnames(zilgm_simulated_data) <- paste0('gene', 1:(num_genes))
    rownames(zilgm_simulated_data) <- paste0('cell', 1:(num_cells))
    write.csv(zilgm_simulated_data, sprintf("%s-ZILGM_%.1fpi-data_mat.csv", file_name, pi))
  }
  # -----
  for (rho in c(0.07, 0.1, 0.13, 0.16)) {
    scLink_simulated_data <- scLinkSimulate(net, num_cells = num_cells, rho = rho)
    print(sprintf("[rho = %f] Sparsity = %f", rho, mean(sign(scLink_simulated_data))))
    par(mfrow = c(1, 2))
    hist(scLink_simulated_data)
    hist(log10(scLink_simulated_data + 1))
    # Save data
    colnames(scLink_simulated_data) <- paste0('gene', 1:(num_genes))
    rownames(scLink_simulated_data) <- paste0('cell', 1:(num_cells))
    write.csv(scLink_simulated_data, sprintf("%s-scLink_%.2frho-data_mat.csv", file_name, rho))
  }
}