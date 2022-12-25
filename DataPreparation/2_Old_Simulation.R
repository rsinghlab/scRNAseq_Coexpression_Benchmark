# Description: ZI-Gaussian and ZI-Poisson simulation processes.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(igraph)
  library(Matrix) # symmetric matrix
  library(Metrics)
  library(MASS) # multivariate Gaussian
  library(matrixcalc) # for SPD check
  library(mc2d)
  source("./DataPreparation/SimulationUtils.R") #TODO: util functions for scLink
})

# ======================================
# The simulation codes are taken from:
#   https://github.com/Vivianstats/scLink/blob/master/inst/docs/scripts/simulation.R

ZIGaussianSim <- function(net, num_cells, rho) {
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
  dem[which(is.na(dem))] <- 0.0
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

# ======================================
# The simulation codes are taken from:
#   https://github.com/bbeomjin/ZILGM/blob/304a638c63755b28ffbca2709b8c95dbf0b1b681/R/network_gen.R

ZIPoissonSim <- function(net, num_cells, pi) {
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
#TODO: file path

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
  # ZI-Poisson simulation
  for (pi in c(1.0, 0.9, 0.8)) {
    zipoisson_simulated_data <- ZIPoissonSim(net, num_cells = num_cells, pi = pi)
    print(sprintf("[pi = %f] Sparsity = %f", pi, mean(sign(zipoisson_simulated_data))))
    par(mfrow = c(1, 2))
    hist(zipoisson_simulated_data)
    hist(log10(zipoisson_simulated_data + 1))
    # Save data
    colnames(zipoisson_simulated_data) <- paste0('gene', 1:(num_genes))
    rownames(zipoisson_simulated_data) <- paste0('cell', 1:(num_cells))
    write.csv(zipoisson_simulated_data, sprintf("%s-ZILGM_%.1fpi-data_mat.csv", file_name, pi))
  }
  # -----
  # Zi-Gaussian simulation
  for (rho in c(0.07, 0.1, 0.13, 0.16)) {
    zigaussian_simulated_data <- ZIGaussianSim(net, num_cells = num_cells, rho = rho)
    print(sprintf("[rho = %f] Sparsity = %f", rho, mean(sign(zigaussian_simulated_data))))
    par(mfrow = c(1, 2))
    hist(zigaussian_simulated_data)
    hist(log10(zigaussian_simulated_data + 1))
    # Save data
    colnames(zigaussian_simulated_data) <- paste0('gene', 1:(num_genes))
    rownames(zigaussian_simulated_data) <- paste0('cell', 1:(num_cells))
    write.csv(zigaussian_simulated_data, sprintf("%s-scLink_%.2frho-data_mat.csv", file_name, rho))
  }
}