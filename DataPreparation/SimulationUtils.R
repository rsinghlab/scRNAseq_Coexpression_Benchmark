# Description: Util function for data simulation.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>

suppressPackageStartupMessages({
  library(Matrix)
  library(MASS) # multivariate Gaussian
  library(matrixcalc) # for SPD check
  library(mc2d)
})

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
  filename <- sprintf("%s/%s-%s-data_mat.csv", data_path, data_name, sim_name)
  expr_mat <- read.csv(filename, header = 1, row.names = 1)
  return(expr_mat)
}