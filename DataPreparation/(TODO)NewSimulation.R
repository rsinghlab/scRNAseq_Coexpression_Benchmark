suppressPackageStartupMessages({
  library(SimCorMultRes) # NORTA algorithm
  library(fitdistrplus) # distribution fitting
  # library(gamlss) # for ZIP
  library(ZIM) # for ZIP
  library(VGAM) # for ZINB
  library(igraph) # network visualization
  library(Matrix) # symmetric matrix
  library(Metrics)
  library(netSmooth) # netSmooth imputation
  source("./util/SimulationUtils.R")
})

# ======================================
#    Marginal Distribution Parameters
# ======================================

fitMarginal <- function(data, dist_name, visualize = FALSE) {
  if (dist_name == "nbinom") {
    fit_res <- fitdist(data, 'nbinom', start = list(mu = 1, size = 0.1), method = "mme")
  } else if (dist_name == "zinb") {
    fit_res <- zim(data ~ 1, control = zim.control(dist = "zinb", type = "solve")) # solve, ginv
  } else if (dist_name == "zip") {
    # stop("Under development...")
    fit_res <- zim(data ~ 1, control = zim.control(dist = "zip", type = "solve")) # solve, ginv
  } else if (dist_name == "exp") {
    fit_res <- fitdist(data, 'exp', start = list(rate = 1), method = "mme")
  }
  if (visualize) {
    if (dist_name == "nbinom") {
      # par(mar = c(2, 2, 2, 2))
      par(mar = c(1, 1, 1, 1))
      plot(fit_res)
    } else if (dist_name == "exp") {
      par(mar = c(0, 0, 0, 0))
      plot(fit_res)
    } else if (dist_name == "zinb") {
      ecdf_func <- ecdf(data)
      plot(sort(data), lapply(sort(data), ecdf_func))
      x_list <- seq(from = min(data), to = max(data), length.out = 1000)
      y_list <- lapply(x_list, function(x) { return(pzinb(x, omega = fit_res$omega[1], k = fit_res$k, lambda = fit_res$lambda[1])) })
      lines(x_list, y_list, col = "red")
    } else if (dist_name == "zip") {
      ecdf_func <- ecdf(data)
      plot(sort(data), lapply(sort(data), ecdf_func))
      x_list <- seq(from = min(data), to = max(data), length.out = 1000)
      y_list <- lapply(x_list, function(x) { return(pzip(x, omega = fit_res$omega[1], lambda = fit_res$lambda[1])) })
      lines(x_list, y_list, col = "red")
    }
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


netSmoothSimulation <- function(ref_data, net) {
  net <- as.matrix(net)
  net[which(abs(net) < 1e-4)] <- 0.0
  net[which(is.na(net))] <- 0.0
  diag(net) <- 1.0
  net <- sign(abs(net))
  rownames(net) <- colnames(net) <- colnames(ref_data)
  # Imputation
  print("Start imputation...")
  imputed_data <- netSmooth(t(ref_data), net, alpha = 0.9)
  imputed_data[imputed_data < 0.1] <- 0
  imputed_data <- round(imputed_data)
  return(t(imputed_data))
}


# ======================================

# Simulate with pars learned from mouse cortex data
if (FALSE) {
  # Settings
  dist_name <- "nbinom" # zinb, nbinom, zip
  num_cells_ratio <- 1.0
  noise_loc <- NA
  # -----
  # Load mouse cortex data
  dir_path <- "./data/experimental/mouse_cortex/processed/expr/"
  exp_type <- "Cortex1" # Cortex1, Cortex2
  protocal_type <- "10xChromium" # 10xChromium, Smart_seq2
  num_genes <- "100hvg" # 50hvg, 100hvg, 500hvg, 1000hvg
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
  } else if (num_genes == "500hvg") {
    gene_per_clust <- 20
    num_clust <- 25
    num_hubs <- 10
    hub_degree <- 4
    num_other_edges <- 4
  } else if (num_genes == "1000hvg") {
    gene_per_clust <- 20
    num_clust <- 50
    num_hubs <- 5
    hub_degree <- 4
    num_other_edges <- 2
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
  # Simulate data
  #TODO: add noise
  print("Start generating data (NORTA)...")
  # Fit marginal distributions
  print("Start fitting marginal distributions...")
  fitted_pars_list <- list()
  fitted_dist_list <- list()
  for (i in 1:dim(sc_data)[2]) {
    # for (i in 1:5) {
    gene_data <- as.vector(sc_data[, i])
    gene_dist <- fitMarginal(gene_data, dist_name, visualize = FALSE)
    gene_pars <- as.list(gene_dist$estimate)
    if (dist_name == "nbinom") {
      fitted_pars_list[[i]] <- list(size = gene_pars$size, prob = gene_pars$size / (gene_pars$size + gene_pars$mu))
    } else if (dist_name == "zinb") {
      fitted_pars_list[[i]] <- list(omega = gene_pars$omega[1], k = gene_pars$k, lambda = gene_pars$lambda[1])
    } else if (dist_name == "zip") {
      fitted_pars_list[[i]] <- list(omega = gene_pars$omega[1], lambda = gene_pars$lambda[1])
    } else {
      stop(sprintf("Unimplemented parameter parsing for %s!", dist_name))
    }
    fitted_dist_list[[i]] <- gene_dist
  }
  norta_simulated_data <- NORTASimulation(
    num_cells, net, rep(sprintf("q%s", dist_name), dim(sc_data)[2]), fitted_pars_list,
    noise = noise_loc, seed = 1
  )
  evaluateSim(sc_data, norta_simulated_data)
  visDataHist(sc_data, norta_simulated_data)
  # -----
  # Evaluation
  netSmooth_simulated_data <- netSmoothSimulation(sc_data, net)
  evaluateSim(sc_data, netSmooth_simulated_data)
  visDataHist(sc_data, netSmooth_simulated_data)
  # Save data
  colnames(norta_simulated_data) <- colnames(netSmooth_simulated_data) <- colnames(sc_data)
  rownames(net) <- colnames(net) <- colnames(sc_data)
  rownames(norta_simulated_data) <- rownames(netSmooth_simulated_data) <- paste0('cell', 1:(num_cells))
  print("Start saving data...")
  file_name <- sprintf("./data/simulated/new/%s-%s-%s", exp_type, protocal_type, num_genes)
  write.csv(norta_simulated_data, sprintf("%s-NORTA-data_mat.csv", file_name))
  write.csv(netSmooth_simulated_data, sprintf("%s-netSmooth-data_mat.csv", file_name))
  write.csv(net, sprintf("%s-net_mat.csv", file_name))
}


# Simulate with pars learned from PBMC data
if (FALSE) {
  # Settings
  dist_name <- "nbinom" # zinb, nbinom, zip
  num_cells_ratio <- 1.0
  noise_loc <- NA
  # -----
  # Load mouse cortex data
  dir_path <- "./data/experimental/PBMC/processed/expr/"
  exp_type <- "pbmc1" # pbmc1, pbmc2
  protocal_type <- "inDrops" # Drop, inDrops
  num_genes <- "50hvg" # 50hvg, 100hvg, 500hvg, 1000hvg
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
  } else if (num_genes == "500hvg") {
    gene_per_clust <- 20
    num_clust <- 25
    num_hubs <- 10
    hub_degree <- 4
    num_other_edges <- 4
  } else if (num_genes == "1000hvg") {
    gene_per_clust <- 20
    num_clust <- 50
    num_hubs <- 5
    hub_degree <- 4
    num_other_edges <- 2
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
  # Simulate data
  #TODO: add noise
  print("Start generating data (NORTA)...")
  # Fit marginal distributions
  print("Start fitting marginal distributions...")
  fitted_pars_list <- list()
  fitted_dist_list <- list()
  for (i in 1:dim(sc_data)[2]) {
    # for (i in 1:5) {
    gene_data <- as.vector(sc_data[, i])
    gene_dist <- fitMarginal(gene_data, dist_name, visualize = FALSE)
    gene_pars <- as.list(gene_dist$estimate)
    if (dist_name == "nbinom") {
      fitted_pars_list[[i]] <- list(size = gene_pars$size, prob = gene_pars$size / (gene_pars$size + gene_pars$mu))
    } else if (dist_name == "zinb") {
      fitted_pars_list[[i]] <- list(omega = gene_pars$omega[1], k = gene_pars$k, lambda = gene_pars$lambda[1])
    } else if (dist_name == "zip") {
      fitted_pars_list[[i]] <- list(omega = gene_pars$omega[1], lambda = gene_pars$lambda[1])
    } else {
      stop(sprintf("Unimplemented parameter parsing for %s!", dist_name))
    }
    fitted_dist_list[[i]] <- gene_dist
  }
  norta_simulated_data <- NORTASimulation(
    num_cells, net, rep(sprintf("q%s", dist_name), dim(sc_data)[2]), fitted_pars_list,
    noise = noise_loc, seed = 1
  )
  evaluateSim(sc_data, norta_simulated_data)
  visDataHist(sc_data, norta_simulated_data)
  # -----
  # Evaluation
  netSmooth_simulated_data <- netSmoothSimulation(sc_data, net)
  evaluateSim(sc_data, netSmooth_simulated_data)
  visDataHist(sc_data, netSmooth_simulated_data)
  # Save data
  colnames(norta_simulated_data) <- colnames(netSmooth_simulated_data) <- colnames(sc_data)
  rownames(net) <- colnames(net) <- colnames(sc_data)
  rownames(norta_simulated_data) <- rownames(netSmooth_simulated_data) <- paste0('cell', 1:(num_cells))
  print("Start saving data...")
  file_name <- sprintf("./data/simulated/new/%s-%s-%s", exp_type, protocal_type, num_genes)
  write.csv(norta_simulated_data, sprintf("%s-NORTA-data_mat.csv", file_name))
  write.csv(netSmooth_simulated_data, sprintf("%s-netSmooth-data_mat.csv", file_name))
  write.csv(net, sprintf("%s-net_mat.csv", file_name))
}


