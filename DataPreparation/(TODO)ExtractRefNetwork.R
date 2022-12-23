library(Matrix)
library(MASS)
library(readxl)

# ===============================

loadGeneList <- function(dir_path) {
  data_file_list <- list.files(dir_path)
  temp_data_file_list <- c()
  for (each in data_file_list) {
    if (grepl("full_genes", each, fixed = TRUE) & grepl("Smart_seq2", each, fixed = TRUE)) { # 10xChromium, Smart_seq2
      temp_data_file_list <- c(temp_data_file_list, each)
    }
  }
  data_file_list <- temp_data_file_list
  gene_list <- lapply(data_file_list, function(filename) {
    con <- file(sprintf("%s/%s", dir_path, filename), "r")
    header_line <- readLines(con, n = 1)
    close(con)
    header_line <- strsplit(header_line, ",")[[1]]
    header_line <- header_line[2:length(header_line)]
    tmp_gene_list <- unlist(lapply(header_line, function(s) { return(strsplit(s, "_")[[1]][2]) }))
    return(list(strsplit(filename, '[.]')[[1]][1], tmp_gene_list))
  })
  return(gene_list)
}


loadTabulaGeneList <- function(dir_path, cell_type) {
  con <- file(sprintf("%s/%s.csv", dir_path, cell_type), "r")
  header_line <- readLines(con, n = 1)
  close(con)
  header_line <- strsplit(header_line, ",")[[1]]
  gene_list <- header_line[2:length(header_line)]
  return(gene_list)
}

# ===============================

loadTRRUSTv2 <- function(gene_list, species) {
  if (species == "mouse") {
    net_table <- read.table(file = './data/experimental/mouse_cortex/network_reference/TRRUSTv2/trrust_rawdata.mouse.tsv', sep = '\t', header = FALSE)
  } else if (species == "human") {
    net_table <- read.table(file = './data/experimental/PBMC/network_reference/TRRUSTv2/trrust_rawdata.human.tsv', sep = '\t', header = FALSE)
  }
  # net_table <- read.table(file = './data/experimental/mouse_cortex/network_reference/TRRUSTv2/trrust_rawdata.mouse.tsv', sep = '\t', header = FALSE)
  net_table <- net_table[, c(1, 2)]
  net_table <- net_table[which(net_table$V1 %in% gene_list),]
  net_table <- net_table[which(net_table$V2 %in% gene_list),]
  return(net_table)
}


loadMouseNetv2 <- function(gene_list) {
  net_table <- read.table(file = './data/experimental/mouse_cortex/network_reference/MouseNetv2/MouseNetV2_symbol_rat.txt', sep = '\t', header = FALSE)
  net_table <- net_table[, c(1, 2)]
  net_table <- net_table[which(net_table$V1 %in% gene_list),]
  net_table <- net_table[which(net_table$V2 %in% gene_list),]
  return(net_table)
}


loadSTRING <- function(gene_list, species) {
  if (species == "mouse") {
    link_table <- read.table("./data/experimental/mouse_cortex/network_reference/STRING/protein_links.txt", sep = " ", header = 1)
    protein_map <- read.csv("./data/experimental/mouse_cortex/network_reference/STRING/protein_info.txt", sep = "\t", colClasses = c(NA, NA, NA, "NULL"))
  } else if (species == "human") {
    link_table <- read.table("./data/experimental/PBMC/network_reference/STRING/protein_links.txt", sep = " ", header = 1)
    protein_map <- read.csv("./data/experimental/PBMC/network_reference/STRING/protein_info.txt", sep = "\t", colClasses = c(NA, NA, NA, "NULL"))
  } else {
    stop(sprintf("Unknown species %s!", species))
  }
  intersected_gene <- intersect(gene_list, protein_map[, "preferred_name"])
  selected_protein_map <- protein_map[which(protein_map$preferred_name %in% intersected_gene),]
  link_table <- link_table[which(link_table$protein1 == selected_protein_map$X.string_protein_id),]
  name_link_table <- link_table
  print("Start mapping...")
  for (i in 1:dim(name_link_table)[1]) {
    each <- name_link_table[i,]
    node1_id <- each$protein1
    node2_id <- each$protein2
    node1_name <- protein_map[which(protein_map$X.string_protein_id == node1_id),]
    if (dim(node1_name)[1] == 0) {
      node1_name <- NA
    } else {
      node1_name <- node1_name$preferred_name
    }
    node2_name <- protein_map[which(protein_map$X.string_protein_id == node2_id),]
    if (dim(node2_name)[1] == 0) {
      node2_name <- NA
    } else {
      node2_name <- node2_name$preferred_name
    }
    name_link_table[i,]$protein1 <- node1_name
    name_link_table[i,]$protein2 <- node2_name
  }
  name_link_table <- name_link_table[which(!is.na(name_link_table$protein2)),]
  name_link_table <- name_link_table[which(name_link_table$protein2 %in% intersected_gene),]
  return(name_link_table)
}

# ===============================

convert2Mat <- function(gene_list, edge_table) {
  adj_mat <- diag(length(gene_list))
  colnames(adj_mat) <- gene_list
  rownames(adj_mat) <- gene_list
  if (length(edge_table[, 1]) == 0) {
    return(adj_mat)
  }
  for (i in 1:dim(edge_table)[1]) {
    node1 <- edge_table[i, 1]
    node2 <- edge_table[i, 2]
    adj_mat[node1, node2] <- 1
    adj_mat[node2, node1] <- 1
  }
  return(adj_mat)
}

# ===============================

if (FALSE) {
  species <- "mouse" # human, mouse
  if (species == "mouse") {
    dir_path <- "./data/experimental/mouse_cortex/processed/"
  } else if (species == "human") {
    dir_path <- "./data/experimental/PBMC/processed/"
  } else {
    stop(sprintf("Unknown species %s!", species))
  }
  gene_list <- loadGeneList(sprintf("%s/expr/", dir_path))
  print(sprintf("[ %s ] %s", species, dir_path))
  # -----
  for (each_file in gene_list) {
    print("--------------------------------")
    filename <- each_file[[1]]
    data_gene_list <- each_file[[2]]
    print(each_file[[1]])
    if (species == "mouse") {
      mouse_net <- loadMouseNetv2(data_gene_list)
      print(sprintf("The num of edges in MouseNetv2 : %d", dim(mouse_net)[1]))
      # Convert to network adjacency matrix
      mouse_mat <- convert2Mat(data_gene_list, mouse_net)
      print(sprintf("Num of elements in MouseNetv2 mat : %d", sum(mouse_mat) - length(data_gene_list)))
      writeMM(as(as.matrix(mouse_mat), "sparseMatrix"), sprintf("%s/network/%s-MouseNet_net.mtx", dir_path, filename))
    }
    trust_net <- loadTRRUSTv2(data_gene_list, species)
    print(sprintf("The num of edges in TRUSTv2 : %d", dim(trust_net)[1]))
    trust_mat <- convert2Mat(data_gene_list, trust_net)
    print(sprintf("Num of elements in TRUSTv2 mat : %d", sum(trust_mat) - length(data_gene_list)))
    writeMM(as(as.matrix(trust_mat), "sparseMatrix"), sprintf("%s/network/%s-TRRUST_net.mtx", dir_path, filename))

    string_net <- loadSTRING(data_gene_list, species)
    print(sprintf("The num of edges in STRING : %d", dim(string_net)[1]))
    string_mat <- convert2Mat(data_gene_list, string_net)
    print(sprintf("Num of elements in STRING mat : %d", sum(string_mat) - length(data_gene_list)))
    writeMM(as(as.matrix(string_mat), "sparseMatrix"), sprintf("%s/network/%s-STRING_net.mtx", dir_path, filename))
  }
}


# For Tabula Muris data
if (FALSE) {
  dir_path <- "./data/Tabula_Muris/500hvg/"
  cell_type_list <- c("T_cell", "skeletal_muscle_satellite_stem_cell", "type_B_pancreatic_cell")
  for (t in cell_type_list) {
    print(sprintf("[ %s ] %s", dir_path, t))
    gene_list <- loadTabulaGeneList(dir_path, t)
    # -----
    string_net <- loadSTRING(gene_list, species="mouse")
    print(sprintf("The num of edges in STRING : %d", dim(string_net)[1]))
    string_mat <- convert2Mat(gene_list, string_net)
    print(sprintf("Num of elements in STRING mat : %d", sum(string_mat) - length(gene_list)))
    writeMM(as(as.matrix(string_mat), "sparseMatrix"), sprintf("%s/%s-STRING_net.mtx", dir_path, t))
    # -----
    trust_net <- loadTRRUSTv2(gene_list, species="mouse")
    print(sprintf("The num of edges in TRUSTv2 : %d", dim(trust_net)[1]))
    trust_mat <- convert2Mat(gene_list, trust_net)
    print(sprintf("Num of elements in TRUSTv2 mat : %d", sum(trust_mat) - length(gene_list)))
    writeMM(as(as.matrix(trust_mat), "sparseMatrix"), sprintf("%s/%s-TRRUST_net.mtx", dir_path, t))
  }
}