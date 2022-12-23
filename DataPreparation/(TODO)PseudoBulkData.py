import torch
import itertools
import pandas as pd
import numpy as np
import scanpy



def loadData(data_path, data_name, sim_name):
    filename = "{}/{}-{}-data_mat.csv".format(data_path, data_name, sim_name)
    expr_mat = pd.read_csv(filename, header=0, index_col=0)
    return expr_mat


def loadSERGIOData(data_path, data_name):
    filename = "{}/{}-{}sparsity.csv".format(data_path, data_name[0], data_name[1])
    expr_mat = pd.read_csv(filename, header=0, index_col=0)
    return expr_mat


def pseudoBulk(X):
    X = scanpy.AnnData(X)
    scanpy.pp.neighbors(X, n_neighbors=10, n_pcs=50)
    scanpy.tl.umap(X)
    scanpy.tl.leiden(X, resolution=0.5, key_added="type")
    # scanpy.pl.umap(X, color=["type"])
    # -----
    cluster_label = X.obs['type'].values
    unique_cluster_label = np.unique(cluster_label)
    res = pd.DataFrame(columns=X.var_names, index=unique_cluster_label)
    for clust in unique_cluster_label:
        # res.loc[clust] = X[np.where(cluster_label == clust)[0], :].X.mean(0)
        res.loc[clust] = X[np.where(cluster_label == clust)[0], :].X.sum(0)
    return res

# =======================================

def pseudoPBMCSim():
    exp_types = ["pbmc1", "pbmc2"]
    protocol_types = ["Drop", "inDrops"]
    num_genes = ["50hvg", "100hvg"]
    sim_name = "NORTA"
    data_path = "./simulated/new/"
    # save_path = "./pseudo_bulk/"
    save_path = "./pseudo_bulk_sum/"
    data_names = itertools.product(exp_types, protocol_types, num_genes)
    for each in data_names:
        data_name = "-".join(each)
        print("[ {} ]".format(data_name))
        expr_mat = loadData(data_path, data_name, sim_name)
        print("Sparsity (raw) = {} | # cells = {}".format(np.count_nonzero(expr_mat) / np.prod(expr_mat.shape), expr_mat.shape[0]))
        # -----
        pseudo_data = pseudoBulk(expr_mat)
        print("Sparsity (pseudo-bulk) = {} | # cells = {}".format(np.count_nonzero(pseudo_data) / np.prod(pseudo_data.shape), pseudo_data.shape[0]))
        pseudo_data.to_csv("{}/{}-{}-data_mat.csv".format(save_path, data_name, sim_name))


def pseudoCortexSim():
    exp_types = ["Cortex1", "Cortex2"]
    protocol_types = ["10xChromium", "Smart_seq2"]
    num_genes = ["50hvg", "100hvg"]
    sim_name = "NORTA"
    data_path = "./simulated/new/"
    save_path = "./pseudo_bulk_sum/"
    data_names = itertools.product(exp_types, protocol_types, num_genes)
    for each in data_names:
        data_name = "-".join(each)
        print("[ {} ]".format(data_name))
        expr_mat = loadData(data_path, data_name, sim_name)
        print("Sparsity (raw) = {} | # cells = {}".format(np.count_nonzero(expr_mat) / np.prod(expr_mat.shape), expr_mat.shape[0]))
        # -----
        pseudo_data = pseudoBulk(expr_mat)
        print("Sparsity (pseudo-bulk) = {} | # cells = {}".format(np.count_nonzero(pseudo_data) / np.prod(pseudo_data.shape), pseudo_data.shape[0]))
        pseudo_data.to_csv("{}/{}-{}-data_mat.csv".format(save_path, data_name, sim_name))


def pseudoSERGIOSim():
    data_path = "./SERGIO_simulation_all/"
    save_path = "./SERGIO_pseudobulk/"
    sparsity_list = [1, 5, 10, 15, 20]
    data_names = itertools.product(["100gene-9groups"], sparsity_list)
    for each in data_names:
        data_name = each
        print("[ {} ]".format(data_name))
        expr_mat = loadSERGIOData(data_path, data_name)
        print("Sparsity (raw) = {} | # cells = {}".format(np.count_nonzero(expr_mat) / np.prod(expr_mat.shape), expr_mat.shape[0]))
        # -----
        pseudo_data = pseudoBulk(expr_mat)
        print("Sparsity (pseudo-bulk) = {} | # cells = {}".format(np.count_nonzero(pseudo_data) / np.prod(pseudo_data.shape), pseudo_data.shape[0]))
        pseudo_data.to_csv("{}/{}-{}.csv".format(save_path, data_name[0], data_name[1]))


if __name__ == '__main__':
    # pseudoPBMCSim()
    # pseudoCortexSim()
    pseudoSERGIOSim()