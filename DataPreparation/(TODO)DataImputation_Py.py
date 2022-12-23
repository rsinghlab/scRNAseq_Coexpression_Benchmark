import itertools
import magic
import pandas as pd
import numpy as np
import scanpy

# import sys
# sys.path.append("../util/dca/")
# from api import *


def loadData(data_path, data_name, sim_name):
    filename = "{}/{}-{}-data_mat.csv".format(data_path, data_name, sim_name)
    expr_mat = pd.read_csv(filename, header=0, index_col=0)
    return expr_mat


def loadSERGIOData(data_path, data_name):
    filename = "{}/{}-{}sparsity.csv".format(data_path, data_name[0], data_name[1])
    expr_mat = pd.read_csv(filename, header=0, index_col=0)
    return expr_mat


def MAGICImputation(X):
    magic_operator = magic.MAGIC(verbose=False)
    X_magic = magic_operator.fit_transform(X)
    # output_data = X_magic.round()
    return X_magic

# =======================================

def imputePBMCSim():
    exp_types = ["pbmc1", "pbmc2"]
    protocol_types = ["Drop", "inDrops"]
    num_genes = ["50hvg", "100hvg"]
    sim_name = "NORTA"
    data_path = "./simulated/new/"
    save_path = "./imputed_simulation/"
    data_names = itertools.product(exp_types, protocol_types, num_genes)
    for each in data_names:
        data_name = "-".join(each)
        print("[ {} ]".format(data_name))
        expr_mat = loadData(data_path, data_name, sim_name)
        print("Sparsity (before) = {}".format(np.count_nonzero(expr_mat) / np.prod(expr_mat.shape)))
        # -----
        imputed_data = MAGICImputation(expr_mat)
        print("Sparsity (after MAGIC) = {}".format(np.count_nonzero(imputed_data) / np.prod(imputed_data.shape)))
        imputed_data.to_csv("{}/{}-{}-MAGIC-data_mat.csv".format(save_path, data_name, sim_name))
        # -----
        # imputed_data = DCAImputation(expr_mat)
        # print("Sparsity (after DCA) = {}".format(np.count_nonzero(imputed_data) / np.prod(imputed_data.shape)))
        # imputed_data.to_csv("{}/{}-{}-DCA-data_mat.csv".format(save_path, data_name, sim_name))


def imputeCortexSim():
    exp_types = ["Cortex1", "Cortex2"]
    protocol_types = ["10xChromium", "Smart_seq2"]
    num_genes = ["50hvg", "100hvg"]
    sim_name = "NORTA"
    data_path = "./simulated/new/"
    save_path = "./imputed_simulation/"
    data_names = itertools.product(exp_types, protocol_types, num_genes)
    for each in data_names:
        data_name = "-".join(each)
        print("[ {} ]".format(data_name))
        expr_mat = loadData(data_path, data_name, sim_name)
        print("Sparsity (before) = {}".format(np.count_nonzero(expr_mat) / np.prod(expr_mat.shape)))
        imputed_data = MAGICImputation(expr_mat)
        print("Sparsity (after) = {}".format(np.count_nonzero(imputed_data) / np.prod(imputed_data.shape)))
        imputed_data.to_csv("{}/{}-{}-MAGIC-data_mat.csv".format(save_path, data_name, sim_name))


def imputeSERGIOSim():
    sparsity_list = [1, 5, 10, 15, 20]
    data_names = itertools.product(["100gene-9groups"], sparsity_list)
    data_path = "./SERGIO_simulation_all"
    save_path = "./SERGIO_imputation/"
    for each in data_names:
        data_name =each
        print("[ {} ]".format(data_name))
        expr_mat = loadSERGIOData(data_path, data_name)
        print("Sparsity (before) = {}".format(np.count_nonzero(expr_mat) / np.prod(expr_mat.shape)))
        # -----
        imputed_data = MAGICImputation(expr_mat)
        print("Sparsity (after MAGIC) = {}".format(np.count_nonzero(imputed_data) / np.prod(imputed_data.shape)))
        imputed_data.to_csv("{}/{}-{}-MAGIC-data_mat.csv".format(save_path, data_name[0], data_name[1]))


if __name__ == '__main__':
    # imputePBMCSim()
    # imputeCortexSim()
    imputeSERGIOSim()