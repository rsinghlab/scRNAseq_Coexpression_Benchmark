'''
Description:
    MAGIC imputation.

Author:
    Jiaqi Zhang <jiaqi_zhang2@brown.edu>
'''
import itertools
import magic
import pandas as pd
import numpy as np
import scanpy
from SimulationUtils import loadData

# ======================================

def MAGICImputation(X):
    magic_operator = magic.MAGIC(verbose=False)
    X_magic = magic_operator.fit_transform(X)
    return X_magic

# =======================================
#TODO: file path

def imputePBMCSim():
    exp_types = ["pbmc1", "pbmc2"]
    protocol_types = ["Drop", "inDrops"]
    num_genes = ["100hvg"]
    sim_name = "NORTA"
    data_path = "./simulated/new/"
    save_path = "./imputed_simulation/"
    data_names = itertools.product(exp_types, protocol_types, num_genes)
    for each in data_names:
        data_name = "-".join(each)
        print("[ {} ]".format(data_name))
        expr_mat = loadData("{}/{}-{}-data_mat.csv".format(data_path, data_name, sim_name))
        print("Sparsity (before) = {}".format(np.count_nonzero(expr_mat) / np.prod(expr_mat.shape)))
        # -----
        imputed_data = MAGICImputation(expr_mat)
        print("Sparsity (after MAGIC) = {}".format(np.count_nonzero(imputed_data) / np.prod(imputed_data.shape)))
        imputed_data.to_csv("{}/{}-{}-MAGIC-data_mat.csv".format(save_path, data_name, sim_name))


def imputeCortexSim():
    exp_types = ["Cortex1", "Cortex2"]
    protocol_types = ["10xChromium", "Smart_seq2"]
    num_genes = ["100hvg"]
    sim_name = "NORTA"
    data_path = "./simulated/new/"
    save_path = "./imputed_simulation/"
    data_names = itertools.product(exp_types, protocol_types, num_genes)
    for each in data_names:
        data_name = "-".join(each)
        print("[ {} ]".format(data_name))
        expr_mat = loadData("{}/{}-{}-data_mat.csv".format(data_path, data_name, sim_name))
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
        expr_mat = loadData("{}/{}-{}sparsity.csv".format(data_path, data_name[0], data_name[1]))
        print("Sparsity (before) = {}".format(np.count_nonzero(expr_mat) / np.prod(expr_mat.shape)))
        # -----
        imputed_data = MAGICImputation(expr_mat)
        print("Sparsity (after MAGIC) = {}".format(np.count_nonzero(imputed_data) / np.prod(imputed_data.shape)))
        imputed_data.to_csv("{}/{}-{}-MAGIC-data_mat.csv".format(save_path, data_name[0], data_name[1]))


if __name__ == '__main__':
    imputePBMCSim()
    imputeCortexSim()
    imputeSERGIOSim()