'''
Description:
    DCA imputation.

Author:
    Jiaqi Zhang <jiaqi_zhang2@brown.edu>
'''
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "1"

import itertools
import pandas as pd
import numpy as np
import scanpy

import sys
sys.path.append("../util/dca/") #TODO: dca
from api import *
from SimulationUtils import loadData

# =======================================

def dcaImputation(data_mat):
    original_shape = data_mat.shape
    original_columns = data_mat.columns
    gene_sum = np.sum(data_mat, axis=0)
    non_zero_idx = np.where(gene_sum > 0)[0]
    data_mat = data_mat.iloc[:, non_zero_idx]
    train_data = scanpy.AnnData(data_mat) if not isinstance(data_mat, scanpy.AnnData) else data_mat
    denoised_data, dca_model = dca(
        train_data,
        return_model=True,
        copy=True,
        normalize_per_cell=True,
        scale=True,
        log1p=True
    )
    denoised_data = denoised_data.to_df().values
    imputed_data = np.zeros((denoised_data.shape[0], original_shape[1]))
    imputed_data[:, non_zero_idx] = denoised_data
    imputed_data = pd.DataFrame(data=imputed_data, index=data_mat.index[:denoised_data.shape[0]], columns=original_columns)
    return imputed_data

# =======================================

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
        imputed_data = dcaImputation(expr_mat)
        print("Sparsity (after DCA) = {}".format(np.count_nonzero(imputed_data) / np.prod(imputed_data.shape)))
        imputed_data.to_csv("{}/{}-{}-DCA-data_mat.csv".format(save_path, data_name, sim_name))


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
        # -----
        imputed_data = dcaImputation(expr_mat)
        print("Sparsity (after DCA) = {}".format(np.count_nonzero(imputed_data) / np.prod(imputed_data.shape)))
        imputed_data.to_csv("{}/{}-{}-DCA-data_mat.csv".format(save_path, data_name, sim_name))


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
        imputed_data = dcaImputation(expr_mat)
        print("Sparsity (after DCA) = {}".format(np.count_nonzero(imputed_data) / np.prod(imputed_data.shape)))
        imputed_data.to_csv("{}/{}-{}-DCA-data_mat.csv".format(save_path, data_name[0], data_name[1]))


if __name__ == '__main__':
    #TODO: file path
    imputePBMCSim()
    imputeCortexSim()
    imputeSERGIOSim()
