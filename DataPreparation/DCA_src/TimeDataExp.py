'''
Description:
    Model running for PBMC data.

Author:
    Jiaqi Zhang
'''
import sys
sys.path.append("/")
sys.path.append("../")
from api import *
import numpy as np
import scanpy
from scipy.sparse import coo_matrix
from scipy.io import mmwrite



train_set = np.load("../../Data/dyn_true.npy")
print("Data shape : {}".format(train_set.shape))
train_data = scanpy.AnnData(train_set) if not isinstance(train_set, scanpy.AnnData) else train_set
filters = scanpy.pp.filter_genes(train_data, min_counts=1, inplace=False)
filtered_train_data = train_data[:, filters[0]] # cells x genes

denoised_data, dca_model = dca(
        filtered_train_data,
        return_model=True,
        copy=True,
        normalize_per_cell=False,
        scale=False,
        log1p=False
)
print("Denoised data shape : {}".format(denoised_data))
sparse_predictions = denoised_data.to_df().values
np.save("../../Data/denoised_dyn.npy", sparse_predictions)
np.save("../../Data/dyn_true_filter.npy", filters[0])



