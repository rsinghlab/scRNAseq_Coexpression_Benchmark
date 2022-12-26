'''
Description:
    Train and evaluate the DCA model.

Authot:
    Jiaqi Zhang
'''
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "1"

import sys
sys.path.append("/")
sys.path.append("../")
sys.path.append("./dca/")
from api import *
import numpy as np
import scanpy
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
import pickle as pkl
# from Models.Py_Utils import readMtx, readH5ad
from Py_Utils import readMtx, readH5ad


def normalizeData(data):
    '''
    Description:
        Data normalization.
    '''
    data = scanpy.AnnData(data) if not isinstance(data, scanpy.AnnData) else data
    scanpy.pp.filter_genes(data, min_counts=1)
    scanpy.pp.normalize_per_cell(data)
    scanpy.pp.log1p(data)
    return data


def train_model(train_set):
    '''
    Description:
        Train DCA model.
    '''
    train_data = scanpy.AnnData(train_set) if not isinstance(train_set, scanpy.AnnData) else train_set
    denoised_data, dca_model = dca(
        train_data,
        return_model=True,
        copy=True,
        normalize_per_cell=True,
        scale=True,
        log1p=True
    )
    return denoised_data, dca_model


def eval_model(model, sample_set):
    sample_set = scanpy.AnnData(sample_set) if not isinstance(sample_set, scanpy.AnnData) else sample_set
    sample_set = sample_set.copy()
    scanpy.pp.normalize_per_cell(sample_set)
    sample_set.obs['size_factors'] = sample_set.obs.n_counts / np.median(sample_set.obs.n_counts)
    generated_data = model.predict(sample_set, return_info=False, copy=True)
    return generated_data

# ------------------------------------

def clusterAugmentation(config):
    # Parameters configuration
    need_save = config["need_save"]
    train_size = config["train_size"]
    # Prepare datasets
    print("=" * 70)
    print("START LOADING DATA...")
    if config["train_data"].split(".")[-1] == "mtx":
        train_data = readMtx(config["train_data"])  # gene by cell
        train_data = scanpy.AnnData(train_data)
    elif config["train_data"].split(".")[-1] == "h5ad":
        train_data = readH5ad(config["train_data"])
        train_data = scanpy.AnnData(train_data)
    else:
        raise ValueError("Files of {} format cannot be processed!".format(config["train_data"].split(".")[-1]))
    # -----
    # Split data into two parts
    num_cells = train_data.shape[0]
    sampled_idx = np.random.choice(np.arange(num_cells), int(num_cells*train_size), replace=False)
    sampled_data = train_data[sampled_idx, :]
    removed_data = train_data[~sampled_idx, :]
    # Train the model
    print("=" * 70)
    print("START TRAINING...")
    filters = scanpy.pp.filter_genes(sampled_data, min_counts=1, inplace=False)
    filtered_sampled_data = sampled_data[:, filters[0]]  # cells x genes
    filtered_removed_data = removed_data[:, filters[0]]  # cells x genes
    print("Num of cells = {} | Num of genes = {} ".format(filtered_sampled_data.shape[0], filtered_sampled_data.shape[1]))
    filtered_sampled_set = scanpy.AnnData(filtered_sampled_data) if not isinstance(filtered_sampled_data, scanpy.AnnData) else filtered_sampled_data
    _, dca_model = train_model(filtered_sampled_set)
    # Data simulation with trained model
    print("=" * 70)
    print("START GENERATING...")
    num_new_cells = filtered_removed_data.shape[0]
    for_prediction_idx = np.random.choice(np.arange(filtered_sampled_data.shape[0]), num_new_cells, replace=True)
    for_prediction_data = filtered_sampled_data[for_prediction_idx]
    generated_data = eval_model(dca_model, for_prediction_data)
    print("Generated data shape : {}".format(generated_data))
    # Save records
    summary = {}
    summary["generated_data"] = generated_data.to_df().values
    summary["sampled_data"] = filtered_sampled_data.to_df().values
    summary["removed_data"] = filtered_removed_data.to_df().values
    summary["filters"] = filters[0]
    if need_save:
        np.save(config["prediction_save_path"], summary)
        print("Finished saving records.")
    return summary


def normalAugmentation(config):
    need_save = config["need_save"]
    # Prepare datasets
    print("=" * 70)
    print("START LOADING DATA...")
    if config["train_data"].split(".")[-1] == "mtx":
        train_data = readMtx(config["train_data"])  # cell by gene
        train_data = scanpy.AnnData(train_data)
    elif config["train_data"].split(".")[-1] == "h5ad":
        train_data = readH5ad(config["train_data"])  # cell by gene
        train_data = scanpy.AnnData(train_data)
    else:
        raise ValueError("Files of {} format cannot be processed!".format(config["train_data"].split(".")[-1]))
    # train_data = train_data[:50, :20]  # TODO: SERGIO_simulation_all
    filters = scanpy.pp.filter_genes(train_data, min_counts=1, inplace=False)
    filtered_train_data = train_data[:, filters[0]] # cells x genes

    print("Num of cells = {} | Num of genes = {} ".format(filtered_train_data.shape[0], filtered_train_data.shape[1]))
    # -----
    # Train the model (DCA has no SERGIO_simulation_all)
    print("=" * 70)
    print("START TRAINING...")
    denoised_data, dca_model = train_model(filtered_train_data)
    #TODO save model
    # dca_model.save(config["model_save_path"])
    # with open(config["model_save_path"], "wb") as file:
    #     pkl.dump(dca_model, file)
    summary = {}
    summary["denoised_data"] = coo_matrix(denoised_data.to_df().values)
    summary["train_data"] = coo_matrix(train_data.to_df().values)
    summary["filters"] = filters
    if need_save:
        # np.savez_compressed(config["prediction_save_path"], **summary)
        mmwrite(config["prediction_save_path"], summary["denoised_data"])
        mmwrite(config["train_data_save_path"], summary["train_data"])
        np.save(config["filter_save_path"], summary["filters"])
        print("Finished saving records.")
    return summary


# ------------------------------------

def train_model_w_par(train_set, dropout, layer_size):
    '''
    Description:
        Train DCA model.
    '''
    train_data = scanpy.AnnData(train_set) if not isinstance(train_set, scanpy.AnnData) else train_set
    denoised_data, dca_model = dca(
        train_data,
        return_model=True,
        copy=True,
        normalize_per_cell=True,
        scale=True,
        log1p=True,
        hidden_size=layer_size,  # network args
        hidden_dropout=dropout,
        epochs=10,
    )
    return denoised_data, dca_model


def modelTrainForCV(config):
    filtered_train_data = config["train_data"]  # cell by gene
    filtered_validate_data = config["validate_data"]  # cell by gene
    # -----
    # Train the model (DCA has no SERGIO_simulation_all)
    print("=" * 70)
    print("START TRAINING...")
    denoised_data, dca_model = train_model_w_par(filtered_train_data, config["dropout"], config["layer_size"])
    return filtered_validate_data, denoised_data


if __name__ == '__main__':
    from scipy.stats import pearsonr, spearmanr
    import matplotlib.pyplot as plt
    import umap
    # --
    config = {
        "train_data": "../../Data/splat_simulation/wo_preprocess/testing_all.mtx",
        # ------------------------------------------
        "need_save": False,
        "train_size": 0.5,
        # ------------------------------------------
        # "model_save_path": "../../Prediction/DCA/PBMC-DCA_model.h5",
        # "prediction_save_path": "../../Prediction/DCA/PBMC-DCA_estimation.npy",
    }
    # res = normalAugmentation(config)
    res = normalAugmentation(config)
    denoised_data = res["denoised_data"]
    train_data = res["train_data"]
    filters = res["filters"]
    # ---
    mean_gene_exp_pred = np.mean(denoised_data, axis=0)
    mean_gene_exp_labels = np.mean(train_data, axis=0)
    test_pcc, _ = pearsonr(mean_gene_exp_pred, mean_gene_exp_labels)
    test_scc, _ = spearmanr(mean_gene_exp_pred, mean_gene_exp_labels)
    print('Test PCC:', test_pcc)
    print('Test SCC:', test_scc)

    test_pcc_str = str(np.round(test_pcc, 3))
    test_scc_str = str(np.round(test_scc, 3))
    cc_str = (
        f"Test PCC = {test_pcc_str}\n"
        f"Test SCC = {test_scc_str}"
    )

    scores_and_labels = np.concatenate((train_data, denoised_data), axis=0)
    model = umap.UMAP().fit(train_data)
    emedding = model.transform(scores_and_labels)
    plt.scatter(emedding[:1500, 0], emedding[:1500, 1], color='blue', alpha=0.5, s=5, label='True')
    plt.scatter(emedding[1500:, 0], emedding[1500:, 1], color='orange', alpha=0.5, s=5, label='Simulated')
    plt.legend()
    plt.title('UMAP', fontsize=15)
    plt.xticks([])
    plt.yticks([])
    plt.xlabel('Dimension 1', fontsize=15)
    plt.ylabel('Dimension 2', fontsize=15)
    plt.show()
