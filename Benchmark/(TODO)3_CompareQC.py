import scanpy
import torch
import pandas as pd
import numpy as np
import itertools
import scipy.stats
import os
from util.BenchmarkUtils import loadExprMat, loadNetMat, loadMetric
from util.Visualization import visMetric4QC



def loadSimData(data_dir_path, data_name, sim_name, num_genes):
    data_filename = "{}/{}-{}-data_mat.csv".format(data_dir_path, data_name, sim_name)
    expr_mat = loadExprMat(data_filename)
    return expr_mat


def loadQCData(data_dir_path, data_name, sim_name, num_genes, percentile):
    data_filename = "{}/{}-{}-{}percentile-data_mat.csv".format(data_dir_path, data_name, sim_name, percentile)
    expr_mat = loadExprMat(data_filename)
    return expr_mat


# ======================================

def dataSparsity(data):
    data_vec = (data.reshape(-1) != 0).astype(int)
    sparsity = np.mean(data_vec)
    return sparsity


def percentileLibSizeList(data, p_list):
    library_size = np.sum(data, axis=1)
    # percent_size = np.asarray([int(np.percentile(library_size, p)) for p in p_list])
    percent_size = np.asarray([int(np.max(library_size) * p * 0.01) for p in p_list])
    return percent_size


def qcExp(data_dir_path, data_name):
    # protocol_list = ["Smart_seq2", "10xChromium"] if data_name == "Cortex" else ["Drop", "inDrops"]
    protocol_list = ["Smart_seq2", "10xChromium"] if data_name == "Cortex" else ["Drop", "inDrops"]
    for num_genes in ["50hvg", "100hvg"]:
        for exp_name in ["1", "2"]:
            for protocol_name in protocol_list:
                print("=" * 70)
                print("[ {} | {} | {} ] Loading metrics...".format(data_name+exp_name, protocol_name, num_genes))
                # expr_mat = loadRealData(data_dir_path, data_name+exp_name, protocol_name, num_genes)
                expr_mat = loadSimData(data_dir_path, data_name+exp_name, protocol_name, "NORTA", num_genes)
                expr_adata = scanpy.AnnData(expr_mat)
                print("Data : sparsity = {} | # cells = {}".format(dataSparsity(expr_adata.to_df().values), expr_adata.shape[0]))
                p_list = [1, 5, 10, 15, 20]
                qc_cnt_list = percentileLibSizeList(expr_adata.to_df().values, p_list)
                for i, x in enumerate(qc_cnt_list):
                    filter_adata = expr_adata[scanpy.pp.filter_cells(expr_adata, inplace=False, min_counts=x)[0]]
                    print("Data (min count = {}) : sparsity = {} | # cells = {}".format(x, dataSparsity(filter_adata.to_df().values), filter_adata.shape[0]))
                    # sace data
                    filter_adata.to_df().to_csv("../data/quality_control/{}-{}-{}-{}-{}percentile-data_mat.csv".format(data_name+exp_name, protocol_name, num_genes, "NORTA", p_list[i]))

# ======================================

def load4DiffSetting(metric_dir_path, data_name, sim_name, num_genes, percentile):
    raw_metric_list = [each for each in os.listdir(metric_dir_path)
                   if (data_name in each) and (sim_name in each) and ("process-metric" not in each)
                       and (num_genes in each) and ("-{}percentil".format(percentile) in each)]
    metric_dict = {}
    metric_dict["raw"] = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(metric_dir_path, each)) for each in raw_metric_list}
    return metric_dict


def loadMEtric4RawSim(metric_dir_path, data_name, sim_name, num_genes):
    raw_metric_list = [each for each in os.listdir(metric_dir_path)
                   if (data_name in each) and (sim_name in each) and ("process-metric" not in each)
                       and (num_genes in each)]
    metric_dict = {}
    metric_dict["raw"] = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(metric_dir_path, each)) for each in raw_metric_list}
    return metric_dict

# ======================================

def getModelPerformance(metric_df, metric_name, func="max"):
    res_df = metric_df.groupby("Model")[metric_name]
    if func == "max":
        return res_df.max()#.reset_index()
    elif func == "min":
        return res_df.min().reset_index()
    elif func == "mean":
        return res_df.mean().reset_index()
    else:
        raise ValueError("Unknown function {}!".format(func))


def process4Vis(metric_dict):
    model_metrics_list = [[metric_dict[each].loc[m] for each in metric_dict] for m in all_models_list]
    return model_metrics_list

# ======================================

all_models_list = ["pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet"]

def _dataSparsity(data):
    data = data.values
    sign_data = (data != 0).astype(int)
    sparsity =np.mean(sign_data)
    return sparsity


def compareModelMetrics4QC(metric_dir_path, wo_qc_metric_dir_path, qc_data_path, raw_data_path, data_name):
    #TODO: 只能处理一个dataset
    num_cell_list = []
    sparsity_list = []
    percent_list = [1, 5, 10, 15, 20]
    metric_list = []
    for num_genes in ["50hvg", "100hvg"]:
        for sim_name in ["NORTA"]:
            for p in percent_list:
                print("=" * 70)
                print("[ {} | {} ] Loading metrics...".format(sim_name, num_genes))
                qc_metric_dict = load4DiffSetting(metric_dir_path, data_name=data_name, sim_name=sim_name,
                                                  num_genes=num_genes, percentile=p)["raw"]
                wo_qc_metric_dict = loadMEtric4RawSim(wo_qc_metric_dir_path, data_name=data_name, sim_name=sim_name,
                                                     num_genes=num_genes)["raw"]
                wo_qc_metric_dict = {each:wo_qc_metric_dict[each] for each in qc_metric_dict}
                raw_data = {each: loadSimData(raw_data_path, each, sim_name, num_genes) for each in qc_metric_dict}
                qc_data = {each: loadQCData(qc_data_path, each, sim_name, num_genes, p) for each in qc_metric_dict}
                # -----
                if len(num_cell_list) == 0:
                    num_cell_list.append([raw_data[each].shape[0] for each in raw_data])
                    sparsity_list.append([_dataSparsity(raw_data[each]) for each in raw_data])
                num_cell_list.append([qc_data[each].shape[0] for each in qc_data])
                sparsity_list.append([_dataSparsity(qc_data[each]) for each in qc_data])
                # Compare F1 score of models
                qc_best_F1 = {each: getModelPerformance(qc_metric_dict[each], "F1", func="max") for each in qc_metric_dict}
                qc_F1_vis = process4Vis(qc_best_F1)
                if len(metric_list) == 0:
                    wo_qc_best_F1 = {each: getModelPerformance(wo_qc_metric_dict[each], "F1", func="max") for each in wo_qc_metric_dict}
                    wo_qc_F1_vis = process4Vis(wo_qc_best_F1)
                    metric_list.append(wo_qc_F1_vis)
                metric_list.append(qc_F1_vis)
        visMetric4QC(metric_list, num_cell_list, sparsity_list, percent_list, all_models_list,
                     title="{}-{}-{}-imputation".format(data_name, sim_name, num_genes), metric_name="F1 score")


if __name__ == '__main__':
    # # Generate filtered data
    # data_name = "pbmc" # Cortex, pbmc
    # # data_dir_path = "../data/experimental/{}/processed/expr/".format("mouse_cortex" if data_name == "Cortex" else "PBMC")
    # data_dir_path = "../data/simulated/new/"
    # qcExp(data_dir_path, data_name)
    # # -----
    # Compare model performance
    metric_dir_path = "../res/logs/evaluation/quality_control/"
    wo_qc_metric_dir_path = "../res/logs/evaluation/simulation/"
    qc_data_dir_path = "../data/quality_control/"
    raw_data_dir_path = "../data/simulated/new/"
    compareModelMetrics4QC(metric_dir_path, wo_qc_metric_dir_path, qc_data_dir_path, raw_data_dir_path, data_name="Cortex")
    pass
