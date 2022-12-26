import torch
import pandas as pd
import numpy as np
import itertools
import scipy.stats
import os
from sklearn.metrics import f1_score, roc_auc_score, precision_recall_curve, auc
from util.BenchmarkUtils import loadExprMat, loadNetMat, loadMetric
from util.Visualization import compareF14DiffCellNum, compareF14DiffGraph, compareF14DiffSparsity, compareCorr4DiffSparsity

# ======================================

def load4DiffSetting(metric_dir_path, data_name, setting_name):
    raw_metric_list = [each for each in os.listdir(metric_dir_path)
                   if (data_name in each) and (setting_name in each) and ("process-metric" not in each)]
    metric_dict = {}
    metric_dict["raw"] = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(metric_dir_path, each)) for each in raw_metric_list}
    return metric_dict


def loadSimData(data_dir_path, data_name, protocal_name, sim_name, num_genes):
    data_filename = "{}/{}-{}-{}-{}-data_mat.csv".format(data_dir_path, data_name, protocal_name, num_genes, sim_name)
    expr_mat = loadExprMat(data_filename)
    return expr_mat

def loadNet(data_dir_path, data_name, protocal_name, sim_name, num_genes):
    data_filename = "{}/{}-{}-{}-net_mat.csv".format(data_dir_path, data_name, protocal_name, num_genes)
    expr_mat = loadExprMat(data_filename)
    return expr_mat


# ======================================

def computeCorr(data_mat):
    corr_dict = {
        "PCC": np.corrcoef(data_mat.T),
        "SCC": scipy.stats.spearmanr(data_mat).correlation
    }
    return corr_dict

def _CMD(corr1, corr2):
    CMD = 1 - np.trace(corr1 @ corr2) / (np.linalg.norm(corr1, ord="fro") * np.linalg.norm(corr2, ord="fro"))
    return CMD


def _MSE(corr1, corr2):
    MSE = np.linalg.norm(corr1 - corr2, ord="fro") / np.prod(corr1.shape)
    return MSE


def computeCorrSim(corr1, corr2, ref):
    # Compute similarity score between correlation matrix
    # -----
    if isinstance(corr2, float) or isinstance(corr1, float):
        return None, None
    ref[np.isnan(ref)] = 0.0
    corr1[np.isnan(corr1)] = 0.0
    corr2[np.isnan(corr2)] = 0.0
    # Correlation matrix distance
    before_CMD = _CMD(corr1, ref)
    after_CMD = _CMD(corr2, ref)
    CD = after_CMD - before_CMD
    # Mean squared error
    before_MSE = _MSE(corr1, ref)
    after_MSE = _MSE(corr2, ref)
    diff_MSE = after_MSE - before_MSE
    return CD, diff_MSE, (before_CMD, after_CMD), (before_MSE, after_MSE)

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


def getModelTprFpr(metric_df):
    metric_df["FDR"] = metric_df.FP / (metric_df.FP + metric_df.TP)
    res_df = (
        metric_df[metric_df.groupby("Model")["F1"].transform(max) == metric_df["F1"]]
            .drop_duplicates(subset="Model", keep="first")
        [["Model", "Sensitivity", "FDR"]]
    )
    tpr_df = res_df.groupby("Model")["Sensitivity"].mean()
    fdr_df = res_df.groupby("Model")["FDR"].mean()
    return tpr_df, fdr_df


def process4Vis(metric_dict):
    model_metrics_list = [[metric_dict[each].loc[m] for each in metric_dict] for m in all_models_list]
    return model_metrics_list


def process4TprFprVis(metric_dict):
    model_tpr_list = [[metric_dict[each][0].loc[m] if m in metric_dict[each][0] else np.nan for each in metric_dict] for m in all_models_list]
    model_fpr_list = [[metric_dict[each][1].loc[m] if m in metric_dict[each][0] else np.nan for each in metric_dict] for m in all_models_list]
    # model_tpr_list.append([each[0] for each in list(random_metric_dict.values())])
    # model_fpr_list.append([each[1] for each in list(random_metric_dict.values())])
    return model_tpr_list, model_fpr_list

# ======================================

# all_models_list = ["pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet"]
all_models_list = ["pearson", "spearman", "glasso", "scLink", "GENIE3", "scDesign2", "minet"]


def compareModelMetricsOnDiffCell(metric_dir_path, data_name):
    F1_list = []
    tpr_fpr_list = []
    cell_num_list = ["0.1", "0.5", "1.0", "1.5", "2.0"]
    for cell_ratio in ["0.1cell", "0.5cell", "1.0cell", "1.5cell", "2.0cell"]:
        print("=" * 70)
        print("[ {} | {} ] Loading metrics...".format(data_name, cell_ratio))
        metric_dict = load4DiffSetting(metric_dir_path, data_name=data_name, setting_name=cell_ratio)
        # -----
        best_F1 = {each: getModelPerformance(metric_dict["raw"][each], "F1", func="max") for each in
                             metric_dict["raw"]}
        F1_vis = process4Vis(best_F1)
        F1_list.append(F1_vis)
        # -----
        best_tpr_fpr = {each: getModelTprFpr(metric_dict["raw"][each]) for each in metric_dict["raw"]}
        tpr_fpr_vis = process4TprFprVis(best_tpr_fpr)
        tpr_fpr_list.append(tpr_fpr_vis)
    # -----
    F1_list = np.asarray(F1_list).squeeze()
    tpr_fpr_list = np.asarray(tpr_fpr_list).squeeze()
    compareF14DiffCellNum(F1_list, all_models_list, cell_num_list, title="{}-diff_cell_num".format(data_name))


def compareModelMetricsOnDiffGraph(metric_dir_path, data_name):
    F1_list = []
    tpr_fpr_list = []
    graph_type_list = ["scale", "power", "GRN"]
    for g_t in graph_type_list:
        print("=" * 70)
        print("[ {} | {} ] Loading metrics...".format(data_name, g_t))
        metric_dict = load4DiffSetting(metric_dir_path, data_name=data_name, setting_name=g_t)
        # -----
        best_F1 = {each: getModelPerformance(metric_dict["raw"][each], "F1", func="max") for each in metric_dict["raw"]}
        F1_vis = process4Vis(best_F1)
        F1_list.append(F1_vis)
        # -----
        best_tpr_fpr = {each: getModelTprFpr(metric_dict["raw"][each]) for each in metric_dict["raw"]}
        tpr_fpr_vis = process4TprFprVis(best_tpr_fpr)
        tpr_fpr_list.append(tpr_fpr_vis)
    # -----
    F1_list = np.asarray(F1_list).squeeze()
    tpr_fpr_list = np.asarray(tpr_fpr_list).squeeze()
    compareF14DiffGraph(F1_list, all_models_list, graph_type_list, title="{}-diffgraph".format(data_name))


def compareModelMetricsOnDiffSparsity(metric_dir_path, data_name):
    F1_list = []
    tpr_fpr_list = []
    sparsity_list = ["{:.1f}".format(i) for i in np.arange(0.1, 1.0, 0.1)]
    for s in ["{}sparsity".format(i) for i in sparsity_list]:
        print("=" * 70)
        print("[ {} | {} ] Loading metrics...".format(data_name, s))
        metric_dict = load4DiffSetting(metric_dir_path, data_name=data_name, setting_name=s)
        # -----
        best_F1 = {each: getModelPerformance(metric_dict["raw"][each], "F1", func="max") for each in metric_dict["raw"]}
        F1_vis = process4Vis(best_F1)
        F1_list.append(F1_vis)
        # -----
        best_tpr_fpr = {each: getModelTprFpr(metric_dict["raw"][each]) for each in metric_dict["raw"]}
        tpr_fpr_vis = process4TprFprVis(best_tpr_fpr)
        tpr_fpr_list.append(tpr_fpr_vis)
    # -----
    F1_list = np.asarray(F1_list).squeeze()
    tpr_fpr_list = np.asarray(tpr_fpr_list).squeeze()
    compareF14DiffSparsity(F1_list, all_models_list, sparsity_list, title="{}-diff_sparsity".format(data_name))

# ======================================

def _takeUpper(mat):
    mat = np.asarray(mat)
    return mat[np.triu_indices(mat.shape[0], k=1)]


def compareDiffSparsityData(data_dir_path, data_name):
    net_mat = loadNetMat("{}/{}-net_mat.csv".format(data_dir_path, data_name))
    sparsity_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    data_list = [loadExprMat("{}/{}-{}sparsity-NORTA-data_mat.csv".format(data_dir_path, data_name, x)) for x in sparsity_list]
    full_data = loadExprMat("{}/{}-full-NORTA-data_mat.csv".format(data_dir_path, data_name))
    data_list = data_list + [full_data]
    sparsity_list = sparsity_list + ["full"]
    # -----
    corr_list = [computeCorr(each)["PCC"] for each in data_list]
    triu_corr_list = [_takeUpper(each) for each in corr_list]
    CMD_list=  [_CMD(net_mat, each) for each in corr_list]
    compareCorr4DiffSparsity(_takeUpper(net_mat), triu_corr_list, CMD_list, sparsity_list)


if __name__ == '__main__':
    # metric_dir_path = "../res/logs/evaluation/diff_cell/"
    # compareModelMetricsOnDiffCell(metric_dir_path, data_name="Cortex1-10xChromium-100hvg")

    # # -----
    # metric_dir_path = "../res/logs/evaluation/diff_graph/"
    # compareModelMetricsOnDiffGraph(metric_dir_path, data_name="Cortex1-10xChromium-100hvg")

    # # -----
    # metric_dir_path = "../res/logs/evaluation/diff_sparsity/"
    # compareModelMetricsOnDiffSparsity(metric_dir_path, data_name="100gene-hub")

    data_dir_path = "../data/diff_sparsity/"
    compareDiffSparsityData(data_dir_path, data_name="100gene-hub")

