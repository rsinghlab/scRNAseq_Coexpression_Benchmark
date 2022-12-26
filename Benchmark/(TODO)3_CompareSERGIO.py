import torch
import pandas as pd
import numpy as np
import itertools
import scipy.stats
import os
from sklearn.metrics import f1_score, roc_auc_score, precision_recall_curve, auc
from util.BenchmarkUtils import loadExprMat, loadNetMat, loadMetric
from util.Visualization import visBox4Metric, visTprFpr

# ======================================

def load4DiffSetting(metric_dir_path, data_name):
    raw_metric_list = [each for each in os.listdir(metric_dir_path) if (data_name in each) and ("process" not in each)]
    process_metric_list = [each for each in os.listdir(metric_dir_path) if (data_name in each) and ("process" in each)]
    metric_dict = {}
    metric_dict["raw"] = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(metric_dir_path, each)) for each in raw_metric_list}
    metric_dict["process"] = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(metric_dir_path, each)) for each in process_metric_list}
    return metric_dict


def loadSimData(data_dir_path, data_name, sparsity):
    data_filename = "{}/{}-{}.csv".format(data_dir_path, data_name, sparsity)
    expr_mat = loadExprMat(data_filename)
    return expr_mat


def loadNet():
    data_filename = "../data/SERGIO_simulation_all/true_network.csv"
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
    return model_tpr_list, model_fpr_list

# ======================================

# all_models_list = ["pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet"]
all_models_list = ["pearson", "spearman", "glasso", "scLink", "GENIE3", "scDesign2", "minet"]


def compareModelMetricsOnSERGIO(metric_dir_path, data_name):
    print("=" * 70)
    print("[ {} ] Loading metrics...".format(data_name))
    metric_dict = load4DiffSetting(metric_dir_path, data_name=data_name)
    # Compare F1 score of models
    raw_best_F1 = {each: getModelPerformance(metric_dict["raw"][each], "F1", func="max") for each in metric_dict["raw"]}
    raw_F1_vis = process4Vis(raw_best_F1)
    visBox4Metric(raw_F1_vis, all_models_list, title="{}-SERGIO".format(data_name), metric_name="F1 score")
    process_best_F1 = {each: getModelPerformance(metric_dict["process"][each], "F1", func="max") for each in metric_dict["process"]}
    process_F1_vis = process4Vis(process_best_F1)
    visBox4Metric(process_F1_vis, all_models_list, title="{}-SERGIO-process".format(data_name), metric_name="F1 score")
    # Compare TPR and FPR
    raw_best_tpr_fpr = {each: getModelTprFpr(metric_dict["raw"][each]) for each in metric_dict["raw"]}
    raw_tpr_fpr_vis = process4TprFprVis(raw_best_tpr_fpr)
    processed_best_tpr_fpr = {each: getModelTprFpr(metric_dict["process"][each]) for each in metric_dict["process"]}
    processed_tpr_fpr_vis = process4TprFprVis(processed_best_tpr_fpr)
    visTprFpr(raw_tpr_fpr_vis, all_models_list, title="{}-SERGIO-raw".format(data_name))
    visTprFpr(processed_tpr_fpr_vis, all_models_list, title="{}-SERGIO-processed".format(data_name))



if __name__ == '__main__':
    metric_dir_path = "../res/logs/evaluation/SERGIO_simulation_all/"
    compareModelMetricsOnSERGIO(metric_dir_path, data_name="100gene-9group")
