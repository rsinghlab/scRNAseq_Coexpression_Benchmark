import torch
import pandas as pd
import numpy as np
import itertools
import scipy.stats
import os
from sklearn.metrics import f1_score, roc_auc_score, precision_recall_curve, auc
from util.BenchmarkUtils import loadExprMat, loadNetMat, loadMetric
from util.Visualization import visBox4Metric, visTprFpr, visMetric4Imputation, compareGraphImputation

# ======================================

def loadMetric4Impute(metric_dir_path, data_name, sim_name, num_genes, n):
    raw_metric_list = [each for each in os.listdir(metric_dir_path)
                   if (data_name in each) and (sim_name in each) and ("process-metric" not in each) and (num_genes in each)
                       and ("n{}".format(n) in each) and ("netSmooth" in each)]
    metric_dict = {}
    metric_dict["raw"] = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(metric_dir_path, each)) for each in raw_metric_list}
    return metric_dict

def loadMetric4WoImpute(metric_dir_path, data_name, sim_name, num_genes):
    raw_metric_list = [each for each in os.listdir(metric_dir_path)
                   if (data_name in each) and (sim_name in each) and ("process-metric" not in each) and (num_genes in each)]
    metric_dict = {}
    metric_dict["raw"] = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(metric_dir_path, each)) for each in raw_metric_list}
    return metric_dict


def loadSimData(data_dir_path, data_name, protocal_name, sim_name, num_genes, n):
    data_filename = "{}/{}-{}-{}-{}-n{}-data_mat.csv".format(data_dir_path, data_name, protocal_name, num_genes, sim_name, n)
    expr_mat = loadExprMat(data_filename)
    return expr_mat

def loadNet(data_dir_path, data_name, protocal_name, sim_name, num_genes, n):
    data_filename = "{}/{}-{}-{}-n{}-net_mat.csv".format(data_dir_path, data_name, protocal_name, num_genes, n)
    expr_mat = loadExprMat(data_filename)
    return expr_mat


def loadImputedData(data_dir_path, data_name, protocal_name, sim_name, num_genes, impute_name, n):
    data_filename = "{}/{}-{}-{}-{}-{}-n{}-data_mat.csv".format(data_dir_path, data_name, protocal_name, num_genes, sim_name, impute_name, n)
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
    model_metrics_list = [metric_dict.loc[m] for m in all_models_list]
    return model_metrics_list


def process4TprFprVis(metric_dict, random_metric_dict):
    model_tpr_list = [[metric_dict[each][0].loc[m] if m in metric_dict[each][0] else np.nan for each in metric_dict] for m in all_models_list]
    model_fpr_list = [[metric_dict[each][1].loc[m] if m in metric_dict[each][0] else np.nan for each in metric_dict] for m in all_models_list]
    model_tpr_list.append([each[0] for each in list(random_metric_dict.values())])
    model_fpr_list.append([each[1] for each in list(random_metric_dict.values())])
    return model_tpr_list, model_fpr_list

# ======================================

# all_models_list = ["pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet"]
all_models_list = ["pearson", "spearman", "glasso", "scLink", "GENIE3", "scDesign2", "minet"]


def compareModelMetricsOnImputation(metric_dir_path, wo_impute_metric_dir_path, data_name):
    # data_name = "Cortex1-Smart_seq2-"
    data_name = "pbmc1-Drop-"
    sim_name = "NORTA"
    num_genes = "50hvg"
    # -----
    wo_impute_metric_dict = loadMetric4WoImpute(wo_impute_metric_dir_path, data_name=data_name, sim_name=sim_name, num_genes=num_genes)
    wo_impute_metric_df = wo_impute_metric_dict["raw"][data_name + num_genes]
    wo_impute_best_F1 = getModelPerformance(wo_impute_metric_df, "F1", func="max")
    wo_impute_F1_vis = process4Vis(wo_impute_best_F1)
    # ----
    impute_metric_list = []
    for n in ["full", "10", "20", "30"]:
        print("=" * 70)
        print("[ {} | {} ] Loading metrics...".format(sim_name, num_genes))
        impute_metric_dict = loadMetric4Impute(metric_dir_path, data_name=data_name, sim_name=sim_name, num_genes=num_genes, n=n)
        impute_metric_df = impute_metric_dict["raw"][data_name+num_genes]
        impute_best_F1 = getModelPerformance(impute_metric_df, "F1", func="max")
        impute_F1_vis = process4Vis(impute_best_F1)
        impute_metric_list.append(impute_F1_vis)
    compareGraphImputation(wo_impute_F1_vis, impute_metric_list, all_models_list,
                           edge_list=["full", "10", "20", "30"], title="{}{}-graph_based".format(data_name, num_genes),
                           metric_name="F1-score")



def compareCorrMetricsOnImputation(raw_data_dir, imputed_data_dir, data_name):
    protocol_list = ["Smart_seq2", "10xChromium"] if data_name == "Cortex" else ["Drop", "inDrops"]
    sim_name = "NORTA"
    data_name_list = []
    cmd_list = []
    for num_genes in ["50hvg", "100hvg"]:
        for protocol_name in protocol_list:
            for exp_name in ["1", "2"]:
                print("=" * 70)
                print("[ {} | {} ] Loading metrics...".format(sim_name, num_genes))
                raw_data = loadSimData(raw_data_dir, data_name + exp_name, protocol_name, "NORTA", num_genes)
                true_net = loadNet(raw_data_dir, data_name + exp_name, protocol_name, "NORTA", num_genes)
                saver_imputed_data = loadImputedData(imputed_data_dir, data_name + exp_name, protocol_name, "NORTA", num_genes, "SAVER")
                magic_imputed_data = loadImputedData(imputed_data_dir, data_name + exp_name, protocol_name, "NORTA", num_genes, "MAGIC")
                raw_corr = computeCorr(raw_data)["PCC"]
                saver_corr = computeCorr(saver_imputed_data)["PCC"]
                magic_corr = computeCorr(magic_imputed_data)["PCC"]
                _, _, (saver_before_CMD, saver_after_CMD), (saver_before_MSE, saver_after_MSE) = computeCorrSim(raw_corr, saver_corr, true_net)
                _, _, (magic_before_CMD, magic_after_CMD), (magic_before_MSE, magic_after_MSE) = computeCorrSim(magic_corr, saver_corr, true_net)
                print() #TODO: 看起来impute之后CMD减小了





if __name__ == '__main__':
    metric_dir_path = "../res/logs/evaluation/graph_imputed_simulation/"
    wo_imputation_metric_dir_path = "../res/logs/evaluation/simulation/"

    imputed_data_dir = "../data/graph_imputed_simulation/"
    raw_data_dir = "../data/graph_imputed_simulation/"
    compareModelMetricsOnImputation(metric_dir_path, wo_imputation_metric_dir_path, data_name="Cortex")

    # compareCorrMetricsOnImputation(raw_data_dir, imputed_data_dir, data_name="Cortex")
    # compareCorrMetricsOnImputation(raw_data_dir, imputed_data_dir, data_name="pbmc")
    pass