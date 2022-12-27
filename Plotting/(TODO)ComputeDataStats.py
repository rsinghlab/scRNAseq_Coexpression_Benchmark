import torch
import pandas as pd
import numpy as np
import itertools
import scipy.stats
import os
from sklearn.metrics import auc
from util.BenchmarkUtils import loadExprMat, loadNetMat, loadMetric

# ======================================

def computeBasicStats(data_mat):
    # Compute library size, cell/gene zero percent, cell/gene mean/std expression
    ind_data_mat = (data_mat != 0).astype(int)
    stat_dict = {}
    stat_dict["library_size"] = np.sum(data_mat, axis=1)
    stat_dict["cell_zero_percent"] = np.mean(ind_data_mat, axis=1)
    stat_dict["gene_zero_percent"] = np.mean(ind_data_mat, axis=0)
    stat_dict["cell_mean"] = np.mean(data_mat, axis=1)
    stat_dict["gene_mean"] = np.mean(data_mat, axis=0)
    stat_dict["cell_std"] = np.std(data_mat, axis=1)
    stat_dict["gene_std"] = np.std(data_mat, axis=0)
    return stat_dict


def computeCorr(data_mat):
    corr_dict = {
        "PCC": np.corrcoef(data_mat.T),
        "SCC": scipy.stats.spearmanr(data_mat).correlation
    }
    return corr_dict

# ======================================

def _minMaxScale(data_mat):
    scaled_data_mat = (data_mat - data_mat.min()) / (data_mat.max() - data_mat.min())
    return scaled_data_mat


def computeHist(data_mat, bins):
    data_hist = np.histogram(data_mat.reshape(-1), bins=bins)[0]
    return data_hist

# ======================================


def _CMD(corr1, corr2):
    CMD = 1 - np.trace(corr1 @ corr2) / (np.linalg.norm(corr1, ord="fro") * np.linalg.norm(corr2, ord="fro"))
    return CMD


def _MSE(corr1, corr2):
    MSE = np.linalg.norm(corr1 - corr2, ord="fro") / np.prod(corr1.shape)
    return MSE


def computeCorrSim(corr_dict, type):
    # Compute similarity score between correlation matrix
    # -----
    ref_cor = corr_dict["ref"][type]
    sim_cor = corr_dict["NORTA"][type]
    true_net = corr_dict["sim_net"]
    if isinstance(ref_cor, float) or isinstance(sim_cor, float):
        return None, None
    ref_cor[np.isnan(ref_cor)] = 0.0
    sim_cor[np.isnan(sim_cor)] = 0.0
    # Correlation matrix distance
    before_CMD = _CMD(ref_cor, true_net)
    after_CMD = _CMD(sim_cor, true_net)
    CD = after_CMD - before_CMD
    # Mean squared error
    before_MSE = _MSE(ref_cor, true_net)
    after_MSE = _MSE(sim_cor, true_net)
    diff_MSE = after_MSE - before_MSE
    return CD, diff_MSE, (before_CMD, after_CMD), (before_MSE, after_MSE)

# ======================================

new_simulation_list = ["NORTA", "SERGIO"]
old_simulation_list = ["scLink", "ZILGM"]
experimental_list = ["mouse_cortex", "PBMC"]


def experimentalStats():
    data_dir_path = "../data/experimental/"
    save_dir_path = "./res/"
    num_genes_list = ["50hvg", "100hvg"]
    bins = 20
    file_stats_dict = {}
    file_hist_dict = {}
    file_corr_dict = {}
    for num_genes in num_genes_list:
        stats_dict = {}
        hist_dict = {}
        corr_dict = {}
        print("=" * 70)
        print("[ {} ]".format(num_genes))
        # Experimental data
        for data_name in experimental_list:
            if data_name == "PBMC":
                t_list = ["pbmc1-Drop", "pbmc2-Drop", "pbmc1-inDrops", "pbmc2-inDrops"]
            else:
                t_list = ["Cortex1-10xChromium", "Cortex2-10xChromium", "Cortex1-Smart_seq2", "Cortex2-Smart_seq2"]
            print("{} experimental data...".format(data_name))
            for t in t_list:
                exp_data_mat = loadExprMat("{}/{}/processed/expr/{}-{}.csv".format(data_dir_path, data_name, t, num_genes)).values
                stats_dict["{}-{}".format(t, num_genes)] = computeBasicStats(exp_data_mat)
                hist_dict["{}-{}".format(t, num_genes)] = computeHist(_minMaxScale(exp_data_mat), bins)
                corr_dict["{}-{}".format(t, num_genes)] = computeCorr(exp_data_mat)
        file_stats_dict[num_genes] = stats_dict
        file_hist_dict[num_genes] = hist_dict
        file_corr_dict[num_genes] = corr_dict
    # Save data
    np.save("{}/experimental-basic_stats.npy".format(save_dir_path), file_stats_dict)
    np.save("{}/experimental-hist.npy".format(save_dir_path), file_hist_dict)
    np.save("{}/experimental-corr.npy".format(save_dir_path), file_corr_dict)


def oldSimStats():
    sim_dir_path = "../data/simulated/old/"
    save_dir_path = "./res/"
    num_genes_list = ["50hvg", "100hvg"]
    bins = 20
    file_stats_dict = {}
    file_hist_dict = {}
    for num_genes in num_genes_list:
        stats_dict = {}
        hist_dict = {}
        print("=" * 70)
        print("[ {} ]".format(num_genes))
        # Simulated data
        for simulation_type in old_simulation_list:
            t_list = ["0.07", "0.10", "0.13", "0.16"] if simulation_type == "scLink" else ["0.8", "0.9", "1.0"]
            print("{} simulation data...".format(simulation_type))
            for t in t_list:
                if simulation_type == "scLink":
                    log_sim_data_mat = loadExprMat(
                        "{}/{}-{}_{}rho-data_mat.csv".format(sim_dir_path, num_genes, simulation_type, t)).values
                    non_zero_idx = np.where(log_sim_data_mat != 0)
                    tmp = np.power(10, log_sim_data_mat[non_zero_idx]) - 1
                    sim_data_mat = np.zeros_like(log_sim_data_mat)
                    sim_data_mat[non_zero_idx] = tmp
                else:
                    sim_data_mat = loadExprMat(
                        "{}/{}-{}_{}pi-data_mat.csv".format(sim_dir_path, num_genes, simulation_type, t)).values
                    log_sim_data_mat = np.log10(sim_data_mat + 1)
                stats_dict["{}_{}".format(simulation_type, t)] = computeBasicStats(sim_data_mat)
                hist_dict["{}_{}".format(simulation_type, t)] = computeHist(_minMaxScale(sim_data_mat), bins)
        file_stats_dict[num_genes] = stats_dict
        file_hist_dict[num_genes] = hist_dict
    # Save data
    np.save("{}/old_sim-basic_stats.npy".format(save_dir_path), file_stats_dict)
    np.save("{}/old_sim-hist.npy".format(save_dir_path), file_hist_dict)


def newSimStats():
    save_dir_path = "./res/"
    num_genes_list = ["50hvg", "100hvg"]
    bins = 20
    file_stats_dict = {}
    file_hist_dict = {}
    file_corr_dict = {}
    for num_genes in num_genes_list:
        stats_dict = {}
        hist_dict = {}
        corr_dict = {}
        print("=" * 70)
        print("[ {} ]".format(num_genes))
        # Experimental data
        for data_name in new_simulation_list:
            print("{} experimental data...".format(data_name))
            if num_genes == "50hvg" and data_name == "SERGIO":
                continue
            if data_name == "NORTA":
                t_list = ["pbmc1-Drop", "pbmc2-Drop", "pbmc1-inDrops", "pbmc2-inDrops"] + \
                         ["Cortex1-10xChromium", "Cortex2-10xChromium", "Cortex1-Smart_seq2", "Cortex2-Smart_seq2"]
                for t in t_list:
                    exp_data_mat = loadExprMat("../data/simulated/new/{}-{}-NORTA-data_mat.csv".format(t, num_genes)).values
                    stats_dict["NORTA-{}-{}".format(t, num_genes)] = computeBasicStats(exp_data_mat)
                    hist_dict["NORTA-{}-{}".format(t, num_genes)] = computeHist(_minMaxScale(exp_data_mat), bins)
                    corr_dict["NORTA-{}-{}".format(t, num_genes)] = computeCorr(exp_data_mat)
            else:
                t_list = [1, 5, 10, 15, 20]
                for t in t_list:
                    exp_data_mat = loadExprMat("../data/SERGIO_simulation_all/100gene-9groups-{}sparsity.csv".format(t)).values
                    stats_dict["SERGIO-{}-{}".format(t, num_genes)] = computeBasicStats(exp_data_mat)
                    hist_dict["SERGIO-{}-{}".format(t, num_genes)] = computeHist(_minMaxScale(exp_data_mat), bins)
                    corr_dict["SERGIO-{}-{}".format(t, num_genes)] = computeCorr(exp_data_mat)
        file_stats_dict[num_genes] = stats_dict
        file_hist_dict[num_genes] = hist_dict
        file_corr_dict[num_genes] = corr_dict
    # Save data
    np.save("{}/new_sim-basic_stats.npy".format(save_dir_path), file_stats_dict)
    np.save("{}/new_sim-hist.npy".format(save_dir_path), file_hist_dict)
    np.save("{}/new_sim-corr.npy".format(save_dir_path), file_corr_dict)

# ======================================

all_models_list = ["pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet", "ZILGM"]

def _getF1(metric_df):
    res_df = metric_df.groupby("Model")["F1"]
    return res_df.max()  # .reset_index()


def _getFDR(metric_df):
    metric_df["FDR"] =  metric_df.FP / ( metric_df.FP +  metric_df.TP)
    res_df = (
        metric_df[metric_df.groupby("Model")["F1"].transform(max) == metric_df["F1"]]
            .drop_duplicates(subset="Model",keep="first")
            [["Model", "Sensitivity", "FDR"]]
    )
    fdr_df = res_df.groupby("Model")["FDR"].mean()
    return fdr_df


def _auc(df, xaxis, yaxis):
    df = df.reset_index(drop=True)
    if xaxis == "FPR":
        df.loc[len(df.index)] = [0.0, 0.0]
        df.loc[len(df.index)] = [1.0, 1.0]
    else:
        df.loc[len(df.index)] = [0.0, 1.0]
        df.loc[len(df.index)] = [1.0, 0.0]
    df = df.sort_values(by=xaxis, ascending=True)
    auc_val = auc(df[xaxis], df[yaxis])
    return auc_val


def _getAUROC(metric_df):
    metric_df = metric_df[~pd.isna(metric_df.F1)]
    metric_df["TPR"] =  metric_df.TP / ( metric_df.FN +  metric_df.TP)
    metric_df["FPR"] =  metric_df.FP / ( metric_df.FP +  metric_df.TN)
    auroc = metric_df.groupby("Model")["TPR", "FPR"].apply(lambda x: _auc(x, "FPR", "TPR"))
    return auroc


def _getPRC(metric_df):
    metric_df = metric_df[~pd.isna(metric_df.F1)]
    metric_df["precision"] = metric_df.TP / (metric_df.FP + metric_df.TP)
    metric_df["recall"] = metric_df.TP / (metric_df.TP + metric_df.FN)
    auprc = metric_df.groupby("Model")["precision", "recall"].apply(lambda x: _auc(x, "recall", "precision"))
    return auprc


def _loadZILGMRes(dir_path, filename):
    if filename in os.listdir(dir_path):
        return loadMetric("{}/{}".format(dir_path, filename))
    else:
        return pd.DataFrame()


def _loadPIDCMRes(dir_path, filename):
    if filename in os.listdir(dir_path):
        return loadMetric("{}/{}".format(dir_path, filename))
    else:
        return pd.DataFrame()


def methodPerformance():
    save_dir_path = "./res/"
    norta_res_dir_path = "../res/logs/evaluation/simulation"
    zilgm_res_dir_path = "../res/logs/evaluation/ZILGM_new_sim"
    norta_metric = {}
    for num_genes in ["50hvg", "100hvg"]:
        norta_metric_file_list = [each for each in os.listdir(norta_res_dir_path)
                                  if ("NORTA" in each) and ("-process" not in each) and (num_genes in each) and ("PIDC" not in each)]
        norta_metric_dict = {"-".join(each.split("-")[:3]):
                                 pd.concat([
                                     loadMetric("{}/{}".format(norta_res_dir_path, each)),
                                     _loadZILGMRes(zilgm_res_dir_path, each),
                                 ], axis=0)
                             for each in norta_metric_file_list}
        norta_raw_F1 = {each: _getF1(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_FDR = {each: _getFDR(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_AUROC = {each: _getAUROC(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_AUPRC = {each: _getPRC(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_metric[num_genes] = {"F1": norta_raw_F1, "FDR": norta_raw_FDR, "AUROC": norta_raw_AUROC, "AUPRC": norta_raw_AUPRC}
    # -----
    sergio_metric = {}
    sergio_res_dir_path = "../res/logs/evaluation/SERGIO_simulation_all"
    zilgm_res_dir_path = "../res/logs/evaluation/ZILGM_model_SERGIO"
    sergio_metric_file_list = [each for each in os.listdir(sergio_res_dir_path) if ("-process" not in each) and ("PIDC" not in each)]
    sergio_metric_dict = {"-".join(each.split("-")[:3]):
                                 pd.concat([
                                     loadMetric("{}/{}".format(sergio_res_dir_path, each)),
                                     _loadZILGMRes(zilgm_res_dir_path, "-".join(each.split("-")[:3])+"sparsity-metrics.csv"),
                                     loadMetric("{}/{}".format(sergio_res_dir_path, "-".join(each.split("-")[:3] + ["PIDC", "metrics.csv"])))
                                 ], axis=0)
                             for each in sergio_metric_file_list}
    sergio_raw_F1 = {each: _getF1(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_FDR = {each: _getFDR(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_AUROC = {each: _getAUROC(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_AUPRC = {each: _getPRC(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_metric["100hvg"] = {"F1": sergio_raw_F1, "FDR": sergio_raw_FDR, "AUROC": sergio_raw_AUROC, "AUPRC": sergio_raw_AUPRC}
    # -----
    np.save("{}/NORTA_sim_metric.npy".format(save_dir_path), norta_metric)
    np.save("{}/SERGIO_sim_metric.npy".format(save_dir_path), sergio_metric)


def realPerformance():
    save_dir_path = "./res/"
    res_dir_path = "../res/logs/evaluation/real"
    metric_file_list = [each for each in os.listdir(res_dir_path) if ("PIDC" not in each)]
    # metric_dict = {each: loadMetric("{}/{}".format(res_dir_path, each)) for each in metric_file_list}
    metric_dict = {each: pd.concat([
        loadMetric("{}/{}".format(res_dir_path, each)),
        loadMetric("{}/{}".format(res_dir_path, "-".join(each.split("-")[:4]) + "-PIDC-metrics.csv")),
        ], axis=0)
        for each in metric_file_list}
    raw_F1 = {each: _getF1(metric_dict[each]) for each in metric_dict}
    raw_FDR = {each: _getFDR(metric_dict[each]) for each in metric_dict}
    raw_AUROC = {each: _getAUROC(metric_dict[each]) for each in metric_dict}
    raw_AUPRC = {each: _getPRC(metric_dict[each]) for each in metric_dict}
    metric = {"F1": raw_F1, "FDR": raw_FDR, "AUROC": raw_AUROC, "AUPRC": raw_AUPRC}
    # metric = {"F1": raw_F1, "FDR": raw_FDR}
    # -----
    np.save("{}/experimental_metric.npy".format(save_dir_path), metric)


def _getTPR(metric_df):
    res_df = (
        metric_df[metric_df.groupby("Model")["F1"].transform(max) == metric_df["F1"]]
            .drop_duplicates(subset="Model", keep="first")
        [["Model", "TPR", "Precision", "F1"]]
    )
    tpr_df = res_df.groupby("Model")["TPR"].mean()
    pre_df = res_df.groupby("Model")["Precision"].mean()
    return tpr_df


def _getPrecision(metric_df):
    res_df = (
        metric_df[metric_df.groupby("Model")["F1"].transform(max) == metric_df["F1"]]
            .drop_duplicates(subset="Model", keep="first")
        [["Model", "TPR", "Precision", "F1"]]
    )
    pre_df = res_df.groupby("Model")["Precision"].mean()
    return pre_df

def _addPrecisionAndTPR(metric_df):
    metric_df["TPR"] = metric_df.TP / (metric_df.TP + metric_df.FN)
    metric_df["Precision"] = metric_df.TP / (metric_df.TP + metric_df.FP)
    return metric_df


def realPerformanceExtra():
    save_dir_path = "./res/"
    res_dir_path = "../res/logs/evaluation/real"
    metric_file_list = [each for each in os.listdir(res_dir_path) if ("PIDC" not in each)]
    # metric_dict = {each: loadMetric("{}/{}".format(res_dir_path, each)) for each in metric_file_list}
    metric_dict = {each: pd.concat([
        loadMetric("{}/{}".format(res_dir_path, each)),
        _loadPIDCMRes(
            res_dir_path,
            "-".join(each.split("-")[:4]) + "-PIDC-metrics.csv" if "pbmc" in each or "Cortex" in each
            else "-".join(each.split("-")[:2]) + "-PIDC-metrics.csv"),
        ], axis=0)
        for each in metric_file_list}
    metric_dict = {each: _addPrecisionAndTPR(metric_dict[each]) for each in metric_dict}
    raw_TPR = {each: _getTPR(metric_dict[each]) for each in metric_dict}
    raw_precision = {each: _getPrecision(metric_dict[each]) for each in metric_dict}
    metric = {"TPR": raw_TPR, "Precision": raw_precision}
    # -----
    np.save("{}/experimental_other_metric.npy".format(save_dir_path), metric)


def highDimPerformance():
    save_dir_path = "./res/"
    res_dir_path = "../res/logs/evaluation/high_dim"
    zilgm_dir_path = "../res/logs/evaluation/ZILGM_model_highdim"
    metric_file_list = [each for each in os.listdir(res_dir_path) if ("PIDC" not in each)]
    metric_dict = {each: pd.concat([
                         loadMetric("{}/{}".format(res_dir_path, each)),
                         loadMetric("{}/{}".format(res_dir_path, "-".join(each.split("-")[:4])+"-PIDC-metrics.csv")),
                         _loadZILGMRes(zilgm_dir_path, each)
                         ], axis=0)
        for each in metric_file_list}
    raw_F1 = {each: _getF1(metric_dict[each]) for each in metric_dict}
    raw_FDR = {each: _getFDR(metric_dict[each]) for each in metric_dict}
    raw_AUROC = {each: _getAUROC(metric_dict[each]) for each in metric_dict}
    raw_AUPRC = {each: _getPRC(metric_dict[each]) for each in metric_dict}
    metric = {"F1": raw_F1, "FDR": raw_FDR, "AUROC": raw_AUROC, "AUPRC": raw_AUPRC}
    # metric = {"F1": raw_F1, "FDR": raw_FDR}
    # -----
    np.save("{}/high_dim_metric.npy".format(save_dir_path), metric)


def SERGIO400GenePerformance():
    save_dir_path = "./res/"
    res_dir_path = "../res/logs/evaluation/SERGIO_400genes"
    metric_file_list = [each for each in os.listdir(res_dir_path) if ("PIDC" not in each) and ("ZILGM" not in each)]
    metric_dict = {each: pd.concat([
                         loadMetric("{}/{}".format(res_dir_path, each)),
                         _loadZILGMRes(res_dir_path, "{}".format("-".join(each.split("-")[:4])+"-PIDC-metrics.csv")),
                         _loadZILGMRes(res_dir_path, "{}".format("-".join(each.split("-")[:4])+"-ZILGM-metrics.csv"))
                    ], axis=0)
        for each in metric_file_list}
    raw_F1 = {each: _getF1(metric_dict[each]) for each in metric_dict}
    raw_FDR = {each: _getFDR(metric_dict[each]) for each in metric_dict}
    raw_AUROC = {each: _getAUROC(metric_dict[each]) for each in metric_dict}
    raw_AUPRC = {each: _getPRC(metric_dict[each]) for each in metric_dict}
    metric = {"F1": raw_F1, "FDR": raw_FDR, "AUROC": raw_AUROC, "AUPRC": raw_AUPRC}
    # -----
    np.save("{}/SERGIO_400genes_metric.npy".format(save_dir_path), metric)


def oldSimPerformance():
    save_dir_path = "./res/"
    res_dir_path = "../res/logs/evaluation/scLink_simulation"
    scLink_metric = {}
    for num_genes in ["50hvg", "100hvg"]:
        metric_file_list = [each for each in os.listdir(res_dir_path) if (num_genes in each) and ("PIDC" not in each)]
        metric_dict = {each: pd.concat([
            loadMetric("{}/{}".format(res_dir_path, each)),
            loadMetric("{}/{}".format(res_dir_path, "-".join(each.split("-")[:2] + ["PIDC", "metrics.csv"]))),
            loadMetric("{}/{}".format(res_dir_path, "-".join(each.split("-")[:2] + ["ZILGM", "metrics.csv"])))
            ], axis = 0)
            for each in metric_file_list}
        raw_F1 = {each: _getF1(metric_dict[each]) for each in metric_dict}
        raw_FDR = {each: _getFDR(metric_dict[each]) for each in metric_dict}
        raw_AUROC = {each: _getAUROC(metric_dict[each]) for each in metric_dict}
        raw_AUPRC = {each: _getPRC(metric_dict[each]) for each in metric_dict}
        scLink_metric[num_genes] = {"F1": raw_F1, "FDR": raw_FDR, "AUROC": raw_AUROC, "AUPRC": raw_AUPRC}
    # metric = {"F1": raw_F1, "FDR": raw_FDR}
    # -----
    res_dir_path = "../res/logs/evaluation/ZILGM_simulation"
    zilgm_metric = {}
    for num_genes in ["50hvg", "100hvg"]:
        metric_file_list = [each for each in os.listdir(res_dir_path) if (num_genes in each) and ("PIDC" not in each)]
        # metric_dict = {each: loadMetric("{}/{}".format(res_dir_path, each)) for each in metric_file_list}
        metric_dict = {each: pd.concat([
            loadMetric("{}/{}".format(res_dir_path, each)),
            loadMetric("{}/{}".format(res_dir_path, "-".join(each.split("-")[:2] + ["PIDC", "metrics.csv"]))),
            loadMetric("{}/{}".format(res_dir_path, "-".join(each.split("-")[:2] + ["ZILGM", "metrics.csv"])))
            ], axis=0)
            for each in metric_file_list}
        raw_F1 = {each: _getF1(metric_dict[each]) for each in metric_dict}
        raw_FDR = {each: _getFDR(metric_dict[each]) for each in metric_dict}
        raw_AUROC = {each: _getAUROC(metric_dict[each]) for each in metric_dict}
        raw_AUPRC = {each: _getPRC(metric_dict[each]) for each in metric_dict}
        zilgm_metric[num_genes] = {"F1": raw_F1, "FDR": raw_FDR, "AUROC": raw_AUROC, "AUPRC": raw_AUPRC}
    # -----
    np.save("{}/ZIGaussian_sim_metric.npy".format(save_dir_path), scLink_metric)
    np.save("{}/ZIPoisson_sim_metric.npy".format(save_dir_path), zilgm_metric)


def scLinkExperimentalPerformance():
    save_dir_path = "./res/"
    res_dir_path = "../res/logs/evaluation/real"
    cell_name = ["T_cell", "skeletal_muscle_satellite_stem_cell", "type_B_pancreatic_cell"]
    net_name = ["STRING", "TRRUST"]
    metric_file_list = ["{}-{}-metrics.csv".format(each[0], each[1]) for each in itertools.product(cell_name, net_name)]
    # metric_file_list = [each for each in os.listdir(res_dir_path) if ("PIDC" not in each)]
    # metric_dict = {each: loadMetric("{}/{}".format(res_dir_path, each)) for each in metric_file_list}
    metric_dict = {each: loadMetric("{}/{}".format(res_dir_path, each)) for each in metric_file_list}
    raw_F1 = {each: _getF1(metric_dict[each]) for each in metric_dict}
    raw_FDR = {each: _getFDR(metric_dict[each]) for each in metric_dict}
    raw_AUROC = {each: _getAUROC(metric_dict[each]) for each in metric_dict}
    raw_AUPRC = {each: _getPRC(metric_dict[each]) for each in metric_dict}
    metric = {"F1": raw_F1, "FDR": raw_FDR, "AUROC": raw_AUROC, "AUPRC": raw_AUPRC}
    # metric = {"F1": raw_F1, "FDR": raw_FDR}
    # -----
    np.save("{}/scLink_paper_experimental_metric.npy".format(save_dir_path), metric)

# ======================================

def preprocessingPerformance():
    save_dir_path = "./res/"
    norta_res_dir_path = "../res/logs/evaluation/simulation"
    zilgm_res_dir_path = "../res/logs/evaluation/ZILGM_new_sim"
    norta_metric = {}
    for num_genes in ["50hvg", "100hvg"]:
        # norta_metric_file_list = [each for each in os.listdir(norta_res_dir_path)
        #                           if ("NORTA" in each) and ("-process" in each) and (num_genes) in each]
        # norta_metric_dict = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(norta_res_dir_path, each)) for each in norta_metric_file_list}
        norta_metric_file_list = [each for each in os.listdir(norta_res_dir_path)
                                  if ("NORTA" in each) and ("-process" in each) and (num_genes in each) and (
                                              "PIDC" not in each)]
        norta_metric_dict = {"-".join(each.split("-")[:3]):
            pd.concat([
                loadMetric("{}/{}".format(norta_res_dir_path, each)),
                _loadZILGMRes(zilgm_res_dir_path, each),
            ], axis=0)
            for each in norta_metric_file_list}
        norta_raw_F1 = {each: _getF1(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_FDR = {each: _getFDR(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_AUROC = {each: _getAUROC(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_AUPRC = {each: _getPRC(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_metric[num_genes] = {"F1": norta_raw_F1, "FDR": norta_raw_FDR, "AUROC": norta_raw_AUROC, "AUPRC": norta_raw_AUPRC}
    # -----
    sergio_metric = {}
    sergio_res_dir_path = "../res/logs/evaluation/SERGIO_simulation_all"
    zilgm_res_dir_path = "../res/logs/evaluation/ZILGM_model_SERGIO"
    # sergio_metric_file_list = [each for each in os.listdir(sergio_res_dir_path) if ("-process" in each)]
    # sergio_metric_dict = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(sergio_res_dir_path, each)) for each in
    #                      sergio_metric_file_list}
    sergio_metric_file_list = [each for each in os.listdir(sergio_res_dir_path) if
                               ("-process" in each) and ("PIDC" not in each)]
    sergio_metric_dict = {"-".join(each.split("-")[:3]):
        pd.concat([
            loadMetric("{}/{}".format(sergio_res_dir_path, each)),
            _loadZILGMRes(zilgm_res_dir_path, "-".join(each.split("-")[:3]) + "sparsity-metrics.csv"),
            loadMetric("{}/{}".format(sergio_res_dir_path, "-".join(each.split("-")[:3] + ["PIDC", "metrics.csv"])))
        ], axis=0)
        for each in sergio_metric_file_list}
    sergio_raw_F1 = {each: _getF1(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_FDR = {each: _getFDR(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_AUROC = {each: _getAUROC(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_AUPRC = {each: _getPRC(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_metric["100hvg"] = {"F1": sergio_raw_F1, "FDR": sergio_raw_FDR, "AUROC": sergio_raw_AUROC, "AUPRC": sergio_raw_AUPRC}
    # -----
    np.save("{}/NORTA_preprocess_metric.npy".format(save_dir_path), norta_metric)
    np.save("{}/SERGIO_preprocess_metric.npy".format(save_dir_path), sergio_metric)


def loadNormRes(dir_path, filename):
    if filename in os.listdir(dir_path):
        return loadMetric("{}/{}".format(dir_path, filename))
    else:
        return pd.DataFrame()


def sctransformPerformance():
    save_dir_path = "./res/"
    sctransform_res_dir_path = "../res/logs/evaluation/other_norm"
    sctransform_metric = {}
    for num_genes in ["100hvg"]:
        sctransform_metric_file_list = [each for each in os.listdir(sctransform_res_dir_path)
                                  if ("sctransform" in each) and ("100gene" not in each)]
        # -----
        sctransform_dict = {"-".join(each.split("-")[:3]): loadNormRes(sctransform_res_dir_path,
                                                                       "{}-sctransform-metrics.csv".format(
                                                                           "-".join(each.split("-")[:3]))) for each in
                            sctransform_metric_file_list}
        sctransform_F1 = {each: _getF1(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_FDR = {each: _getFDR(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_AUROC = {each: _getAUROC(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_AUPRC = {each: _getPRC(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_metric[num_genes] = {"F1": sctransform_F1, "FDR": sctransform_FDR, "AUROC": sctransform_AUROC,
                                   "AUPRC": sctransform_AUPRC}
    np.save("{}/NORTA_sctransform_metric.npy".format(save_dir_path), sctransform_metric)
    # -----
    sctransform_res_dir_path = "../res/logs/evaluation/other_norm"
    sctransform_metric = {}
    for num_genes in ["100hvg"]:
        sctransform_metric_file_list = [each for each in os.listdir(sctransform_res_dir_path)
                                        if ("sctransform" in each) and ("100gene" in each)]
        # -----
        sctransform_dict = {"-".join(each.split("-")[:3]): loadNormRes(sctransform_res_dir_path,
                                                                       "{}-sctransform-metrics.csv".format(
                                                                           "-".join(each.split("-")[:3]))) for each in
                            sctransform_metric_file_list}
        sctransform_F1 = {each: _getF1(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_FDR = {each: _getFDR(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_AUROC = {each: _getAUROC(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_AUPRC = {each: _getPRC(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_metric[num_genes] = {"F1": sctransform_F1, "FDR": sctransform_FDR, "AUROC": sctransform_AUROC,
                                         "AUPRC": sctransform_AUPRC}
    np.save("{}/SERGIO_sctransform_metric.npy".format(save_dir_path), sctransform_metric)


def psiNormPerformance():
    save_dir_path = "./res/"
    sctransform_res_dir_path = "../res/logs/evaluation/other_norm"
    sctransform_metric = {}
    for num_genes in ["100hvg"]:
        sctransform_metric_file_list = [each for each in os.listdir(sctransform_res_dir_path)
                                  if ("psiNorm" in each) and ("100gene" not in each) and ("ZILGM" not in each)]
        # -----
        sctransform_dict = {"-".join(each.split("-")[:3]): pd.concat([
                loadNormRes(sctransform_res_dir_path, "{}-psiNorm-metrics.csv".format("-".join(each.split("-")[:3]))),
                _loadZILGMRes(sctransform_res_dir_path, "{}-psiNorm-ZILGM-metrics.csv".format("-".join(each.split("-")[:3]))),
            ], axis=0)
            for each in sctransform_metric_file_list}
        sctransform_F1 = {each: _getF1(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_FDR = {each: _getFDR(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_AUROC = {each: _getAUROC(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_AUPRC = {each: _getPRC(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_metric[num_genes] = {"F1": sctransform_F1, "FDR": sctransform_FDR, "AUROC": sctransform_AUROC,
                                   "AUPRC": sctransform_AUPRC}
    np.save("{}/NORTA_psiNorm_metric.npy".format(save_dir_path), sctransform_metric)
    # -----
    sctransform_res_dir_path = "../res/logs/evaluation/other_norm"
    sctransform_metric = {}
    for num_genes in ["100hvg"]:
        sctransform_metric_file_list = [each for each in os.listdir(sctransform_res_dir_path)
                                        if ("psiNorm" in each) and ("100gene" in each) and ("ZILGM" not in each)]
        # -----
        sctransform_dict = {"-".join(each.split("-")[:3]): pd.concat([
            loadNormRes(sctransform_res_dir_path, "{}-psiNorm-metrics.csv".format("-".join(each.split("-")[:3]))),
            _loadZILGMRes(sctransform_res_dir_path,
                          "{}-psiNorm-ZILGM-metrics.csv".format("-".join(each.split("-")[:3]))),
        ], axis=0)
            for each in sctransform_metric_file_list}
        sctransform_F1 = {each: _getF1(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_FDR = {each: _getFDR(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_AUROC = {each: _getAUROC(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_AUPRC = {each: _getPRC(sctransform_dict[each]) for each in sctransform_dict}
        sctransform_metric[num_genes] = {"F1": sctransform_F1, "FDR": sctransform_FDR, "AUROC": sctransform_AUROC,
                                         "AUPRC": sctransform_AUPRC}
    np.save("{}/SERGIO_psiNorm_metric.npy".format(save_dir_path), sctransform_metric)


def pseudobulkPerformance():
    save_dir_path = "./res/"
    norta_res_dir_path = "../res/logs/evaluation/pseudo_bulk"
    zilgm_dir_path = "../res/logs/evaluation/ZILGM_model_pseudobulk"
    norta_metric = {}
    for num_genes in ["50hvg", "100hvg"]:
        # norta_metric_file_list = [each for each in os.listdir(norta_res_dir_path)
        #                           if ("NORTA" in each) and ("-process" not in each) and (num_genes) in each]
        # norta_metric_dict = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(norta_res_dir_path, each)) for each in norta_metric_file_list}
        norta_metric_file_list = [each for each in os.listdir(norta_res_dir_path)
                                  if ("NORTA" in each) and ("-process" not in each) and (num_genes in each) and (
                                              "PIDC" not in each)]
        norta_metric_dict = {"-".join(each.split("-")[:3]):
            pd.concat([
                loadMetric("{}/{}".format(norta_res_dir_path, each)),
                loadMetric("{}/{}".format(norta_res_dir_path, "-".join(each.split("-")[:4])+"-PIDC-metrics.csv")),
                _loadZILGMRes(zilgm_dir_path, "-".join(each.split("-")[:-1])+"metrics.csv"),
            ], axis=0)
            for each in norta_metric_file_list}
        norta_raw_F1 = {each: _getF1(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_FDR = {each: _getFDR(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_AUROC = {each: _getAUROC(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_AUPRC = {each: _getPRC(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_metric[num_genes] = {"F1": norta_raw_F1, "FDR": norta_raw_FDR, "AUROC": norta_raw_AUROC, "AUPRC": norta_raw_AUPRC}
    # -----
    sergio_metric = {}
    sergio_res_dir_path = "../res/logs/evaluation/SERGIO_pseudobulk"
    zilgm_dir_path = "../res/logs/evaluation/ZILGM_model_SERGIO_pseudobulk"
    # sergio_metric_file_list = [each for each in os.listdir(sergio_res_dir_path) if ("-process" not in each)]
    # sergio_metric_dict = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(sergio_res_dir_path, each)) for each in
    #                      sergio_metric_file_list}
    sergio_metric_file_list = [each for each in os.listdir(sergio_res_dir_path) if
                               ("-process" not in each) and ("PIDC" not in each)]
    sergio_metric_dict = {"-".join(each.split("-")[:3]):
        pd.concat([
            loadMetric("{}/{}".format(sergio_res_dir_path, each)),
            _loadZILGMRes(zilgm_dir_path, "-".join(each.split("-")[:3]) + "-metrics.csv"),
            loadMetric("{}/{}".format(sergio_res_dir_path, "-".join(each.split("-")[:3])+"-PIDC-metrics.csv")) if "20" not in each else pd.DataFrame() # 20 sparsity has only 3 clusters
        ], axis=0)
        for each in sergio_metric_file_list}
    sergio_raw_F1 = {each: _getF1(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_FDR = {each: _getFDR(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_AUROC = {each: _getAUROC(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_AUPRC = {each: _getPRC(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_metric["100hvg"] = {"F1": sergio_raw_F1, "FDR": sergio_raw_FDR, "AUROC": sergio_raw_AUROC, "AUPRC": sergio_raw_AUPRC}
    # -----
    np.save("{}/NORTA_pseudobulk_metric.npy".format(save_dir_path), norta_metric)
    np.save("{}/SERGIO_pseudobulk_metric.npy".format(save_dir_path), sergio_metric)


def imputationPerformance():
    save_dir_path = "./res/"
    norta_res_dir_path = "../res/logs/evaluation/imputed_simulation"
    zilgm_dir_path = "../res/logs/evaluation/ZILGM_model_imputation"
    norta_metric = {}
    for num_genes in ["50hvg", "100hvg"]:
        # norta_metric_file_list = [each for each in os.listdir(norta_res_dir_path)
        #                           if ("NORTA" in each) and ("-process" not in each) and (num_genes) in each]
        # norta_metric_dict = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(norta_res_dir_path, each)) for each in norta_metric_file_list}
        norta_metric_file_list = [each for each in os.listdir(norta_res_dir_path)
                                  if ("NORTA" in each) and ("-process" not in each) and (num_genes in each) and (
                                          "PIDC" not in each)]
        norta_metric_dict = {each:
            pd.concat([
                loadMetric("{}/{}".format(norta_res_dir_path, each)),
                # loadMetric("{}/{}".format(norta_res_dir_path, "-".join(each.split("-")[:5]) + "-PIDC-metrics.csv")),
                _loadPIDCMRes(norta_res_dir_path, "-".join(each.split("-")[:5]) + "-PIDC-metrics.csv"),
                _loadZILGMRes(zilgm_dir_path, each),
            ], axis=0)
            for each in norta_metric_file_list}
        norta_raw_F1 = {each: _getF1(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_FDR = {each: _getFDR(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_AUROC = {each: _getAUROC(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_raw_AUPRC = {each: _getPRC(norta_metric_dict[each]) for each in norta_metric_dict}
        norta_metric[num_genes] = {"F1": norta_raw_F1, "FDR": norta_raw_FDR, "AUROC": norta_raw_AUROC, "AUPRC": norta_raw_AUPRC}
    # -----
    sergio_metric = {}
    sergio_res_dir_path = "../res/logs/evaluation/SERGIO_imputation"
    zilgm_dir_path = "../res/logs/evaluation/ZILGM_model_SERGIO_imputation"
    # sergio_metric_file_list = [each for each in os.listdir(sergio_res_dir_path) if ("-process" not in each)]
    # sergio_metric_dict = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(sergio_res_dir_path, each)) for each in
    #                      sergio_metric_file_list}
    sergio_metric_file_list = [each for each in os.listdir(sergio_res_dir_path) if
                               ("-process" not in each) and ("PIDC" not in each)]
    # sergio_metric_dict = {"-".join(each.split("-")[:3]):
    sergio_metric_dict = {each:
        pd.concat([
            loadMetric("{}/{}".format(sergio_res_dir_path, each)),
            _loadZILGMRes(zilgm_dir_path, "-".join(each.split("-")[:4]) + "-metrics.csv"),
            # loadMetric("{}/{}".format(sergio_res_dir_path, "-".join(each.split("-")[:4]) + "-PIDC-metrics.csv"))
            _loadPIDCMRes(norta_res_dir_path, "-".join(each.split("-")[:4]) + "-PIDC-metrics.csv"),
        ], axis=0)
        for each in sergio_metric_file_list}
    sergio_raw_F1 = {each: _getF1(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_FDR = {each: _getFDR(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_AUROC = {each: _getAUROC(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_raw_AUPRC = {each: _getPRC(sergio_metric_dict[each]) for each in sergio_metric_dict}
    sergio_metric["100hvg"] = {"F1": sergio_raw_F1, "FDR": sergio_raw_FDR, "AUROC": sergio_raw_AUROC, "AUPRC": sergio_raw_AUPRC}
    # -----
    np.save("{}/NORTA_imputation_metric.npy".format(save_dir_path), norta_metric)
    np.save("{}/SERGIO_imputation_metric.npy".format(save_dir_path), sergio_metric)


def diffSettingPerformance():
    save_dir_path = "./res/"
    cell_res_dir_path = "../res/logs/evaluation/diff_cell"
    # norta_metric_dict = {each: loadMetric("{}/{}".format(cell_res_dir_path, each)) for each in
    #                      norta_metric_file_list}
    norta_metric_file_list = [each for each in os.listdir(cell_res_dir_path)
                              if ("-process" not in each) and ("PIDC" not in each) and ("ZILGM" not in each)]
    norta_metric_dict = {each:
        pd.concat([
            loadMetric("{}/{}".format(cell_res_dir_path, each)),
            _loadPIDCMRes(cell_res_dir_path, "-".join(each.split("-")[:4]) + "-PIDC-metrics.csv"),
            _loadZILGMRes(cell_res_dir_path, "-".join(each.split("-")[:4]) + "-ZILGM-metrics.csv")
        ], axis=0)
        for each in norta_metric_file_list}
    norta_raw_F1 = {each: _getF1(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_raw_FDR = {each: _getFDR(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_raw_AUROC = {each: _getAUROC(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_raw_AUPRC = {each: _getPRC(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_cell_metric = {"F1": norta_raw_F1, "FDR": norta_raw_FDR, "AUROC": norta_raw_AUROC, "AUPRC": norta_raw_AUPRC}
    # -----
    graph_res_dir_path = "../res/logs/evaluation/diff_graph"
    # norta_metric_file_list = [each for each in os.listdir(graph_res_dir_path) if ("-process" not in each)]
    # norta_metric_dict = {each: loadMetric("{}/{}".format(graph_res_dir_path, each)) for each in
    #                      norta_metric_file_list}
    norta_metric_file_list = [each for each in os.listdir(graph_res_dir_path)
                              if ("-process" not in each) and ("PIDC" not in each) and ("ZILGM" not in each)]
    norta_metric_dict = {each:
        pd.concat([
            loadMetric("{}/{}".format(graph_res_dir_path, each)),
            _loadPIDCMRes(graph_res_dir_path, "-".join(each.split("-")[:4]) + "-PIDC-metrics.csv"),
            _loadZILGMRes(graph_res_dir_path, "-".join(each.split("-")[:4]) + "-ZILGM-metrics.csv")
        ], axis=0)
        for each in norta_metric_file_list}
    norta_raw_F1 = {each: _getF1(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_raw_FDR = {each: _getFDR(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_raw_AUROC = {each: _getAUROC(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_raw_AUPRC = {each: _getPRC(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_graph_metric = {"F1": norta_raw_F1, "FDR": norta_raw_FDR, "AUROC": norta_raw_AUROC, "AUPRC": norta_raw_AUPRC}
    # -----
    sparsity_res_dir_path = "../res/logs/evaluation/diff_sparsity"
    # norta_metric_file_list = [each for each in os.listdir(sparsity_res_dir_path) if ("-process" not in each)]
    # norta_metric_dict = {each: loadMetric("{}/{}".format(sparsity_res_dir_path, each)) for each in
    #                      norta_metric_file_list}
    norta_metric_file_list = [each for each in os.listdir(sparsity_res_dir_path)
                              if ("-process" not in each) and ("PIDC" not in each) and ("ZILGM" not in each)]
    norta_metric_dict = {each:
        pd.concat([
            loadMetric("{}/{}".format(sparsity_res_dir_path, each)),
            _loadPIDCMRes(sparsity_res_dir_path, "-".join(each.split("-")[:3]) + "-PIDC-metrics.csv"),
            _loadZILGMRes(sparsity_res_dir_path, "-".join(each.split("-")[:3]) + "-ZILGM-metrics.csv")
        ], axis=0)
        for each in norta_metric_file_list}
    norta_raw_F1 = {each: _getF1(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_raw_FDR = {each: _getFDR(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_raw_AUROC = {each: _getAUROC(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_raw_AUPRC = {each: _getPRC(norta_metric_dict[each]) for each in norta_metric_dict}
    norta_sparsity_metric = {"F1": norta_raw_F1, "FDR": norta_raw_FDR, "AUROC": norta_raw_AUROC, "AUPRC": norta_raw_AUPRC}
    # -----
    np.save("{}/diff_cell_metric.npy".format(save_dir_path), norta_cell_metric)
    np.save("{}/diff_graph_metric.npy".format(save_dir_path), norta_graph_metric)
    np.save("{}/diff_sparsity_metric.npy".format(save_dir_path), norta_sparsity_metric)


def qcPerformance():
    save_dir_path = "./res/"
    res_dir_path = "../res/logs/evaluation/quality_control"
    zilgm_res_dir_path = "../res/logs/evaluation/ZILGM_model_qc"
    qc_metric = {}
    for num_genes in ["50hvg", "100hvg"]:
        metric_file_list = [each for each in os.listdir(res_dir_path) if (num_genes in each)]
        # metric_dict = {each: loadMetric("{}/{}".format(res_dir_path, each)) for each in metric_file_list}
        metric_dict = {each:
            pd.concat([
                loadMetric("{}/{}".format(res_dir_path, each)),
                _loadZILGMRes(zilgm_res_dir_path, each)
            ], axis=0)
            for each in metric_file_list}
        raw_F1 = {each: _getF1(metric_dict[each]) for each in metric_dict}
        raw_FDR = {each: _getFDR(metric_dict[each]) for each in metric_dict}
        raw_AUROC = {each: _getAUROC(metric_dict[each]) for each in metric_dict}
        raw_AUPRC = {each: _getPRC(metric_dict[each]) for each in metric_dict}
        qc_metric[num_genes] = {"F1": raw_F1, "FDR": raw_FDR, "AUROC": raw_AUROC, "AUPRC": raw_AUPRC}
    # -----
    np.save("{}/diff_qc_metric.npy".format(save_dir_path), qc_metric)


# ======================================

def _takeUpper(mat):
    mat = np.asarray(mat)
    return mat[np.triu_indices(mat.shape[0], k=1)]


def diffSparsityCMD():
    save_dir_path = "./res/"
    data_dir_path = "../data/diff_sparsity/"
    sparsity_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, "full"]
    diff_sparsity_data_list = {
        i: loadExprMat("{}/100gene-hub-{}-NORTA-data_mat.csv".format(data_dir_path, "{}sparsity".format(i) if i != "full" else i))
        for i in sparsity_list
    }
    diff_sparsity_corr = {i: computeCorr(diff_sparsity_data_list[i]) for i in diff_sparsity_data_list}
    net_mat = loadNetMat("../data/diff_sparsity/100gene-hub-net_mat.csv").values
    CMD_list = {each: _CMD(diff_sparsity_corr[each]["PCC"], net_mat) for each in diff_sparsity_corr}
    # -----
    triu_net_mat = _takeUpper(net_mat)
    triu_corr_list = [_takeUpper(diff_sparsity_corr[each]["PCC"]) for each in diff_sparsity_corr]
    pos_idx = np.where(triu_net_mat != 0)[0]
    neg_idx = np.where(triu_net_mat == 0)[0]
    element_abs_diff = [np.abs(triu_net_mat - each) for each in triu_corr_list]
    pos_abs_diff = [each[pos_idx] for each in element_abs_diff]
    neg_abs_diff = [each[neg_idx] for each in element_abs_diff]
    metric = {
        "CMD": CMD_list,
        "pos_abs_diff": pos_abs_diff,
        "neg_abs_diff": neg_abs_diff
    }
    np.save("{}/diff_sparsity_corr.npy".format(save_dir_path), metric)

# ======================================

def _preprocess(data):
    data = data / np.sum(data, axis=1)[:, np.newaxis] * 1e4
    data[np.isnan(data)] = 0.0
    data = np.log10(data + 1)
    return data


def preprocessStats():
    save_dir_path = "./res/"
    raw_data_dir_path = "../data/simulated/new/"
    data_stats = {}
    bins = 20
    for num_genes in ["50hvg", "100hvg"]:
        stats_dict = {}
        hist_dict = {}
        corr_dict = {}
        # -----
        raw_file_list = [each for each in os.listdir(raw_data_dir_path) if ("NORTA" in each) and (num_genes) in each]
        raw_data_dict = {"-".join(each.split("-")[:3]): loadExprMat("{}/{}".format(raw_data_dir_path, each)).values
                         for each in raw_file_list}
        processed_data_dict = {each: _preprocess(raw_data_dict[each]) for each in raw_data_dict}
        # -----
        stats_dict["raw"] = {each: computeBasicStats(raw_data_dict[each]) for each in raw_data_dict}
        stats_dict["processed"] = {each: computeBasicStats(processed_data_dict[each]) for each in processed_data_dict}
        hist_dict["raw"] = {each: computeHist(_minMaxScale(raw_data_dict[each]), bins) for each in raw_data_dict}
        hist_dict["processed"] = {each: computeHist(_minMaxScale(processed_data_dict[each]), bins) for each in processed_data_dict}
        corr_dict["raw"] = {each: computeCorr(raw_data_dict[each]) for each in raw_data_dict}
        corr_dict["processed"] = {each: computeCorr(processed_data_dict[each]) for each in processed_data_dict}
        # -----
        data_stats[num_genes] = {"stats": stats_dict, "hist": hist_dict, "corr": corr_dict}
    np.save("{}/preprocess_data_metric.npy".format(save_dir_path), data_stats)


def pseudobulkStats():
    save_dir_path = "./res/"
    raw_data_dir_path = "../data/simulated/new/"
    pseudobulk_data_dir_path = "../data/pseudo_bulk/"
    data_stats = {}
    bins = 20
    for num_genes in ["50hvg", "100hvg"]:
        stats_dict = {}
        hist_dict = {}
        corr_dict = {}
        # -----
        raw_file_list = [each for each in os.listdir(raw_data_dir_path) if ("NORTA" in each) and (num_genes in each)]
        pseudo_file_list = [each for each in os.listdir(pseudobulk_data_dir_path) if ("NORTA" in each) and (num_genes in each)]
        raw_data_dict = {"-".join(each.split("-")[:3]): loadExprMat("{}/{}".format(raw_data_dir_path, each)).values
                         for each in raw_file_list}
        processed_data_dict = {"-".join(each.split("-")[:3]): loadExprMat("{}/{}".format(pseudobulk_data_dir_path, each)).values
                         for each in pseudo_file_list}
        # -----
        stats_dict["raw"] = {each: computeBasicStats(raw_data_dict[each]) for each in raw_data_dict}
        stats_dict["pseudobulk"] = {each: computeBasicStats(processed_data_dict[each]) for each in processed_data_dict}
        hist_dict["raw"] = {each: computeHist(_minMaxScale(raw_data_dict[each]), bins) for each in raw_data_dict}
        hist_dict["pseudobulk"] = {each: computeHist(_minMaxScale(processed_data_dict[each]), bins) for each in processed_data_dict}
        corr_dict["raw"] = {each: computeCorr(raw_data_dict[each]) for each in raw_data_dict}
        corr_dict["pseudobulk"] = {each: computeCorr(processed_data_dict[each]) for each in processed_data_dict}
        # -----
        data_stats[num_genes] = {"stats": stats_dict, "hist": hist_dict, "corr": corr_dict}
    np.save("{}/pseudobulk_data_metric.npy".format(save_dir_path), data_stats)


def imputationStats():
    save_dir_path = "./res/"
    raw_data_dir_path = "../data/simulated/new/"
    imputation_data_dir_path = "../data/imputed_simulation/"
    data_stats = {}
    bins = 20
    for num_genes in ["50hvg", "100hvg"]:
        stats_dict = {}
        hist_dict = {}
        corr_dict = {}
        # -----
        raw_file_list = [each for each in os.listdir(raw_data_dir_path) if ("NORTA" in each) and (num_genes in each)]
        SAVER_file_list = [each for each in os.listdir(imputation_data_dir_path) if ("NORTA" in each) and (num_genes in each) and ("SAVER" in each)]
        MAGIC_file_list = [each for each in os.listdir(imputation_data_dir_path) if ("NORTA" in each) and (num_genes in each) and ("MAGIC" in each)]
        raw_data_dict = {"-".join(each.split("-")[:3]): loadExprMat("{}/{}".format(raw_data_dir_path, each)).values
                         for each in raw_file_list}
        SAVER_data_dict = {"-".join(each.split("-")[:3]): loadExprMat("{}/{}".format(imputation_data_dir_path, each)).values
                         for each in SAVER_file_list}
        MAGIC_data_dict = {
            "-".join(each.split("-")[:3]): loadExprMat("{}/{}".format(imputation_data_dir_path, each)).values
            for each in MAGIC_file_list}
        # -----
        stats_dict["raw"] = {each: computeBasicStats(raw_data_dict[each]) for each in raw_data_dict}
        stats_dict["SAVER"] = {each: computeBasicStats(SAVER_data_dict[each]) for each in SAVER_data_dict}
        stats_dict["MAGIC"] = {each: computeBasicStats(MAGIC_data_dict[each]) for each in MAGIC_data_dict}
        hist_dict["raw"] = {each: computeHist(_minMaxScale(raw_data_dict[each]), bins) for each in raw_data_dict}
        hist_dict["SAVER"] = {each: computeHist(_minMaxScale(SAVER_data_dict[each]), bins) for each in SAVER_data_dict}
        hist_dict["MAGIC"] = {each: computeHist(_minMaxScale(MAGIC_data_dict[each]), bins) for each in MAGIC_data_dict}
        corr_dict["raw"] = {each: computeCorr(raw_data_dict[each]) for each in raw_data_dict}
        corr_dict["SAVER"] = {each: computeCorr(SAVER_data_dict[each]) for each in SAVER_data_dict}
        corr_dict["MAGIC"] = {each: computeCorr(MAGIC_data_dict[each]) for each in MAGIC_data_dict}
        # -----
        data_stats[num_genes] = {"stats": stats_dict, "hist": hist_dict, "corr": corr_dict}
    np.save("{}/imputation_data_metric.npy".format(save_dir_path), data_stats)


# ======================================


def _normalize(data):
    data = (data / np.sum(data, axis=1)[:, np.newaxis]) * 1e4
    data = np.log10(data + 1.0)
    return data


def compareNORTAProcessingStats():
    num_genes = "100hvg"
    # -----
    file_stats_dict = {}
    print("=" * 70)
    t_list = ["pbmc1-Drop", "pbmc2-Drop", "pbmc1-inDrops", "pbmc2-inDrops"] + \
             ["Cortex1-10xChromium", "Cortex2-10xChromium", "Cortex1-Smart_seq2", "Cortex2-Smart_seq2"]
    for t in t_list:
        raw_data_mat = loadExprMat("../data/simulated/new/{}-{}-NORTA-data_mat.csv".format(t, num_genes)).values
        normalize_data_mat = _normalize(raw_data_mat)
        pseudobulk_data_mat = loadExprMat("../data/pseudo_bulk/{}-{}-NORTA-data_mat.csv".format(t, num_genes)).values
        SAVER_data_mat = loadExprMat("../data/imputed_simulation/{}-{}-NORTA-SAVER-data_mat.csv".format(t, num_genes)).values
        MAGIC_data_mat = loadExprMat("../data/imputed_simulation/{}-{}-NORTA-MAGIC-data_mat.csv".format(t, num_genes)).values
        ALRA_data_mat = loadExprMat("../data/imputed_simulation/{}-{}-NORTA-ALRA-data_mat.csv".format(t, num_genes)).values
        DCA_data_mat = loadExprMat("../data/imputed_simulation/{}-{}-NORTA-DCA-data_mat.csv".format(t, num_genes)).values
        file_stats_dict["{}-{}".format(t, num_genes)] = {
            "raw":computeBasicStats(raw_data_mat),
            "normalize":computeBasicStats(normalize_data_mat),
            "pseudobulk":computeBasicStats(pseudobulk_data_mat),
            "SAVER":computeBasicStats(SAVER_data_mat),
            "MAGIC":computeBasicStats(MAGIC_data_mat),
            "ALRA":computeBasicStats(ALRA_data_mat),
            "DCA":computeBasicStats(DCA_data_mat),
        }
    np.save("./res/processing_NORTA_metric.npy", file_stats_dict)


def compareSERGIOProcessingStats():
    file_stats_dict = {}
    print("=" * 70)
    t_list = ["100gene-9groups-1","100gene-9groups-5",
              "100gene-9groups-10","100gene-9groups-15",
              "100gene-9groups-20" ]
    for t in t_list:
        raw_data_mat = loadExprMat("../data/SERGIO_simulation_all/{}sparsity.csv".format(t)).values
        normalize_data_mat = _normalize(raw_data_mat)
        pseudobulk_data_mat = loadExprMat("../data/SERGIO_pseudobulk/{}.csv".format(t)).values
        MAGIC_data_mat = loadExprMat("../data/SERGIO_imputation/{}-MAGIC-data_mat.csv".format(t)).values
        SAVER_data_mat = loadExprMat("../data/SERGIO_imputation/{}-SAVER-data_mat.csv".format(t)).values
        ALRA_data_mat = loadExprMat("../data/SERGIO_imputation/{}-ALRA-data_mat.csv".format(t)).values
        DCA_data_mat = loadExprMat("../data/SERGIO_imputation/{}-DCA-data_mat.csv".format(t)).values
        file_stats_dict["{}".format(t)] = {
            "raw":computeBasicStats(raw_data_mat),
            "normalize":computeBasicStats(normalize_data_mat),
            "pseudobulk":computeBasicStats(pseudobulk_data_mat),
            "SAVER":computeBasicStats(SAVER_data_mat),
            "MAGIC":computeBasicStats(MAGIC_data_mat),
            "ALRA":computeBasicStats(ALRA_data_mat),
            "DCA":computeBasicStats(DCA_data_mat),
        }
    np.save("./res/processing_SERGIO_metric.npy", file_stats_dict)


if __name__ == '__main__':
    # experimentalStats()
    # oldSimStats()
    # newSimStats()
    # -----
    # oldSimPerformance()
    # methodPerformance()
    # realPerformance()
    # realPerformanceExtra()
    # highDimPerformance()
    # preprocessingPerformance()
    # sctransformPerformance()
    # psiNormPerformance()
    # SERGIO400GenePerformance()
    # pseudobulkPerformance()
    # imputationPerformance()
    # diffSettingPerformance()
    # diffSparsityCMD()
    # qcPerformance()
    # -----
    # preprocessStats()
    # pseudobulkStats()
    # imputationStats()
    # -----
    # compareNORTAProcessingStats()
    # compareSERGIOProcessingStats()
    # -----
    # scLinkExperimentalPerformance()
    pass