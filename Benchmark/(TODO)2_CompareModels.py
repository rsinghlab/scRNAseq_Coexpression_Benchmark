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

def _randomBaseline(filename):
    true_net_path = "../data/simulated/new/"
    filename = "-".join(filename.split("-")[:3])
    true_net_filename = "{}/{}-net_mat.csv".format(true_net_path, filename)
    true_net = loadNetMat(true_net_filename).values
    true_net_triu = (true_net[np.triu_indices(true_net.shape[0])] != 0).astype(int)
    # random_prob = np.count_nonzero(true_net_triu) / len(true_net_triu)
    random_prob = 0.5
    random_net_triu = np.random.choice([0, 1], size=len(true_net_triu), p=[1-random_prob, random_prob])
    random_prob_triu = np.asarray([random_prob for _ in range(len(true_net_triu))])
    random_f1 = f1_score(true_net_triu, random_net_triu)
    random_auroc = roc_auc_score(true_net_triu, random_prob_triu)
    random_pr_curve = precision_recall_curve(true_net_triu, random_prob_triu)
    random_aupr = auc(random_pr_curve[1], random_pr_curve[0]) # x=recall , y=precision
    true_p_idx = np.where(true_net_triu == 1)[0]
    true_n_idx = np.where(true_net_triu == 0)[0]
    random_tpr = np.count_nonzero(true_net_triu[true_p_idx] == random_net_triu[true_p_idx]) / len(true_p_idx)
    random_fpr = np.count_nonzero(true_net_triu[true_n_idx] != random_net_triu[true_n_idx]) / len(true_n_idx)
    return random_f1, random_auroc, random_aupr, random_tpr, random_fpr


def load4DiffSetting(metric_dir_path, data_name, sim_name, num_genes):
    raw_metric_list = [each for each in os.listdir(metric_dir_path)
                   if (data_name in each) and (sim_name in each) and ("process-metric" not in each) and (num_genes in each)]
    processed_metric_list = [each for each in os.listdir(metric_dir_path)
                       if (data_name in each) and (sim_name in each) and ("process-metric" in each) and (num_genes in each)]
    metric_dict = {}
    metric_dict["raw"] = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(metric_dir_path, each)) for each in raw_metric_list}
    metric_dict["processed"] = {"-".join(each.split("-")[:3]): loadMetric("{}/{}".format(metric_dir_path, each)) for each in processed_metric_list}
    metric_dict["baseline"] = {"-".join(each.split("-")[:3]): _randomBaseline(each) for each in processed_metric_list}
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


def getModelTprFpr(metric_df):
    metric_df["FDR"] =  metric_df.FP / ( metric_df.FP +  metric_df.TP)
    res_df = (
        metric_df[metric_df.groupby("Model")["F1"].transform(max) == metric_df["F1"]]
            .drop_duplicates(subset="Model",keep="first")
            [["Model", "Sensitivity", "FDR"]]
    )
    tpr_df = res_df.groupby("Model")["Sensitivity"].mean()
    fdr_df = res_df.groupby("Model")["FDR"].mean()
    return tpr_df, fdr_df


def process4Vis(metric_dict, random_metric_dict):
    model_metrics_list = [[metric_dict[each].loc[m] for each in metric_dict] for m in all_models_list]
    model_metrics_list.append(list(random_metric_dict.values()))
    return model_metrics_list


def process4TprFprVis(metric_dict, random_metric_dict):
    model_tpr_list = [[metric_dict[each][0].loc[m] if m in metric_dict[each][0] else np.nan for each in metric_dict] for m in all_models_list]
    model_fpr_list = [[metric_dict[each][1].loc[m] if m in metric_dict[each][0] else np.nan for each in metric_dict] for m in all_models_list]
    model_tpr_list.append([each[0] for each in list(random_metric_dict.values())])
    model_fpr_list.append([each[1] for each in list(random_metric_dict.values())])
    return model_tpr_list, model_fpr_list

# ======================================

all_models_list = ["pearson", "spearman", "glasso", "scLink", "GENIE3", "PIDC", "scDesign2", "minet"]


def compareModelMetrics(metric_dir_path, data_name):
    for num_genes in ["50hvg", "100hvg"]:
        # for sim_name in ["NORTA", "netSmooth"]:
        for sim_name in ["NORTA"]:
            print("=" * 70)
            print("[ {} | {} ] Loading metrics...".format(sim_name, num_genes))
            metric_dict = load4DiffSetting(metric_dir_path, data_name=data_name, sim_name=sim_name, num_genes=num_genes)
            # Compare F1 score
            raw_best_F1 = {each: getModelPerformance(metric_dict["raw"][each], "F1", func="max") for each in metric_dict["raw"]}
            processed_best_F1 = {each: getModelPerformance(metric_dict["processed"][each], "F1", func="max") for each in metric_dict["processed"]}
            random_F1_dict = {each: metric_dict["baseline"][each][0] for each in metric_dict["baseline"]}
            raw_F1_vis = process4Vis(raw_best_F1, random_F1_dict)
            processed_F1_vis = process4Vis(processed_best_F1, random_F1_dict)
            filename_list = list(raw_best_F1.keys())
            visBox4Metric(raw_F1_vis, all_models_list+["random"], title="{}-{}-raw".format(sim_name, num_genes), metric_name="F1 score")
            visBox4Metric(processed_F1_vis, all_models_list+["random"], title="{}-{}-processed".format(sim_name, num_genes), metric_name="F1 score")
            # Compare TPR and FPR
            raw_best_tpr_fpr = {each: getModelTprFpr(metric_dict["raw"][each]) for each in metric_dict["raw"]}
            processed_best_tpr_fpr = {each: getModelTprFpr(metric_dict["processed"][each]) for each in metric_dict["raw"]}
            random_tpr_fpr_dict = {each: metric_dict["baseline"][each][3:] for each in metric_dict["baseline"]}
            raw_tpr_fpr_vis = process4TprFprVis(raw_best_tpr_fpr, random_tpr_fpr_dict)
            processed_tpr_fpr_vis = process4TprFprVis(processed_best_tpr_fpr, random_tpr_fpr_dict)
            visTprFpr(raw_tpr_fpr_vis, all_models_list+["random"], title="{}-{}-raw".format(sim_name, num_genes))
            visTprFpr(processed_tpr_fpr_vis, all_models_list+["random"], title="{}-{}-processed".format(sim_name, num_genes))
            # Compare AUROC
            raw_best_auroc = {each: getModelPerformance(metric_dict["raw"][each], "AUROC", func="max") for each in
                           metric_dict["raw"]}
            processed_best_auroc = {each: getModelPerformance(metric_dict["processed"][each], "AUROC", func="max") for each in
                                 metric_dict["processed"]}
            random_auroc_dict = {each: metric_dict["baseline"][each][1] for each in metric_dict["baseline"]}
            raw_auroc_vis = process4Vis(raw_best_auroc, random_auroc_dict)
            processed_auroc_vis = process4Vis(processed_best_auroc, random_auroc_dict)
            filename_list = list(raw_best_auroc.keys())
            visBox4Metric(raw_auroc_vis, all_models_list + ["random"],
                          title="{}-{}-raw".format(sim_name, num_genes), metric_name="AUROC")
            visBox4Metric(processed_auroc_vis, all_models_list + ["random"],
                          title="{}-{}-processed".format(sim_name, num_genes), metric_name="AUROC")
            # Compare AUPR
            raw_best_aupr = {each: getModelPerformance(metric_dict["raw"][each], "AUPRC", func="max") for each in
                              metric_dict["raw"]}
            processed_best_aupr = {each: getModelPerformance(metric_dict["processed"][each], "AUPRC", func="max") for
                                    each in metric_dict["processed"]}
            random_aupr_dict = {each: metric_dict["baseline"][each][1] for each in metric_dict["baseline"]}
            raw_aupr_vis = process4Vis(raw_best_aupr, random_aupr_dict)
            processed_aupr_vis = process4Vis(processed_best_auroc, random_auroc_dict)
            filename_list = list(raw_best_aupr.keys())
            visBox4Metric(raw_aupr_vis, all_models_list + ["random"],
                          title="{}-{}-raw".format(sim_name, num_genes), metric_name="AUPR")
            visBox4Metric(processed_aupr_vis, all_models_list + ["random"],
                          title="{}-{}-processed".format(sim_name, num_genes), metric_name="AUPR")


if __name__ == '__main__':
    metric_dir_path = "../res/logs/evaluation/simulation/"
    compareModelMetrics(metric_dir_path, data_name="Cortex")
    compareModelMetrics(metric_dir_path, data_name="pbmc")