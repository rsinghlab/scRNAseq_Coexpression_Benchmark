import os.path

import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import rpy2
from rpy2.robjects import numpy2ri, pandas2ri
from sklearn.metrics import confusion_matrix, f1_score

# ===================================

def getGeneList(data_name):
    if "pbmc" in data_name:
        data_dir_path = "../data/experimental/PBMC/processed/expr/"
    elif "Cortex" in data_name:
        data_dir_path = "../data/experimental/mouse_cortex/processed/expr/"
    else:
        raise ValueError("Unknown  data name {}!".format(data_name))
    with open("{}/{}-1000hvg.csv".format(data_dir_path, data_name)) as f:
        first_line = f.readline().strip('\n').strip(",")
    gene_list = first_line.split(",")
    gene_list = [each.split("_")[-1] for each in gene_list]
    return gene_list


def getTabulaGeneList(data_name):
    data_dir_path = "../data/Tabula_Muris/500hvg/"
    with open("{}/{}.csv".format(data_dir_path, data_name)) as f:
        first_line = f.readline().strip('\n').strip(",")
    gene_list = first_line.split(",")
    # gene_list = [each.split("_")[-1] for each in gene_list]
    return gene_list


def loadRefNet(filename):
    net = pd.read_csv(filename, header=0, index_col=0)
    return net


def loadEstNetList(filename, model_name):
    if not os.path.isfile(filename):
        return None
    numpy2ri.activate()
    pandas2ri.activate()
    r = robjects.r
    r['source']('./loadRDS.R')
    model_function_r = robjects.globalenv['loadPred']
    net_list = model_function_r(filename, model_name)
    net_list = [each for each in net_list]
    numpy2ri.deactivate()
    pandas2ri.deactivate()
    return net_list

# ===================================

def _evaluate(est_net, ref_net):
    assert est_net.shape == ref_net.shape
    num_genes = est_net.shape[0]
    abs_ref_net = np.sign(np.abs(ref_net))
    abs_est_net = np.sign(np.abs(est_net))
    triu_idx = np.triu_indices(num_genes, 1)
    est_triu = abs_est_net[triu_idx]
    ref_triu = abs_ref_net[triu_idx]
    # -----
    f1 = f1_score(ref_triu, est_triu)
    conf_mat = confusion_matrix(ref_triu, est_triu)
    TP = conf_mat[1, 1]
    FP = conf_mat[0, 1]
    TN = conf_mat[0, 0]
    FN = conf_mat[1, 0]
    TPR = TP/(TP+FN)
    precision=TP/(TP+FP)
    FDR=1-precision
    return {
        "conf_mat": conf_mat, "TPR": TPR, "precision": precision, "FDR": FDR, "f1": f1,
        "edge_num":np.count_nonzero(ref_triu), "est_num": np.count_nonzero(est_triu)}


def seqEvalFull(est_net_list, ref_net, ref_thr):
    ref_net = ref_net.values
    ref_net[np.abs(ref_net) < ref_thr] = 0.0
    seq_metric = [_evaluate(each, ref_net) for each in est_net_list]
    return seq_metric


def findBestRes(seq_metric):
    seq_metric = [each for each in seq_metric if not np.isnan(each["precision"])]
    best_precision_idx = np.argmax([each["precision"] for each in seq_metric])
    best_f1_idx = np.argmax([each["f1"] for each in seq_metric])
    best_tpr_idx = np.argmax([each["TPR"] for each in seq_metric])
    # -----
    best_precision = seq_metric[best_precision_idx]
    best_f1 = seq_metric[best_f1_idx]
    best_tpr = seq_metric[best_tpr_idx]
    return best_precision, best_f1, best_tpr


if __name__ == '__main__':
    pbmc_sim_list = [
        "pbmc1-Drop", "pbmc2-Drop", "pbmc1-inDrops", "pbmc2-inDrops"
    ]
    cortex_sim_list = [
        "Cortex1-10xChromium", "Cortex2-10xChromium", "Cortex1-Smart_seq2", "Cortex2-Smart_seq2"
    ]
    tabula_data_list = ["skeletal_muscle_satellite_stem_cell", "T_cell", "type_B_pancreatic_cell"]
    model_list = ["pearson", "spearman", "glasso", "scLink", "GENIE3", "scDesign2", "minet", "PIDC"]
    ref_thr = 0.2
    # # -----
    # # Pathway Common similarity network
    # all_data_metric = {}
    # for d in pbmc_sim_list:
    #     print("=" * 70)
    #     print("[ {} ]".format(d))
    #     print("Loading gene list and reference network...")
    #     gene_list = getGeneList(d)
    #     ref_net = loadRefNet("../data/ChIP_seq_data/{}-Pathway_similarity-net_mat.csv".format(d))
    #     model_metric = {}
    #     for m in model_list:
    #         print("Model: {}".format(m))
    #         if m == "PIDC":
    #             est_net_list = loadEstNetList("../res/prediction/real/{}-1000hvg-STRING-PIDC-predictions.rds".format(d), m)
    #         else:
    #             est_net_list = loadEstNetList("../res/prediction/real/{}-1000hvg-STRING-predictions.rds".format(d), m)
    #         if est_net_list is None:
    #             continue
    #         seq_metric = seqEvalFull(est_net_list, ref_net, ref_thr=0.5)
    #         best_precision, best_f1, best_tpr = best_metric = findBestRes(seq_metric)
    #         model_metric[m] = {"precision": best_precision, "f1": best_f1, "TPR": best_tpr}
    #     all_data_metric[d] = model_metric
    # np.save("./Pathway_similarity-pbmc-metric.npy", all_data_metric)
    # -----
    # Mouse TF co-expression network
    all_data_metric = {}
    for d in cortex_sim_list:
        print("=" * 70)
        print("[ {} ]".format(d))
        print("Loading gene list and reference network...")
        gene_list = getGeneList(d)
        ref_net = loadRefNet("../data/ChIP_seq_data/{}-mouse_TF_PCC-sub_net_mat.csv".format(d))
        model_metric = {}
        for m in model_list:
            print("Model: {}".format(m))
            if m == "PIDC":
                est_net_list = loadEstNetList("../res/prediction/real/{}-1000hvg-STRING-PIDC-predictions.rds".format(d), m)
            else:
                est_net_list = loadEstNetList("../res/prediction/real/{}-1000hvg-STRING-predictions.rds".format(d), m)
            if est_net_list is None:
                print("No such file")
                continue
            # -----
            tf_list = ref_net.columns.values
            print("Num of TFs={}".format(len(tf_list)))
            gene_idx = [i for i, g in enumerate(gene_list) if g in tf_list]
            est_net_list = [each[gene_idx, :][:, gene_idx] for each in est_net_list]
            # -----
            seq_metric = seqEvalFull(est_net_list, ref_net, ref_thr=ref_thr)
            best_precision, best_f1, best_tpr = best_metric = findBestRes(seq_metric)
            model_metric[m] = {"precision": best_precision, "f1": best_f1, "TPR": best_tpr}
        all_data_metric[d] = model_metric
    np.save("./moust_TF_PCC-cortex-metric.npy", all_data_metric)
    # -----
    # Human TF co-expression network
    all_data_metric = {}
    for d in pbmc_sim_list:
        print("=" * 70)
        print("[ {} ]".format(d))
        print("Loading gene list and reference network...")
        gene_list = getGeneList(d)
        ref_net = loadRefNet("../data/ChIP_seq_data/{}-human_TF_similarity-sub_net_mat.csv".format(d))
        model_metric = {}
        for m in model_list:
            print("Model: {}".format(m))
            if m == "PIDC":
                est_net_list = loadEstNetList("../res/prediction/real/{}-1000hvg-STRING-PIDC-predictions.rds".format(d), m)
            else:
                est_net_list = loadEstNetList("../res/prediction/real/{}-1000hvg-STRING-predictions.rds".format(d), m)
            if est_net_list is None:
                print("No such file")
                continue
            # -----
            tf_list = ref_net.columns.values
            print("Num of TFs={}".format(len(tf_list)))
            gene_idx = [i for i, g in enumerate(gene_list) if g in tf_list]
            est_net_list = [each[gene_idx, :][:, gene_idx] for each in est_net_list]
            # -----
            seq_metric = seqEvalFull(est_net_list, ref_net, ref_thr=ref_thr)
            best_precision, best_f1, best_tpr = best_metric = findBestRes(seq_metric)
            model_metric[m] = {"precision": best_precision, "f1": best_f1, "TPR": best_tpr}
        all_data_metric[d] = model_metric
    np.save("./human_TF_similarity-cortex-metric.npy", all_data_metric)
    # -----
    # Tabula Muris TF co-expression network
    all_data_metric = {}
    for d in tabula_data_list:
        print("=" * 70)
        print("[ {} ]".format(d))
        print("Loading gene list and reference network...")
        gene_list = getTabulaGeneList(d)
        ref_net = loadRefNet("../data/ChIP_seq_data/{}-mouse_TF_PCC-sub_net_mat.csv".format(d))
        model_metric = {}
        for m in model_list:
            print("Model: {}".format(m))
            # if m == "PIDC":
            #     est_net_list = loadEstNetList("../res/prediction/real/{}-STRING-PIDC-predictions.rds".format(d), m)
            # else:
            est_net_list = loadEstNetList("../res/prediction/real/{}-STRING-predictions.rds".format(d), m)
            if est_net_list is None:
                print("No such file")
                continue
            # -----
            tf_list = ref_net.columns.values
            print("Num of TFs={}".format(len(tf_list)))
            gene_idx = [i for i, g in enumerate(gene_list) if g in tf_list]
            est_net_list = [each[gene_idx, :][:, gene_idx] for each in est_net_list]
            # -----
            seq_metric = seqEvalFull(est_net_list, ref_net, ref_thr=ref_thr)
            best_precision, best_f1, best_tpr = best_metric = findBestRes(seq_metric)
            model_metric[m] = {"precision": best_precision, "f1": best_f1, "TPR": best_tpr}
        all_data_metric[d] = model_metric
    np.save("./tabula_TF_PCC-cortex-metric.npy", all_data_metric)
