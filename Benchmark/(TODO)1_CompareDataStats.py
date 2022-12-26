import matplotlib.pyplot as plt
import torch
import pandas as pd
import numpy as np
import itertools
import scipy.stats
import os
from util.BenchmarkUtils import loadExprMat, loadNetMat
from util.Visualization import visHist4NewSim, vis4NewSim, vis4OldSim, visHist4OldSim, vis4Corr, vis4SERGIOSim, visHist4SERGIOSim


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

def computeStatSim(stat1, stat2):
    # Compute similarity score between basic stats
    pass


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

new_simulation_list = ["NORTA"]  # ["NORTA", "netSmooth"]
old_simulation_list = ["scLink", "ZILGM"]


def compute4NewSim(data_name):
    sim_dir_path = "../data/simulated/new/"
    save_dir_path = "../res/logs/comparison/"
    if data_name == "PBMC":
        ref_dir_path = "../data/experimental/PBMC/processed/expr/"
        filename_list = ["-".join(each) for each in
                         itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"], ["50hvg", "100hvg"])]
    elif data_name == "Cortex":
        ref_dir_path = "../data/experimental/mouse_cortex/processed/expr/"
        filename_list = ["-".join(each) for each in
                         itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"], ["50hvg", "100hvg"])]
    else:
        raise ValueError("Unknown data name {}!".format(data_name))
    bins = 20
    file_stats_dist = {}
    file_corr_dist = {}
    file_hist_dict = {}
    for each_file in filename_list:
        stats_dict = {}
        corr_dict = {}
        hist_dict = {}
        print("=" * 70)
        print("[ {} ]".format(each_file))
        net_mat = loadNetMat("{}/{}-net_mat.csv".format(sim_dir_path, each_file)).values
        corr_dict["sim_net"] = net_mat
        # Experimental data
        print("Reference data...")
        ref_data_mat = loadExprMat("{}/{}.csv".format(ref_dir_path, each_file)).values
        log_ref_data_mat = np.log10(ref_data_mat + 1)
        stats_dict["ref"] = computeBasicStats(ref_data_mat)
        corr_dict["ref"] = computeCorr(ref_data_mat)
        hist_dict["ref"] = computeHist(_minMaxScale(ref_data_mat), bins)
        stats_dict["log_ref"] = computeBasicStats(log_ref_data_mat)
        corr_dict["log_ref"] = computeCorr(log_ref_data_mat)
        hist_dict["log_ref"] = computeHist(_minMaxScale(log_ref_data_mat), bins)
        # Simulated data
        for simulation_type in new_simulation_list:
            print("{} simulation data...".format(simulation_type))
            sim_data_mat = loadExprMat("{}/{}-{}-data_mat.csv".format(sim_dir_path, each_file, simulation_type)).values
            log_sim_data_mat = np.log10(sim_data_mat + 1)
            stats_dict[simulation_type] = computeBasicStats(sim_data_mat)
            corr_dict[simulation_type] = computeCorr(sim_data_mat)
            hist_dict[simulation_type] = computeHist(_minMaxScale(sim_data_mat), bins)
            stats_dict["log_{}".format(simulation_type)] = computeBasicStats(log_sim_data_mat)
            corr_dict["log_{}".format(simulation_type)] = computeCorr(log_sim_data_mat)
            hist_dict["log_{}".format(simulation_type)] = computeHist(_minMaxScale(log_sim_data_mat), bins)
        file_stats_dist[each_file] = stats_dict
        file_corr_dist[each_file] = corr_dict
        file_hist_dict[each_file] = hist_dict
    # Save data
    np.save("{}/{}-new_sim-basic_stats.npy".format(save_dir_path, data_name), file_stats_dist)
    np.save("{}/{}-new_sim-corr.npy".format(save_dir_path, data_name), file_corr_dist)
    np.save("{}/{}-new_sim-hist.npy".format(save_dir_path, data_name), file_hist_dict)


def compute4OldSim():
    sim_dir_path = "../data/simulated/old/"
    save_dir_path = "../res/logs/comparison/"
    num_genes_list = ["50hvg", "100hvg"]
    rho = "0.070000"
    bins = 20
    file_stats_dict = {}
    file_corr_dict = {}
    file_hist_dict = {}
    for num_genes in num_genes_list:
        stats_dict = {}
        corr_dict = {}
        hist_dict = {}
        print("=" * 70)
        print("[ {} ]".format(num_genes))
        # Simulated data
        net_mat = loadNetMat("{}/{}-net_mat.csv".format(sim_dir_path, num_genes)).values
        for simulation_type in old_simulation_list:
            print("{} simulation data...".format(simulation_type))
            if simulation_type == "scLink":
                log_sim_data_mat = loadExprMat(
                    "{}/{}-{}_rho{}-data_mat.csv".format(sim_dir_path, num_genes, simulation_type, rho)).values
                non_zero_idx = np.where(log_sim_data_mat != 0)
                tmp = np.power(10, log_sim_data_mat[non_zero_idx]) - 1
                sim_data_mat = np.zeros_like(log_sim_data_mat)
                sim_data_mat[non_zero_idx] = tmp
            else:
                sim_data_mat = loadExprMat(
                    "{}/{}-{}_rho{}-data_mat.csv".format(sim_dir_path, num_genes, simulation_type, rho)).values
                log_sim_data_mat = np.log10(sim_data_mat + 1)
            stats_dict[simulation_type] = computeBasicStats(sim_data_mat)
            stats_dict["log_{}".format(simulation_type)] = computeBasicStats(log_sim_data_mat)
            corr_dict[simulation_type] = {
                "true_net": net_mat,
                "sim_net": computeCorr(sim_data_mat)
            }
            corr_dict["log_{}".format(simulation_type)] = {
                "true_net": net_mat,
                "sim_net": computeCorr(log_sim_data_mat)
            }
            hist_dict[simulation_type] = computeHist(_minMaxScale(sim_data_mat), bins)
            hist_dict["log_{}".format(simulation_type)] = computeHist(_minMaxScale(log_sim_data_mat), bins)
        file_stats_dict[num_genes] = stats_dict
        file_corr_dict[num_genes] = corr_dict
        file_hist_dict[num_genes] = hist_dict
    # Save data
    np.save("{}/old_sim-basic_stats.npy".format(save_dir_path), file_stats_dict)
    np.save("{}/old_sim-corr.npy".format(save_dir_path), file_corr_dict)
    np.save("{}/old_sim-hist.npy".format(save_dir_path), file_hist_dict)


def compute4SERGIO():
    save_dir_path = "../res/logs/comparison/"
    bins = 20
    # -----
    simulation_type = "SERGIO"
    print("=" * 70)
    # Simulated data
    file_stats_dict = {}
    file_hist_dict = {}
    for sparsity_level in [1, 5, 10, 15, 20]:
        stats_dict = {}
        hist_dict = {}
        # filename = "../data/SERGIO_simulation/100gene-9groups-{}sparsity.csv".format(sparsity_level)
        filename = "../data/SERGIO_simulation_all/100gene-9groups-{}sparsity.csv".format(sparsity_level)
        sim_data_mat = loadExprMat(filename).values
        log_sim_data_mat = np.log10(sim_data_mat + 1)
        stats_dict[simulation_type] = computeBasicStats(sim_data_mat)
        stats_dict["log_{}".format(simulation_type)] = computeBasicStats(log_sim_data_mat)
        hist_dict[simulation_type] = computeHist(_minMaxScale(sim_data_mat), bins)
        hist_dict["log_{}".format(simulation_type)] = computeHist(_minMaxScale(log_sim_data_mat), bins)
        # -----
        file_stats_dict[sparsity_level] = stats_dict
        file_hist_dict[sparsity_level] = hist_dict
    # Save data
    np.save("{}/SERGIO_sim-basic_stats.npy".format(save_dir_path), file_stats_dict)
    np.save("{}/SERGIO_sim-hist.npy".format(save_dir_path), file_hist_dict)

# ======================================

def statsSummary(data_name):
    sim_dir_path = "../data/simulated/new/"
    if data_name == "PBMC":
        ref_dir_path = "../data/experimental/PBMC/processed/expr/"
        filename_list = ["-".join(each) for each in
                         itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"], ["50hvg", "100hvg"])]
    elif data_name == "Cortex":
        ref_dir_path = "../data/experimental/mouse_cortex/processed/expr/"
        filename_list = ["-".join(each) for each in
                         itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"], ["50hvg", "100hvg"])]
    else:
        raise ValueError("Unknown data name {}!".format(data_name))
    for each_file in filename_list:
        print("=" * 70)
        print("[ {} ]".format(each_file))
        net_mat = loadNetMat("{}/{}-net_mat.csv".format(sim_dir_path, each_file)).values
        # Experimental data
        print("Reference data...")
        ref_data_mat = loadExprMat("{}/{}.csv".format(ref_dir_path, each_file)).values.reshape(-1)
        print("Sparsity = {} | Max = {} | Dispersion = {}".format(np.count_nonzero(ref_data_mat) / len(ref_data_mat),
                                                                  np.max(ref_data_mat),
                                                                  np.var(ref_data_mat) / np.mean(ref_data_mat)))
        print()
        # Simulated data
        for simulation_type in new_simulation_list:
            print("{} simulation data...".format(simulation_type))
            sim_data_mat = loadExprMat(
                "{}/{}-{}-data_mat.csv".format(sim_dir_path, each_file, simulation_type)).values.reshape(-1)
            log_sim_data_mat = np.log10(sim_data_mat + 1)
            print(
                "Sparsity = {} | Max = {} | Dispersion = {}".format(np.count_nonzero(sim_data_mat) / len(sim_data_mat),
                                                                    np.max(sim_data_mat),
                                                                    np.var(sim_data_mat) / np.mean(sim_data_mat)))

# ======================================


def compareNewSim(data_name):
    res_dir_path = "../res/logs/comparison/"
    data_dir_path = "../data/simulated/new/"
    stats_dict = np.load("{}/{}-new_sim-basic_stats.npy".format(res_dir_path, data_name), allow_pickle=True).item()
    corr_dict = np.load("{}/{}-new_sim-corr.npy".format(res_dir_path, data_name), allow_pickle=True).item()
    hist_dict = np.load("{}/{}-new_sim-hist.npy".format(res_dir_path, data_name), allow_pickle=True).item()
    # -----
    # exp_data_name = (
    #         ["-".join(list(each)) for each in
    #          itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"], ["50hvg", "100hvg"])] +
    #         ["-".join(list(each)) for each in itertools.product(["pbmc1", "pbmc2"], ["50hvg", "100hvg"])]
    # )
    pcc_CMD = []
    pcc_MSE = []
    for each_file in stats_dict.keys():
        print("[ {} ]".format(each_file))
        # # Compare basic stats
        vis4NewSim(stats_dict[each_file], data_name=each_file)
        # Compare correlation similarity
        CD, diff_MSE, (before_CMD, after_CMD), (before_MSE, after_MSE) = computeCorrSim(corr_dict[each_file],
                                                                                        type="PCC")
        pcc_CMD.append((before_CMD, after_CMD))
        pcc_MSE.append((before_MSE, after_MSE))
        print("| PCC | CD = {}, Diff MSE = {}".format(CD, diff_MSE))
        # Compare hist
        visHist4NewSim(hist_dict[each_file], title=each_file, filename="")
    vis4Corr(pcc_CMD, pcc_MSE, list(stats_dict.keys()))


def compareOldSim():
    # Compare simulation with 8 experimental datasets
    res_dir_path = "../res/logs/comparison/"
    cortex_stats_dict = np.load("{}/{}-new_sim-basic_stats.npy".format(res_dir_path, "Cortex"),
                                allow_pickle=True).item()
    cortex_hist_dict = np.load("{}/{}-new_sim-hist.npy".format(res_dir_path, "Cortex"), allow_pickle=True).item()
    pbmc_stats_dict = np.load("{}/{}-new_sim-basic_stats.npy".format(res_dir_path, "PBMC"), allow_pickle=True).item()
    pbmc_hist_dict = np.load("{}/{}-new_sim-hist.npy".format(res_dir_path, "PBMC"), allow_pickle=True).item()
    sim_stats_dict = np.load("{}/old_sim-basic_stats.npy".format(res_dir_path), allow_pickle=True).item()
    sim_hist_dict = np.load("{}/old_sim-hist.npy".format(res_dir_path), allow_pickle=True).item()
    # -----
    true_stats_dict = {**cortex_stats_dict, **pbmc_stats_dict}
    true_hist_dict = {**cortex_hist_dict, **pbmc_hist_dict}
    for num_genes in ["50hvg", "100hvg"]:
        tmp_true_stats_dict = {each: true_stats_dict[each] for each in true_stats_dict if num_genes in each}
        tmp_true_hist_dict = {each: true_hist_dict[each] for each in true_hist_dict if num_genes in each}
        # Compare basic stats
        vis4OldSim(tmp_true_stats_dict, sim_stats_dict[num_genes], data_name=num_genes)
        # Compare hist
        visHist4OldSim(tmp_true_hist_dict, sim_hist_dict[num_genes], title=num_genes, filename="")


def compareSERGIO():
    res_dir_path = "../res/logs/comparison/"
    cortex_stats_dict = np.load("{}/{}-new_sim-basic_stats.npy".format(res_dir_path, "Cortex"), allow_pickle=True).item()
    cortex_hist_dict = np.load("{}/{}-new_sim-hist.npy".format(res_dir_path, "Cortex"), allow_pickle=True).item()
    pbmc_stats_dict = np.load("{}/{}-new_sim-basic_stats.npy".format(res_dir_path, "PBMC"), allow_pickle=True).item()
    pbmc_hist_dict = np.load("{}/{}-new_sim-hist.npy".format(res_dir_path, "PBMC"), allow_pickle=True).item()
    sim_stats_dict = np.load("{}/SERGIO_sim-basic_stats.npy".format(res_dir_path), allow_pickle=True).item()
    sim_hist_dict = np.load("{}/SERGIO_sim-hist.npy".format(res_dir_path), allow_pickle=True).item()
    # -----
    true_stats_dict = {**cortex_stats_dict, **pbmc_stats_dict}
    true_hist_dict = {**cortex_hist_dict, **pbmc_hist_dict}
    for sparsity_level in [1, 5, 10, 15, 20]:
    # for sparsity_level in [20]:
        tmp_true_stats_dict = {each: true_stats_dict[each] for each in true_stats_dict if "100hvg" in each}
        tmp_true_hist_dict = {each: true_hist_dict[each] for each in true_stats_dict if "100hvg" in each}
        # Compare basic stats
        vis4SERGIOSim(tmp_true_stats_dict, sim_stats_dict[sparsity_level], data_name=sparsity_level)
        # Compare hist
        visHist4SERGIOSim(tmp_true_hist_dict, sim_hist_dict[sparsity_level], title=sparsity_level, filename="")


# ======================================

if __name__ == '__main__':
    # # -----
    # compute4NewSim("PBMC")
    # compute4NewSim("Cortex")
    # compute4OldSim()
    # compute4SERGIO()
    # # -----
    # compareNewSim("Cortex")
    compareNewSim("PBMC")
    # compareSERGIO()
    # # -----
    # compareOldSim()
    # # -----
    # statsSummary("Cortex")
    # statsSummary("PBMC")
    pass
