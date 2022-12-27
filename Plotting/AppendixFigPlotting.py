'''
Description:
    Paper appendix figures plotting.

Author:
    Jiaqi Zhang <jiaqi_zhang2@brown.edu>
'''
import pandas as pd
import numpy as np
import itertools
import scipy.stats
import scipy
from PlottingUtils import *

all_models_list = ["pearson", "spearman", "glasso", "GENIE3", "minet", "PIDC", "scDesign2", "scLink", "ZILGM"]
bulk_model_list = ["pearson", "spearman", "glasso", "GENIE3", "minet"]
sc_model_list = ["PIDC", "scDesign2", "scLink", "ZILGM"]

# ========================================================
# Fig A7: method performances on NORTA simulation (500 genes)
# ========================================================
def NORTA500Genes():
    print("Model Performance on NORTA Simulation w/ 500 HVGS")
    metric = np.load("../fig_res/high_dim_metric.npy", allow_pickle=True).item()
    F1_metric = metric["F1"]
    FDR_metric = metric["FDR"]
    AUROC_metric = metric["AUROC"]
    AUPRC_metric = metric["AUPRC"]
    # -----
    model_F1_list = [
        [F1_metric[x][m] for x in F1_metric if (m in F1_metric[x].index) and (not pd.isna(F1_metric[x][m]))] for m in
        all_models_list]
    model_FDR_list = [
        [FDR_metric[x][m] for x in FDR_metric if (m in FDR_metric[x].index) and (not pd.isna(FDR_metric[x][m]))] for m
        in all_models_list]
    model_AUROC_list = [
        [AUROC_metric[x][m] for x in AUROC_metric if (m in AUROC_metric[x].index) and (not pd.isna(AUROC_metric[x][m]))] for m
        in all_models_list]
    model_AUPRC_list = [
        [AUPRC_metric[x][m] for x in AUPRC_metric if (m in AUPRC_metric[x].index) and (not pd.isna(AUPRC_metric[x][m]))] for m
        in all_models_list]
    # -----
    print(all_models_list)
    print("Avg F1: ", [np.mean(each) for each in model_F1_list])
    print("Std F1: ", [np.std(each) for each in model_F1_list])
    print("Avg FDR: ", [np.mean(each) for each in model_FDR_list])
    print("Std FDR: ", [np.std(each) for each in model_FDR_list])
    print("Avg AUROC: ", [np.mean(each) for each in model_AUROC_list])
    print("Std AUROC: ", [np.std(each) for each in model_AUROC_list])
    print("Avg AUPRC: ", [np.mean(each) for each in model_AUPRC_list])
    print("Std AUPRC: ", [np.std(each) for each in model_AUPRC_list])
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_F1_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1 score")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.5, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03]
    plt.ylim(0.0, 0.5)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA7/500HVG-F1_score.pdf")
    plt.show()
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_AUROC_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("AUROC")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.9, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03]
    plt.ylim(0.4, 0.9)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA7/500HVG-AUROC.pdf")
    plt.show()
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_AUPRC_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("AUPRC")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.5, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03]
    # plt.ylim(0.0, 0.5)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA7/500HVG-AUPRC.pdf")
    plt.show()
    # -----
    model_avg_FDR = np.asarray([np.mean(each) for each in model_FDR_list])
    model_std_FDR = np.asarray([np.std(each) for each in model_FDR_list])
    plt.figure(figsize=(8, 5))
    width = 0.6
    plt.bar(np.arange(len(all_models_list)), model_avg_FDR, label="FDR", width=width, color=BlueRed_12.mpl_colors[0])
    plt.bar(np.arange(len(all_models_list)), 1 - model_avg_FDR, label="PPV", width=width, bottom=model_avg_FDR,
            color=BlueRed_12.mpl_colors[-2])
    plt.xticks(np.arange(len(all_models_list)), [model_name_dict[x] for x in all_models_list], fontsize=13)
    # plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    # plt.legend(bbox_to_anchor=(1.0, 0.9))
    plt.legend(bbox_to_anchor=(0.1, 1.0), ncol=2)
    plt.tight_layout()
    plt.savefig("../fig/FigA7/500HVG-FDR.pdf")
    plt.show()

# ========================================================
# Fig A8: method performances on SERGIO simulation (400 genes)
# ========================================================
def SERGIO400Genes():
    print("Model Performance on SERGIO Simulation w/ 400 HVGS")
    metric = np.load("../fig_res/SERGIO_400genes_metric.npy", allow_pickle=True).item()
    F1_metric = metric["F1"]
    FDR_metric = metric["FDR"]
    AUROC_metric = metric["AUROC"]
    AUPRC_metric = metric["AUPRC"]
    # -----
    model_F1_list = [
        [F1_metric[x][m] for x in F1_metric if (m in F1_metric[x].index) and (not pd.isna(F1_metric[x][m]))] for m in
        all_models_list]
    model_FDR_list = [
        [FDR_metric[x][m] for x in FDR_metric if (m in FDR_metric[x].index) and (not pd.isna(FDR_metric[x][m]))] for m
        in all_models_list]
    model_AUROC_list = [
        [AUROC_metric[x][m] for x in AUROC_metric if (m in AUROC_metric[x].index) and (not pd.isna(AUROC_metric[x][m]))] for m
        in all_models_list]
    model_AUPRC_list = [
        [AUPRC_metric[x][m] for x in AUPRC_metric if (m in AUPRC_metric[x].index) and (not pd.isna(AUPRC_metric[x][m]))] for m
        in all_models_list]
    # -----
    print(all_models_list)
    print("Avg F1: ", [np.mean(each) for each in model_F1_list])
    print("Std F1: ", [np.std(each) for each in model_F1_list])
    print("Avg FDR: ", [np.mean(each) for each in model_FDR_list])
    print("Std FDR: ", [np.std(each) for each in model_FDR_list])
    print("Avg AUROC: ", [np.mean(each) for each in model_AUROC_list])
    print("Std AUROC: ", [np.std(each) for each in model_AUROC_list])
    print("Avg AUPRC: ", [np.mean(each) for each in model_AUPRC_list])
    print("Std AUPRC: ", [np.std(each) for each in model_AUPRC_list])
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_F1_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1 score")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.05, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03]
    plt.ylim(0.02, 0.05)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.001, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.001, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA8/400HVG-F1_score.pdf")
    plt.show()
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_AUROC_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("AUROC")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.6, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03]
    plt.ylim(0.4, 0.6)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA8/400HVG-AUROC.pdf")
    plt.show()
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_AUPRC_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("AUPRC")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.04, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03]
    plt.ylim(0.0, 0.04)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.001, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.001, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA8/400HVG-AUPRC.pdf")
    plt.show()
    # -----
    model_avg_FDR = np.asarray([np.mean(each) for each in model_FDR_list])
    model_std_FDR = np.asarray([np.std(each) for each in model_FDR_list])
    plt.figure(figsize=(8, 5))
    width = 0.6
    plt.bar(np.arange(len(all_models_list)), model_avg_FDR, label="FDR", width=width, color=BlueRed_12.mpl_colors[0])
    plt.bar(np.arange(len(all_models_list)), 1 - model_avg_FDR, label="PPV", width=width, bottom=model_avg_FDR,
            color=BlueRed_12.mpl_colors[-2])
    plt.xticks(np.arange(len(all_models_list)), [model_name_dict[x] for x in all_models_list], fontsize=13)
    # plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    # plt.legend(bbox_to_anchor=(1.0, 0.9))
    plt.legend(bbox_to_anchor=(0.1, 1.0), ncol=2)
    plt.tight_layout()
    plt.savefig("../fig/FigA8/400HVG-FDR.pdf")
    plt.show()

# ========================================================
# Fig A9: method performances on ZI-Gaussian simulations
# ========================================================
def performanceZIGaussian():
    print("Model Performance on Old ZI-Gaussian Data")
    metric = np.load("../fig_res/ZIGaussian_sim_metric.npy", allow_pickle=True).item()["100hvg"]
    F1_metric = metric["F1"]
    FDR_metric = metric["FDR"]
    AUROC_metric = metric["AUROC"]
    AUPRC_metric = metric["AUPRC"]
    # -----
    model_F1_list = [
        [F1_metric[x][m] for x in F1_metric if (m in F1_metric[x].index) and (not pd.isna(F1_metric[x][m]))] for m in
        all_models_list]
    model_FDR_list = [
        [FDR_metric[x][m] for x in FDR_metric if (m in FDR_metric[x].index) and (not pd.isna(FDR_metric[x][m]))] for m
        in all_models_list]
    model_AUROC_list = [
        [AUROC_metric[x][m] for x in AUROC_metric if (m in AUROC_metric[x].index) and (not pd.isna(AUROC_metric[x][m]))]
        for m
        in all_models_list]
    model_AUPRC_list = [
        [AUPRC_metric[x][m] for x in AUPRC_metric if (m in AUPRC_metric[x].index) and (not pd.isna(AUPRC_metric[x][m]))]
        for m
        in all_models_list]
    # -----
    print(all_models_list)
    print("Avg F1: ", [np.mean(each) for each in model_F1_list])
    print("Std F1: ", [np.std(each) for each in model_F1_list])
    print("Avg FDR: ", [np.mean(each) for each in model_FDR_list])
    print("Std FDR: ", [np.std(each) for each in model_FDR_list])
    print("Avg AUROC: ", [np.mean(each) for each in model_AUROC_list])
    print("Std AUROC: ", [np.std(each) for each in model_AUROC_list])
    print("Avg AUPRC: ", [np.mean(each) for each in model_AUPRC_list])
    print("Std AUPRC: ", [np.std(each) for each in model_AUPRC_list])
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_F1_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1 score")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.5, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03])
    plt.ylim(0.0, 0.5)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA9/ZIGaussian-F1_score.pdf")
    plt.show()
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_AUROC_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("AUROC")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.8, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03])
    # plt.ylim(0.0, 0.5)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA9/ZIGaussian-AUROC.pdf")
    plt.show()
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_AUPRC_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("AUPRC")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.4, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03])
    plt.ylim(0.0, 0.4)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA9/ZIGaussian-AUPRC.pdf")
    plt.show()
    # -----
    model_avg_FDR = np.asarray([np.mean(each) for each in model_FDR_list])
    model_std_FDR = np.asarray([np.std(each) for each in model_FDR_list])
    plt.figure(figsize=(8, 5))
    width = 0.6
    plt.bar(np.arange(len(all_models_list)), model_avg_FDR, label="FDR", width=width, color=BlueRed_12.mpl_colors[0])
    plt.bar(np.arange(len(all_models_list)), 1 - model_avg_FDR, label="PPV", width=width, bottom=model_avg_FDR,
            color=BlueRed_12.mpl_colors[-2])
    plt.xticks(np.arange(len(all_models_list)), [model_name_dict[x] for x in all_models_list], fontsize=13)
    # plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    # plt.legend(bbox_to_anchor=(1.0, 0.9))
    plt.legend(bbox_to_anchor=(0.1, 1.0), ncol=2)
    plt.tight_layout()
    plt.savefig("../fig/FigA9/ZIGaussian-FDR.pdf")
    plt.show()

# ========================================================
# Fig A10: method performances on ZI-Poisson simulations
# ========================================================
def performanceZIPoisson():
    print("Model Performance on ZI-Poisson Data")
    metric = np.load("../fig_res/ZIPoisson_sim_metric.npy", allow_pickle=True).item()["100hvg"]
    F1_metric = metric["F1"]
    FDR_metric = metric["FDR"]
    # -----
    model_F1_list = [
        [F1_metric[x][m] for x in F1_metric if (m in F1_metric[x].index) and (not pd.isna(F1_metric[x][m]))] for m in
        all_models_list]
    model_FDR_list = [
        [FDR_metric[x][m] for x in FDR_metric if (m in FDR_metric[x].index) and (not pd.isna(FDR_metric[x][m]))] for m
        in all_models_list]
    AUROC_metric = metric["AUROC"]
    AUPRC_metric = metric["AUPRC"]
    # -----
    model_F1_list = [
        [F1_metric[x][m] for x in F1_metric if (m in F1_metric[x].index) and (not pd.isna(F1_metric[x][m]))] for m in
        all_models_list]
    model_FDR_list = [
        [FDR_metric[x][m] for x in FDR_metric if (m in FDR_metric[x].index) and (not pd.isna(FDR_metric[x][m]))] for m
        in all_models_list]
    model_AUROC_list = [
        [AUROC_metric[x][m] for x in AUROC_metric if (m in AUROC_metric[x].index) and (not pd.isna(AUROC_metric[x][m]))]
        for m
        in all_models_list]
    model_AUPRC_list = [
        [AUPRC_metric[x][m] for x in AUPRC_metric if (m in AUPRC_metric[x].index) and (not pd.isna(AUPRC_metric[x][m]))]
        for m
        in all_models_list]
    # -----
    print(all_models_list)
    print("Avg F1: ", [np.mean(each) for each in model_F1_list])
    print("Std F1: ", [np.std(each) for each in model_F1_list])
    print("Avg FDR: ", [np.mean(each) for each in model_FDR_list])
    print("Std FDR: ", [np.std(each) for each in model_FDR_list])
    print("Avg AUROC: ", [np.mean(each) for each in model_AUROC_list])
    print("Std AUROC: ", [np.std(each) for each in model_AUROC_list])
    print("Avg AUPRC: ", [np.mean(each) for each in model_AUPRC_list])
    print("Std AUPRC: ", [np.std(each) for each in model_AUPRC_list])
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_F1_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1 score")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.8, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03])
    plt.ylim(0.0, 1.0)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] + 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] + 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA10/ZIPoisson-F1_score.pdf")
    plt.show()
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_AUROC_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("AUROC")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=1.0, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03])
    plt.ylim(0.0, 1.0)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] + 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] + 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA10/ZIPoisson-AUROC.pdf")
    plt.show()
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_AUPRC_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("AUPRC")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.8, linestyles="solid", colors=gray_color)
    # plt.yticks([0.0, 0.01, 0.02, 0.03], [0.0, 0.01, 0.02, 0.03])
    plt.ylim(0.0, 1.0)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] + 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] + 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA10/ZIPoisson-AUPRC.pdf")
    plt.show()
    # -----
    model_avg_FDR = np.asarray([np.mean(each) for each in model_FDR_list])
    model_std_FDR = np.asarray([np.std(each) for each in model_FDR_list])
    plt.figure(figsize=(8, 5))
    width = 0.6
    plt.bar(np.arange(len(all_models_list)), model_avg_FDR, label="FDR", width=width, color=BlueRed_12.mpl_colors[0])
    plt.bar(np.arange(len(all_models_list)), 1 - model_avg_FDR, label="PPV", width=width, bottom=model_avg_FDR,
            color=BlueRed_12.mpl_colors[-2])
    plt.xticks(np.arange(len(all_models_list)), [model_name_dict[x] for x in all_models_list], fontsize=13)
    # plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    # plt.legend(bbox_to_anchor=(1.0, 0.9))
    plt.legend(bbox_to_anchor=(0.1, 1.0), ncol=2)
    plt.tight_layout()
    plt.savefig("../fig/FigA10/ZIPoisson-FDR.pdf")
    plt.show()

# ========================================================
# Fig A19: method performances on Tabula Muris data
# ========================================================
def performanceTabulaMuris():
    print("Model Performance Compared with ChIP-seq Network ")
    metric = np.load("../fig_res/tabula_TF_PCC-cortex-metric.npy", allow_pickle=True).item()
    metric = {each: metric[each] for each in metric if len(metric[each]) > 0}
    # -----
    all_models_list = ["pearson", "spearman", "glasso", "GENIE3", "minet", "PIDC", "scDesign2", "scLink"]
    model_TPR_list = [[metric[x][m]["f1"]["TPR"] for x in metric] for m in all_models_list]
    model_F1_list = [[metric[x][m]["f1"]["f1"] for x in metric] for m in all_models_list]
    model_precision_list = [[metric[x][m]["f1"]["precision"] for x in metric] for m in all_models_list]
    model_FDR_list = [[1-metric[x][m]["f1"]["precision"] for x in metric] for m in all_models_list]
    # -----
    print(all_models_list)
    print("Avg TPR: ", [np.mean(each) for each in model_TPR_list])
    print("Std TPR: ", [np.std(each) for each in model_TPR_list])
    print("Avg F1: ", [np.mean(each) for each in model_F1_list])
    print("Std F1: ", [np.std(each) for each in model_F1_list])
    print("Avg Precision: ", [np.mean(each) for each in model_precision_list])
    print("Std Precision: ", [np.std(each) for each in model_precision_list])
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_F1_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.9, linestyles="solid", colors=gray_color)
    # plt.yticks([0.3, 0.5, 0.7, 0.9], [0.3, 0.5, 0.7, 0.9])
    plt.ylim(0.0, 0.6)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=0.61, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=0.61, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA19/tabula_TF_net-F1.pdf")
    plt.show()
    # -----
    plt.figure(figsize=(8, 5))
    bplt = plt.boxplot(model_TPR_list, patch_artist=True)
    plt.xticks(np.arange(1, len(all_models_list) + 1), [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("TPR")
    color_list = [model_color_dict[each] for each in all_models_list]
    for i, patch in enumerate(bplt["boxes"]):
        patch.set_facecolor(color_list[i])
        patch.set_linewidth(2.0)
    for patch in bplt["fliers"]:
        patch.set_markeredgecolor("k")
        patch.set_markeredgewidth(1.5)
        patch.set_markersize(10)
    for patch in bplt['medians']:
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch in bplt['whiskers']:
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch in bplt['caps']:
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.vlines(x=[5.5], ymin=0.0, ymax=1.0, linestyles="solid", colors=gray_color)
    # plt.yticks([0.3, 0.5, 0.7, 0.9], [0.3, 0.5, 0.7, 0.9])
    plt.ylim(0.0, 1.0)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=1.1, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=1.1, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/FigA19/tabula_TF_net-TPR.pdf")
    plt.show()
    # -----
    model_avg_FDR = np.asarray([np.mean(each) for each in model_FDR_list])
    model_std_FDR = np.asarray([np.std(each) for each in model_FDR_list])
    plt.figure(figsize=(8, 5))
    width = 0.6
    plt.bar(np.arange(len(all_models_list)), model_avg_FDR, label="FDR", width=width, color=BlueRed_12.mpl_colors[0])
    plt.bar(np.arange(len(all_models_list)), 1 - model_avg_FDR, label="PPV", width=width, bottom=model_avg_FDR,
            color=BlueRed_12.mpl_colors[-2])
    plt.xticks(np.arange(len(all_models_list)), [model_name_dict[x] for x in all_models_list], fontsize=13)
    # plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    # plt.legend(bbox_to_anchor=(1.0, 0.9))
    plt.legend(bbox_to_anchor=(0.1, 1.0), ncol=2)
    plt.tight_layout()
    plt.savefig("../fig/FigA19/tabula_TF_net-FDR.pdf")
    plt.show()
    # -----


if __name__ == '__main__':
    # -----
    # Fig A7: NORTA simulations w/ 500 genes
    NORTA500Genes()
    # -----
    # Fig A8: SERGIO simulations w/ 400 genes
    SERGIO400Genes()
    # -----
    # Fig A9: ZI-Gaussian simulations
    performanceZIGaussian()
    # -----
    # Fig A10: ZI-Poisson simulations
    performanceZIPoisson()
    # -----
    # Fig A19: Tabula-Muris data
    performanceTabulaMuris()

