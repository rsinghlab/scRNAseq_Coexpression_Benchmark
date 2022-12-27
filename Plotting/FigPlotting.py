'''
Description:
    Paper figures plotting.

Author:
    Jiaqi Zhang <jiaqi_zhang2@brown.edu>
'''
import pandas as pd
import numpy as np
import itertools
import scipy.stats
import scipy
from PlottingUtils import *



# ========================================================
#   Fig 4: experimental data v.s. ZI-Gaussian/ZI-Poisson
# ========================================================

def compareOldSimCellSparsity():
    '''
    Compare cell sparsity between experimental data and old simulations of 100 HVGs..
    '''
    print("=" * 70)
    print("Experimental vs. Old Simulation | Cell Sparsity")
    exp_stats_dict = np.load("../fig_res/experimental-basic_stats.npy", allow_pickle=True).item()["100hvg"]
    sim_stats_dict = np.load("../fig_res/old_sim-basic_stats.npy", allow_pickle=True).item()["100hvg"]
    exp_data_name = (
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"])] +
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"])]
    )
    # -----
    # P-values
    print("-" * 70)
    pval_mat = []
    for each_exp in exp_data_name:
        exp_tmp_sparsity = exp_stats_dict[each_exp]["cell_zero_percent"]
        tmp_pval = []
        for each_sim in sim_stats_dict:
            sim_tmp_sparsity = sim_stats_dict[each_sim]["cell_zero_percent"]
            pval = scipy.stats.ks_2samp(exp_tmp_sparsity, sim_tmp_sparsity).pvalue
            tmp_pval.append(pval)
        pval_mat.append(tmp_pval)
    pval_df = pd.DataFrame(data=pval_mat, index=exp_data_name, columns=list(sim_stats_dict.keys()))
    pval_mat = np.asarray(pval_mat)
    for j in range(len(sim_stats_dict)):
        if np.all(pval_mat[:, j] < 1e-5):
            print("{}: ***".format(list(sim_stats_dict.keys())[j]))
    # Mean (std)
    print("-" * 70)
    exp_sparsity = []
    for each_exp in exp_data_name:
        exp_tmp_sparsity = exp_stats_dict[each_exp]["cell_zero_percent"]
        exp_sparsity.append(np.nanmean(exp_tmp_sparsity))
        print("{}: {}({})".format(each_exp, np.nanmean(exp_tmp_sparsity), np.nanstd(exp_tmp_sparsity)))
    print("***** Exp Data *****: {}".format(np.mean(exp_sparsity)))
    for each_sim in sim_stats_dict:
        sim_tmp_sparsity = sim_stats_dict[each_sim]["cell_zero_percent"]
        print("{}: {}({})".format(each_sim, np.nanmean(sim_tmp_sparsity), np.nanstd(sim_tmp_sparsity)))
    # -----
    # Plotting
    x_ticks = np.arange(len(exp_stats_dict) + len(sim_stats_dict)) + 1
    x_ticks_label = ["#{}".format(i + 1) for i in range(len(exp_stats_dict))] + \
                    ["{}".format(each.split("_")[-1]) for each in sim_stats_dict if "scLink" in each] + \
                    ["{}".format(each.split("_")[-1]) for each in sim_stats_dict if "ZILGM" in each]
    # x_ticks_label = ["#{}".format(i + 1) for i in range(len(exp_stats_dict))] + \
    #                 ["$\\rho=${}".format(each.split("_")[-1]) for each in sim_stats_dict if "scLink" in each] + \
    #                 ["$\pi={}$".format(each.split("_")[-1]) for each in sim_stats_dict if "ZILGM" in each]
    plt.figure(figsize=(10, 5))
    bplt = plt.boxplot(
        [exp_stats_dict[each]["cell_zero_percent"] for each in exp_data_name] +
        [sim_stats_dict[each]["cell_zero_percent"] for each in sim_stats_dict],
        patch_artist=True
    )
    colors = [simulation_color_dict["ref"] for _ in range(len(exp_stats_dict))] + [
        simulation_color_dict[each.split("_")[0]] for each in sim_stats_dict]
    for patch, color in zip(bplt['boxes'], colors):
        patch.set_facecolor(color)
    for patch, color in zip(bplt['medians'], colors):
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch, color in zip(bplt['whiskers'], colors):
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch, color in zip(bplt['caps'], colors):
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.xticks(x_ticks, x_ticks_label, fontsize=18)
    plt.xlabel("Dataset")
    plt.ylabel("Percent of Non-Zeros per Cell", fontsize=20)
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=20)
    plt.ylim(0.0, 1.1)
    plt.vlines(x=[8.5, 12.5], ymin=0.0, ymax=1.1, linestyles="solid", colors=gray_color)
    plt.text(x=4.5, y=plt.gca().get_ylim()[1] - 0.01, s="Experimental \n", color="k", fontsize=20, horizontalalignment="center")
    plt.text(x=10.5, y=plt.gca().get_ylim()[1] - 0.01, s="ZI-Gaussian \n (varying $\\rho$)", color=simulation_color_dict["scLink"], fontsize=20, horizontalalignment="center")
    plt.text(x=14, y=plt.gca().get_ylim()[1] - 0.01, s="ZI-Poisson \n (varying $\pi$)", color=simulation_color_dict["ZILGM"], fontsize=20, horizontalalignment="center")
    # plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig4/cell_sparsity.pdf")
    plt.show()


def compareOldSimGeneSparsity():
    '''
    Compare gene sparsity between experimental data and old simulations of 100 HVGs..
    '''
    print("=" * 70)
    print("Experimental vs. Old Simulation | Gene Sparsity")
    exp_stats_dict = np.load("../fig_res/experimental-basic_stats.npy", allow_pickle=True).item()["100hvg"]
    sim_stats_dict = np.load("../fig_res/old_sim-basic_stats.npy", allow_pickle=True).item()["100hvg"]
    exp_data_name = (
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"])] +
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"])]
    )
    # -----
    # P-values
    print("-" * 70)
    pval_mat = []
    for each_exp in exp_data_name:
        exp_tmp_sparsity = exp_stats_dict[each_exp]["gene_zero_percent"]
        tmp_pval = []
        for each_sim in sim_stats_dict:
            sim_tmp_sparsity = sim_stats_dict[each_sim]["gene_zero_percent"]
            pval = scipy.stats.ks_2samp(exp_tmp_sparsity, sim_tmp_sparsity).pvalue
            tmp_pval.append(pval)
        pval_mat.append(tmp_pval)
    pval_df = pd.DataFrame(data=pval_mat, index=exp_data_name, columns=list(sim_stats_dict.keys()))
    pval_mat = np.asarray(pval_mat)
    for j in range(len(sim_stats_dict)):
        if np.all(pval_mat[:, j] < 1e-5):
            print("{}: ***".format(list(sim_stats_dict.keys())[j]))
    # Mean (std)
    print("-" * 70)
    exp_sparsity = []
    for each_exp in exp_data_name:
        exp_tmp_sparsity = exp_stats_dict[each_exp]["gene_zero_percent"]
        exp_sparsity.append(np.nanmean(exp_tmp_sparsity))
        print("{}: {}({})".format(each_exp, np.nanmean(exp_tmp_sparsity), np.nanstd(exp_tmp_sparsity)))
    print("***** Exp Data *****: {}".format(np.mean(np.asarray(exp_sparsity))))
    for each_sim in sim_stats_dict:
        sim_tmp_sparsity = sim_stats_dict[each_sim]["gene_zero_percent"]
        print("{}: {}({})".format(each_sim, np.nanmean(sim_tmp_sparsity), np.nanstd(sim_tmp_sparsity)))
    # -----
    # Plotting
    x_ticks = np.arange(len(exp_stats_dict) + len(sim_stats_dict)) + 1
    # x_ticks_label = ["#{}".format(i + 1) for i in range(len(exp_stats_dict))] + \
    #                 ["$\\rho=${}".format(each.split("_")[-1]) for each in sim_stats_dict if "scLink" in each] + \
    #                 ["$\pi={}$".format(each.split("_")[-1]) for each in sim_stats_dict if "ZILGM" in each]
    x_ticks_label = ["#{}".format(i + 1) for i in range(len(exp_stats_dict))] + \
                    ["{}".format(each.split("_")[-1]) for each in sim_stats_dict if "scLink" in each] + \
                    ["{}".format(each.split("_")[-1]) for each in sim_stats_dict if "ZILGM" in each]
    plt.figure(figsize=(10, 5))
    bplt = plt.boxplot(
        [exp_stats_dict[each]["gene_zero_percent"] for each in exp_data_name] +
        [sim_stats_dict[each]["gene_zero_percent"] for each in sim_stats_dict],
        patch_artist=True
    )
    colors = [simulation_color_dict["ref"] for _ in range(len(exp_stats_dict))] + [
        simulation_color_dict[each.split("_")[0]] for each in sim_stats_dict]
    for patch, color in zip(bplt['boxes'], colors):
        patch.set_facecolor(color)
    for patch, color in zip(bplt['medians'], colors):
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch, color in zip(bplt['whiskers'], colors):
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch, color in zip(bplt['caps'], colors):
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.xticks(x_ticks, x_ticks_label, fontsize=18)
    plt.xlabel("Dataset")
    plt.ylabel("Percent of Non-Zeros per Gene", fontsize=20)
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=20)
    plt.ylim(0.0, 1.1)
    plt.vlines(x=[8.5, 12.5], ymin=0.0, ymax=1.1, linestyles="solid", colors=gray_color)
    plt.text(x=4.5, y=plt.gca().get_ylim()[1] - 0.01, s="Experimental \n", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=10.5, y=plt.gca().get_ylim()[1] - 0.01, s="ZI-Gaussian \n (varying $\\rho$)", color=simulation_color_dict["scLink"],
             fontsize=20, horizontalalignment="center")
    plt.text(x=14, y=plt.gca().get_ylim()[1] - 0.01, s="ZI-Poisson \n (varying $\pi$)", color=simulation_color_dict["ZILGM"], fontsize=20,
             horizontalalignment="center")
    # plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig4/gene_sparsity.pdf")
    plt.show()


def compareOldSimDistShape():
    '''
    Compare distribution shape between scaled experimental data and old simulations of 100 HVGs..
    '''
    print("=" * 70)
    print("Experimental vs. Old Simulation | Distribution Shape")
    exp_hist_dict = np.load("../fig_res/experimental-hist.npy", allow_pickle=True).item()["100hvg"]
    sim_hist_dict = np.load("../fig_res/old_sim-hist.npy", allow_pickle=True).item()["100hvg"]
    exp_data_name = (
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"])] +
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"])]
    )
    # -----
    pval_list = []  # two sample ks test p-value
    ise_list = []  # integrated squared error
    # Plotting
    # =======================
    # DEPRECATED
    # plt.figure(figsize=(10, 10))
    # x_ticks = np.linspace(0.0, 1.0, 20)
    # xtick_label = ["{:.2f}".format(0.0 + 0.05 * i) for i in range(len(x_ticks))]
    # width = (x_ticks[1] - x_ticks[0]) / 3.0
    # for i, each_e in enumerate(exp_data_name[:4]):
    #     plt.subplot(2, 4, i + 1)
    #     exp_data_hist = exp_hist_dict[each_e]
    #     exp_data_freq = exp_data_hist / np.sum(exp_data_hist)
    #     scLink_data_hist = sim_hist_dict["scLink_0.07"]
    #     scLink_data_freq = scLink_data_hist / np.sum(scLink_data_hist)
    #     zilgm_data_hist = sim_hist_dict["ZILGM_0.8"]
    #     zilgm_data_freq = zilgm_data_hist / np.sum(zilgm_data_hist)
    #     # -----
    #     pval_list.append([scipy.stats.ks_2samp(exp_data_freq, scLink_data_freq).pvalue,
    #                       scipy.stats.ks_2samp(exp_data_freq, zilgm_data_freq).pvalue])
    #     ise_list.append([ISE(exp_data_freq, scLink_data_freq), ISE(exp_data_freq, zilgm_data_freq)])
    #     # -----
    #     plt.bar(x_ticks - 3 * width, exp_data_freq, width=width, label="Experimental",
    #             color=simulation_color_dict["ref"], align="edge", linewidth=1.0, edgecolor="k")
    #     plt.bar(x_ticks - 2 * width, scLink_data_freq, width=width, label="ZI-Gaussian",
    #             color=simulation_color_dict["scLink"], align="edge", linewidth=1.0, edgecolor="k")
    #     plt.bar(x_ticks - 1 * width, zilgm_data_freq, width=width, label="ZI-Poisson",
    #             color=simulation_color_dict["ZILGM"], align="edge", linewidth=1.0, edgecolor="k")
    #     plt.xticks([], [])
    #     plt.yticks([], [])
    #     plt.gca().spines['right'].set_visible(False)
    #     plt.gca().spines['left'].set_visible(False)
    #     plt.gca().spines['top'].set_visible(False)
    #     plt.gca().spines['bottom'].set_visible(False)
    #     if i == 3:
    #         plt.legend(ncol=1)
    #     # -----
    #     plt.subplot(2, 4, i + 5)
    #     plt.bar(x_ticks[1:] - 3 * width, exp_data_freq[1:], width=width, label="Experimental",
    #             color=simulation_color_dict["ref"], align="edge", linewidth=1.0, edgecolor="k")
    #     plt.bar(x_ticks[1:] - 2 * width, scLink_data_freq[1:], width=width, label="ZI-Gaussian",
    #             color=simulation_color_dict["scLink"], align="edge", linewidth=1.0, edgecolor="k")
    #     plt.bar(x_ticks[1:] - 1 * width, zilgm_data_freq[1:], width=width, label="ZI-Poisson",
    #             color=simulation_color_dict["ZILGM"], align="edge", linewidth=1.0, edgecolor="k")
    #     plt.xticks([], [])
    #     plt.yticks([], [])
    #     plt.gca().spines['right'].set_visible(False)
    #     plt.gca().spines['left'].set_visible(False)
    #     plt.gca().spines['top'].set_visible(False)
    #     plt.gca().spines['bottom'].set_visible(False)
    # plt.tight_layout()
    # plt.show()
    # =======================
    fig = plt.figure(figsize=(10, 10))
    x_ticks = np.linspace(0.0, 1.0, 20)
    xtick_label = ["{:.2f}".format(0.0 + 0.05 * i) for i in range(len(x_ticks))]
    width = (x_ticks[1] - x_ticks[0])
    i_dict = {0:2, 1:3, 2:0, 3:1}
    for i, each_e in enumerate(exp_data_name[:4]):
        exp_data_hist = exp_hist_dict[exp_data_name[i_dict[i]]]
        exp_data_freq = exp_data_hist / np.sum(exp_data_hist)
        scLink_data_hist = sim_hist_dict["scLink_0.07"]
        scLink_data_freq = scLink_data_hist / np.sum(scLink_data_hist)
        zilgm_data_hist = sim_hist_dict["ZILGM_0.8"]
        zilgm_data_freq = zilgm_data_hist / np.sum(zilgm_data_hist)
        # -----
        pval_list.append([scipy.stats.ks_2samp(exp_data_freq, scLink_data_freq).pvalue,
                          scipy.stats.ks_2samp(exp_data_freq, zilgm_data_freq).pvalue])
        ise_list.append([ISE(exp_data_freq, scLink_data_freq), ISE(exp_data_freq, zilgm_data_freq)])
        # -----
        outergs = gridspec.GridSpec(1, 1)
        outergs.update(bottom=(i // 2) * .49 + 0.01, left=(i % 2) * .5 + 0.01,
                       top=(1 + i // 2) * .49 - 0.01, right=(1 + i % 2) * .5 - 0.01)
        outerax = fig.add_subplot(outergs[0])
        outerax.tick_params(axis='both', which='both', bottom=0, left=0, labelbottom=0, labelleft=0)
        # Add inner grid
        gs = gridspec.GridSpec(1, 3)
        gs.update(bottom=(i // 2) * .49 + 0.05, left=(i % 2) * .5 + 0.05,
                  top=(1 + i // 2) * .49 - 0.04, right=(1 + i % 2) * .5 - 0.05,
                  wspace=0.0, hspace=0.0)
        # Histograms
        ax = fig.add_subplot(gs[0])
        ax.set_xlabel("Exp.", fontsize=18)
        ax.bar(x_ticks[1:], exp_data_freq[1:], width=width, color=simulation_color_dict["ref"], align="edge", linewidth=1.0, edgecolor="k")
        removeAllBorders()
        plt.xticks([], [])
        plt.yticks([], [])
        ax = fig.add_subplot(gs[1])
        ax.set_title("#{}".format(i_dict[i] + 1), fontsize=18, pad=0.0)
        ax.set_xlabel("ZI-Gaussian", fontsize=18)
        ax.bar(x_ticks[1:], scLink_data_freq[1:], width=width, color=simulation_color_dict["scLink"], align="edge", linewidth=1.0, edgecolor="k")
        removeAllBorders()
        plt.xticks([], [])
        plt.yticks([], [])
        ax = fig.add_subplot(gs[2])
        ax.set_xlabel("ZI-Poisson", fontsize=18)
        ax.bar(x_ticks[1:], zilgm_data_freq[1:], width=width, color=simulation_color_dict["ZILGM"], align="edge", linewidth=1.0, edgecolor="k")
        removeAllBorders()
        plt.xticks([], [])
        plt.yticks([], [])
    plt.tight_layout()
    plt.savefig("../fig/Fig4/distribution_shape.pdf")
    plt.show()
    # -----
    pval_df = pd.DataFrame(data=pval_list, index=exp_data_name[:4], columns=["ZI-Gaussian", "ZILGM"])
    ise_df = pd.DataFrame(data=ise_list, index=exp_data_name[:4], columns=["ZI-Gaussian", "ZILGM"])
    print(pval_df)
    print(ise_df)
    # -----
    print("-" * 70)
    iae_list = []
    for i, each_e in enumerate(exp_data_name):
        exp_data_hist = exp_hist_dict[each_e]
        exp_data_freq = exp_data_hist / np.sum(exp_data_hist)
        for rho in ["0.07", "0.10", "0.13", "0.16"]:
            scLink_data_hist = sim_hist_dict["scLink_{}".format(rho)]
            scLink_data_freq = scLink_data_hist / np.sum(scLink_data_hist)
            iae_list.append(IAE(exp_data_freq, scLink_data_freq))
    print("Avg IAE (scLink): {}".format(np.mean(iae_list)))
    iae_list = []
    for i, each_e in enumerate(exp_data_name):
        exp_data_hist = exp_hist_dict[each_e]
        exp_data_freq = exp_data_hist / np.sum(exp_data_hist)
        for rho in ["0.8", "0.9", "1.0"]:
            zilgm_data_hist = sim_hist_dict["ZILGM_{}".format(rho)]
            zilgm_data_freq = zilgm_data_hist / np.sum(zilgm_data_hist)
            iae_list.append(IAE(exp_data_freq, zilgm_data_freq))
    print("Avg IAE (ZILGM): {}".format(np.mean(iae_list)))
    iae_list = []
    for i, each_e in enumerate(exp_data_name):
        exp_data_hist = exp_hist_dict[each_e]
        exp_data_freq = exp_data_hist / np.sum(exp_data_hist)
        for j, another_e in enumerate(exp_data_name):
            if i == j:
                continue
            another_data_hist = exp_hist_dict[another_e]
            another_data_freq = another_data_hist / np.sum(another_data_hist)
            iae_list.append(IAE(exp_data_freq, another_data_freq))
    print("Avg IAE (exp): {}".format(np.mean(iae_list)))

# ========================================================
#    Fig 5: experimental data v.s. NORTA/SERGIO
# ========================================================

def compareNORTAGeneSparsity():
    '''
    Compare gene sparsity between experimental data and NORTA simulations of 100 HVGs..
    '''
    print("=" * 70)
    print("Experimental vs. NORTA Simulation | Gene Sparsity")
    exp_stats_dict = np.load("../fig_res/experimental-basic_stats.npy", allow_pickle=True).item()["100hvg"]
    sim_stats_dict = np.load("../fig_res/new_sim-basic_stats.npy", allow_pickle=True).item()["100hvg"]
    exp_data_name = (
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"])] +
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"])]
    )
    # -----
    # P-values
    print("-" * 70)
    pval_mat = []
    for each_exp in exp_data_name:
        exp_tmp_sparsity = exp_stats_dict[each_exp]["gene_zero_percent"]
        sim_tmp_sparsity = sim_stats_dict["NORTA-{}".format(each_exp)]["gene_zero_percent"]
        pval = scipy.stats.ks_2samp(exp_tmp_sparsity, sim_tmp_sparsity).pvalue
        pval_mat.append(pval)
    pval_mat = np.asarray(pval_mat)
    print("p-value: ", pval_mat)
    # Mean (std)
    print("-" * 70)
    for each_exp in exp_data_name:
        exp_tmp_sparsity = exp_stats_dict[each_exp]["gene_zero_percent"]
        print("{}: {}({})".format(each_exp, np.nanmean(exp_tmp_sparsity), np.nanstd(exp_tmp_sparsity)))
    sim_sparsity = []
    for each_sim in sim_stats_dict:
        sim_tmp_sparsity = sim_stats_dict[each_sim]["gene_zero_percent"]
        sim_sparsity.append(np.nanmean(sim_tmp_sparsity))
        print("{}: {}({})".format(each_sim, np.nanmean(sim_tmp_sparsity), np.nanstd(sim_tmp_sparsity)))
    print("*****b Avg NORTA *****: {}".format(np.mean(sim_sparsity)))
    # -----
    # Plotting
    x_ticks = np.arange(len(exp_data_name)) + 1
    x_ticks_label = ["#{}".format(i + 1) for i in range(len(exp_data_name))]
    width = 0.2
    plt.figure(figsize=(10, 5))
    bplt1 = plt.boxplot([exp_stats_dict[each]["gene_zero_percent"] for each in exp_data_name], positions=x_ticks - 0.15,
                        widths=width, patch_artist=True)
    bplt2 = plt.boxplot([sim_stats_dict["NORTA-{}".format(each)]["gene_zero_percent"] for each in exp_data_name],
                        positions=x_ticks + 0.15, widths=width, patch_artist=True)
    plt.xticks(x_ticks, x_ticks_label)
    plt.xlabel("Dataset")
    colors = [simulation_color_dict["ref"], simulation_color_dict["NORTA"]]
    for i, bplt in enumerate((bplt1, bplt2)):
        for patch in bplt["boxes"]:
            patch.set_facecolor(colors[i])
            patch.set_linewidth(1.5)
        for patch in bplt["fliers"]:
            patch.set_markeredgecolor(colors[i])
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
    plt.ylabel("Percent of Non-Zeros per Gene", fontsize=20)
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=20)
    plt.ylim(0.0, 1.1)
    plt.gca().spines['top'].set_visible(False)
    plt.bar([1], [np.nan], color=colors[0], label="Experimental")
    plt.bar([1], [np.nan], color=colors[1], label="NORTA")
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.savefig("../fig/Fig5/gene_sparsity_NORTA.pdf")
    plt.show()


def compareNORTACorr():
    '''
    Compare CMD between experimental data and NORTA simulations of 100 HVGs..
    '''
    print("=" * 70)
    print("Experimental vs. NORTA Simulation | CMD")
    exp_data_name = (
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"])] +
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"])]
    )
    exp_corr_dict = np.load("../fig_res/experimental-corr.npy", allow_pickle=True).item()["100hvg"]
    sim_corr_dict = np.load("../fig_res/new_sim-corr.npy", allow_pickle=True).item()["100hvg"]
    net_mat = {each: loadNetMat("../data/simulation/NORTA/{}-net_mat.csv".format(each)) for each in exp_data_name}
    # -----
    pcc_CMD = []
    pcc_MSE = []
    for each in exp_data_name:
        CD, diff_MSE, (before_CMD, after_CMD), (before_MSE, after_MSE) = computeCorrSim(
            exp_corr_dict[each], sim_corr_dict["NORTA-{}".format(each)], net_mat[each].values, type="PCC"
        )
        pcc_CMD.append((before_CMD, after_CMD))
        pcc_MSE.append((before_MSE, after_MSE))
    # -----
    before_cmd_list = [each[0] for each in pcc_CMD]
    after_cmd_list = [each[1] for each in pcc_CMD]
    before_mse_list = [each[0] for each in pcc_MSE]
    after_mse_list = [each[1] for each in pcc_MSE]
    data_list = ["#{}".format(i + 1) for i in range(len(exp_data_name))]
    print("Before CMD: ", before_cmd_list)
    print("After CMD: ", after_cmd_list)
    print("Before MSE: ", before_mse_list)
    print("After MSE: ", after_mse_list)
    # -----
    width = 0.4
    plt.figure(figsize=(10, 5))
    plt.bar(np.arange(len(data_list)) - width, before_cmd_list, width=width, label="Experimental", align="edge",
            color=simulation_color_dict["ref"], edgecolor="k", linewidth=2)
    plt.bar(np.arange(len(data_list)), after_cmd_list, width=width, label="NORTA", align="edge",
            color=simulation_color_dict["NORTA"], edgecolor="k", linewidth=2)
    plt.xticks(np.arange(len(data_list)), data_list)
    plt.xlabel("Dataset")
    plt.ylabel("Correlation Matrix Distance")
    plt.legend(ncol=2)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig5/CMD_NORTA.pdf")
    plt.show()
    # -----
    plt.figure(figsize=(10, 5))
    plt.bar(np.arange(len(data_list)) - width, before_mse_list, width=width, label="Experimental", align="edge",
            color=simulation_color_dict["ref"], edgecolor="k", linewidth=2)
    plt.bar(np.arange(len(data_list)), after_mse_list, width=width, label="NORTA", align="edge",
            color=simulation_color_dict["NORTA"], edgecolor="k", linewidth=2)
    plt.xticks(np.arange(len(data_list)), data_list)
    plt.xlabel("Dataset")
    plt.ylabel("Mean Squared Error")
    plt.legend(ncol=2)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig5/MSE_NORTA.pdf")
    plt.show()


def compareSERGIOGeneSparsity():
    '''
    Compare gene sparsity between experimental data and SERGIO simulations of 100 HVGs..
    '''
    print("=" * 70)
    print("Experimental vs. SERGIO Simulation | Gene Sparsity")
    exp_stats_dict = np.load("../fig_res/experimental-basic_stats.npy", allow_pickle=True).item()["100hvg"]
    sim_stats_dict = np.load("../fig_res/new_sim-basic_stats.npy", allow_pickle=True).item()["100hvg"]
    sim_stats_dict = {each: sim_stats_dict[each] for each in sim_stats_dict if "SERGIO" in each}
    exp_data_name = (
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"])] +
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"])]
    )
    # -----
    # P-values
    print("-" * 70)
    pval_mat = []
    for each_exp in exp_data_name:
        exp_tmp_sparsity = exp_stats_dict[each_exp]["gene_zero_percent"]
        tmp_pval = []
        for each_sim in sim_stats_dict:
            sim_tmp_sparsity = sim_stats_dict[each_sim]["gene_zero_percent"]
            pval = scipy.stats.ks_2samp(exp_tmp_sparsity, sim_tmp_sparsity).pvalue
            tmp_pval.append(pval)
        pval_mat.append(tmp_pval)
    pval_df = pd.DataFrame(data=pval_mat, index=exp_data_name, columns=list(sim_stats_dict.keys()))
    pval_mat = np.asarray(pval_mat)
    for j in range(len(sim_stats_dict)):
        if np.all(pval_mat[:, j] < 1e-5):
            print("{}: ***".format(list(sim_stats_dict.keys())[j]))
        else:
            print("{}: {}".format(list(sim_stats_dict.keys())[j], np.max(pval_mat[:, j])))
    # Mean (std)
    print("-" * 70)
    for each_exp in exp_data_name:
        exp_tmp_sparsity = exp_stats_dict[each_exp]["gene_zero_percent"]
        print("{}: {}({})".format(each_exp, np.nanmean(exp_tmp_sparsity), np.nanstd(exp_tmp_sparsity)))
    sim_sparsity = []
    for each_sim in sim_stats_dict:
        sim_tmp_sparsity = sim_stats_dict[each_sim]["gene_zero_percent"]
        sim_sparsity.append(np.nanmean(sim_tmp_sparsity))
        print("{}: {}({})".format(each_sim, np.nanmean(sim_tmp_sparsity), np.nanstd(sim_tmp_sparsity)))
    print("***** Avg SERGIO *****: {}".format(np.mean(sim_sparsity)))
    # -----
    # Plotting
    x_ticks = np.arange(len(exp_stats_dict) + len(sim_stats_dict)) + 1
    # x_ticks_label = ["#{}".format(i + 1) for i in range(len(exp_stats_dict))] + \
    #                 ["$s=${}".format(each.split("-")[1]) for each in sim_stats_dict]
    x_ticks_label = ["#{}".format(i + 1) for i in range(len(exp_stats_dict))] + \
                    ["{}".format(each.split("-")[1]) for each in sim_stats_dict]
    plt.figure(figsize=(10, 5))
    bplt = plt.boxplot(
        [exp_stats_dict[each]["gene_zero_percent"] for each in exp_data_name] +
        [sim_stats_dict[each]["gene_zero_percent"] for each in sim_stats_dict],
        patch_artist=True
    )
    colors = [simulation_color_dict["ref"] for _ in range(len(exp_stats_dict))] + [
        simulation_color_dict[each.split("-")[0]] for each in sim_stats_dict]
    for patch, color in zip(bplt['boxes'], colors):
        patch.set_facecolor(color)
    for patch, color in zip(bplt['medians'], colors):
        patch.set_color("k")
        patch.set_linewidth(1.0)
        patch.set_linestyle("--")
    for patch, color in zip(bplt['whiskers'], colors):
        patch.set_color("k")
        patch.set_linewidth(2)
    for patch, color in zip(bplt['caps'], colors):
        patch.set_color("k")
        patch.set_linewidth(2)
    plt.xticks(x_ticks, x_ticks_label, fontsize=18)
    plt.xlabel("Dataset")
    plt.ylabel("Percent of Non-Zeros per Gene", fontsize=20)
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=20)
    plt.ylim(0.0, 1.1)
    plt.vlines(x=[8.5], ymin=0.0, ymax=1.1, linestyles="solid", colors=gray_color)
    # plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=4.5, y=plt.gca().get_ylim()[1] - 0.01, s="Experimental \n", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=11, y=plt.gca().get_ylim()[1] - 0.01, s="SERGIO \n (varying s)", color=simulation_color_dict["SERGIO"],
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/Fig5/gene_sparsity_SERGIO.pdf")
    plt.show()


def compareNORTADistShape():
    '''
    Compare distribution shape between scaled experimental data and NORTA simulations of 100 HVGs..
    '''
    print("=" * 70)
    print("Experimental vs. NORTA Simulation | Distribution Shape")
    exp_hist_dict = np.load("../fig_res/experimental-hist.npy", allow_pickle=True).item()["100hvg"]
    sim_hist_dict = np.load("../fig_res/new_sim-hist.npy", allow_pickle=True).item()["100hvg"]
    exp_data_name = (
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"])] +
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"])]
    )
    # -----
    pval_list = []  # two sample ks test p-value
    ise_list = []  # integrated squared error
    # Plotting
    fig = plt.figure(figsize=(10, 8))
    x_ticks = np.linspace(0.0, 1.0, 20)
    xtick_label = ["{:.2f}".format(0.0 + 0.05 * i) for i in range(len(x_ticks))]
    width = (x_ticks[1] - x_ticks[0])
    i_dict = {0: 2, 1: 3, 2: 0, 3: 1}
    for i, each_e in enumerate(exp_data_name[:4]):
        exp_data_hist = exp_hist_dict[exp_data_name[i_dict[i]]]
        exp_data_freq = exp_data_hist / np.sum(exp_data_hist)
        norta_data_hist = sim_hist_dict["NORTA-{}".format(each_e)]
        norta_data_freq = norta_data_hist / np.sum(norta_data_hist)
        # -----
        pval_list.append(scipy.stats.ks_2samp(exp_data_freq, norta_data_freq).pvalue)
        ise_list.append(ISE(exp_data_freq, norta_data_freq))
        # -----
        outergs = gridspec.GridSpec(1, 1)
        outergs.update(bottom=(i // 2) * .49 + 0.01, left=(i % 2) * .5 + 0.01,
                       top=(1 + i // 2) * .49 - 0.01, right=(1 + i % 2) * .5 - 0.01)
        outerax = fig.add_subplot(outergs[0])
        outerax.tick_params(axis='both', which='both', bottom=0, left=0, labelbottom=0, labelleft=0)
        # Add inner grid
        gs = gridspec.GridSpec(1, 2)
        gs.update(bottom=(i // 2) * .49 + 0.05, left=(i % 2) * .5 + 0.05,
                  top=(1 + i // 2) * .49 - 0.04, right=(1 + i % 2) * .5 - 0.05,
                  wspace=0.0, hspace=0.0)
        # Histograms
        ax = fig.add_subplot(gs[0])
        ax.set_title("#{}".format(i_dict[i] + 1), pad=0.0, fontsize=18)
        ax.set_xlabel("Exp.", fontsize=18)
        ax.bar(x_ticks[1:], exp_data_freq[1:], width=width, color=simulation_color_dict["ref"], align="edge",
               linewidth=1.0, edgecolor="k")
        removeAllBorders()
        plt.xticks([], [])
        plt.yticks([], [])
        ax = fig.add_subplot(gs[1])
        ax.set_xlabel("NORTA", fontsize=18)
        ax.bar(x_ticks[1:], norta_data_freq[1:], width=width, color=simulation_color_dict["NORTA"], align="edge",
               linewidth=1.0, edgecolor="k")
        removeAllBorders()
        plt.xticks([], [])
        plt.yticks([], [])
        # ============
        # DEPRECATED
        # plt.subplot(2, 4, i + 1)
        # exp_data_hist = exp_hist_dict[each_e]
        # exp_data_freq = exp_data_hist / np.sum(exp_data_hist)
        # norta_data_hist = sim_hist_dict["NORTA-{}".format(each_e)]
        # norta_data_freq = norta_data_hist / np.sum(norta_data_hist)
        # # -----
        # pval_list.append(scipy.stats.ks_2samp(exp_data_freq, norta_data_freq).pvalue)
        # ise_list.append(ISE(exp_data_freq, norta_data_freq))
        # # -----
        # plt.bar(x_ticks - 1 * width, exp_data_freq, width=width, label="Experimental",
        #         color=simulation_color_dict["ref"], align="edge", linewidth=1.0, edgecolor="k")
        # plt.bar(x_ticks, norta_data_freq, width=width, label="NORTA", color=simulation_color_dict["NORTA"],
        #         align="edge", linewidth=1.0, edgecolor="k")
        # plt.xticks([], [])
        # plt.yticks([], [])
        # plt.gca().spines['right'].set_visible(False)
        # plt.gca().spines['left'].set_visible(False)
        # plt.gca().spines['top'].set_visible(False)
        # plt.gca().spines['bottom'].set_visible(False)
        # if i == 3:
        #     plt.legend(ncol=1, fontsize=14)
        # # -----
        # plt.subplot(2, 4, i + 5)
        # plt.bar(x_ticks[1:] - 1 * width, exp_data_freq[1:], width=width, label="Experimental",
        #         color=simulation_color_dict["ref"], align="edge", linewidth=1.0, edgecolor="k")
        # plt.bar(x_ticks[1:], norta_data_freq[1:], width=width, label="NORTA", color=simulation_color_dict["NORTA"],
        #         align="edge", linewidth=1.0, edgecolor="k")
        # plt.xticks([], [])
        # plt.yticks([], [])
        # plt.gca().spines['right'].set_visible(False)
        # plt.gca().spines['left'].set_visible(False)
        # plt.gca().spines['top'].set_visible(False)
        # plt.gca().spines['bottom'].set_visible(False)
        # ============
    plt.tight_layout()
    plt.savefig("../fig/Fig5/distribution_shape_NORTA.pdf")
    plt.show()
    # -----
    print("p-value: ", pval_list)
    print("ISE: ", ise_list)
    # -----
    iae_list = []
    for i, each_e in enumerate(exp_data_name):
        exp_data_hist = exp_hist_dict[each_e]
        exp_data_freq = exp_data_hist / np.sum(exp_data_hist)
        norta_data_hist = sim_hist_dict["NORTA-{}".format(each_e)]
        norta_data_freq = norta_data_hist / np.sum(norta_data_hist)
        # -----
        iae_list.append(IAE(exp_data_freq, norta_data_freq))
    print("Avg IAE (NORTA) = {}".format(np.mean(iae_list)))


def compareSERGIODistShape():
    '''
    Compare distribution shape between scaled experimental data and SERGIO simulations of 100 HVGs..
    '''
    print("=" * 70)
    print("Experimental vs. SERGIO Simulation | Distribution Shape")
    exp_hist_dict = np.load("../fig_res/experimental-hist.npy", allow_pickle=True).item()["100hvg"]
    sim_hist_dict = np.load("../fig_res/new_sim-hist.npy", allow_pickle=True).item()["100hvg"]
    exp_data_name = (
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"])] +
            ["{}-{}-100hvg".format(each[0], each[1]) for each in
             itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"])]
    )
    # exp_data_name = [exp_data_name[5], exp_data_name[1]]
    pval_list = []  # two sample ks test p-value
    ise_list = []  # integrated squared error
    # Plotting
    fig = plt.figure(figsize=(10, 8))
    x_ticks = np.linspace(0.0, 1.0, 20)
    width = (x_ticks[1] - x_ticks[0])
    i_dict = {0: 2, 1: 3, 2: 0, 3: 1}
    for i, each_e in enumerate(exp_data_name[:4]):
        exp_data_hist = exp_hist_dict[exp_data_name[i_dict[i]]]
        exp_data_freq = exp_data_hist / np.sum(exp_data_hist)
        sergio_data_hist = sim_hist_dict["SERGIO-1-100hvg".format(each_e)] if i == 0 else sim_hist_dict[
            "SERGIO-15-100hvg".format(each_e)]
        sergio_data_freq = sergio_data_hist / np.sum(sergio_data_hist)
        # -----
        pval_list.append(scipy.stats.ks_2samp(exp_data_freq, sergio_data_freq).pvalue)
        ise_list.append(ISE(exp_data_freq, sergio_data_freq))
        # -----
        # -----
        outergs = gridspec.GridSpec(1, 1)
        outergs.update(bottom=(i // 2) * .49 + 0.01, left=(i % 2) * .5 + 0.01,
                       top=(1 + i // 2) * .49 - 0.01, right=(1 + i % 2) * .5 - 0.01)
        outerax = fig.add_subplot(outergs[0])
        outerax.tick_params(axis='both', which='both', bottom=0, left=0, labelbottom=0, labelleft=0)
        # Add inner grid
        gs = gridspec.GridSpec(1, 2)
        gs.update(bottom=(i // 2) * .49 + 0.05, left=(i % 2) * .5 + 0.05,
                  top=(1 + i // 2) * .49 - 0.04, right=(1 + i % 2) * .5 - 0.05,
                  wspace=0.0, hspace=0.0)
        # Histograms
        ax = fig.add_subplot(gs[0])
        ax.set_title("#{}".format(i_dict[i] + 1), pad=0.0, fontsize=18)
        ax.set_xlabel("Exp.", fontsize=18)
        ax.bar(x_ticks[1:], exp_data_freq[1:], width=width, color=simulation_color_dict["ref"], align="edge",
               linewidth=1.0, edgecolor="k")
        removeAllBorders()
        plt.xticks([], [])
        plt.yticks([], [])
        ax = fig.add_subplot(gs[1])
        ax.set_xlabel("SERGIO", fontsize=18)
        ax.bar(x_ticks[1:], sergio_data_freq[1:], width=width, color=simulation_color_dict["SERGIO"], align="edge",
               linewidth=1.0, edgecolor="k")
        removeAllBorders()
        plt.xticks([], [])
        plt.yticks([], [])
        # =============
        # DEPRECATED
        # plt.bar(x_ticks - 1 * width, exp_data_freq, width=width, label="Experimental",
        #         color=simulation_color_dict["ref"], align="edge", linewidth=1.0, edgecolor="k")
        # plt.bar(x_ticks, sergio_data_freq, width=width, label="SERGIO", color=simulation_color_dict["SERGIO"],
        #         align="edge", linewidth=1.0, edgecolor="k")
        # plt.xticks([], [])
        # plt.yticks([], [])
        # plt.gca().spines['right'].set_visible(False)
        # plt.gca().spines['left'].set_visible(False)
        # plt.gca().spines['top'].set_visible(False)
        # plt.gca().spines['bottom'].set_visible(False)
        # # -----
        # plt.subplot(1, 4, i + 3)
        # plt.bar(x_ticks[1:] - 1 * width, exp_data_freq[1:], width=width, label="Experimental",
        #         color=simulation_color_dict["ref"], align="edge", linewidth=1.0, edgecolor="k")
        # plt.bar(x_ticks[1:], sergio_data_freq[1:], width=width, label="SERGIO", color=simulation_color_dict["SERGIO"],
        #         align="edge", linewidth=1.0, edgecolor="k")
        # plt.xticks([], [])
        # plt.yticks([], [])
        # plt.gca().spines['right'].set_visible(False)
        # plt.gca().spines['left'].set_visible(False)
        # plt.gca().spines['top'].set_visible(False)
        # plt.gca().spines['bottom'].set_visible(False)
        # if i == 1:
        #     plt.legend(ncol=1, fontsize=14)
        # =============
    plt.tight_layout()
    plt.savefig("../fig/Fig5/distribution_shape_SERGIO.pdf")
    plt.show()
    # -----
    print("p-value: ", pval_list)
    print("ISE: ", ise_list)
    # -----
    # -----
    iae_list = []
    for i, each_e in enumerate(exp_data_name):
        exp_data_hist = exp_hist_dict[each_e]
        exp_data_freq = exp_data_hist / np.sum(exp_data_hist)
        for each_s in sim_hist_dict:
            sergio_data_hist = sim_hist_dict[each_s]
            sergio_data_freq = sergio_data_hist / np.sum(sergio_data_hist)
            iae_list.append(IAE(exp_data_freq, sergio_data_freq))
    print("Avg IAE (SERGIO) = {}".format(np.mean(iae_list)))

# ========================================================
#          Fig 6: clustering of simulations
# ========================================================

def experimentalVsNORTALowDim():
    num_genes = "100hvg"
    pbmc_path = "../data/experimental/PBMC/processed/expr/"
    cortex_path = "../data/experimental/mouse_cortex/processed/expr/"
    norta_path = "../data/simulation/NORTA/"
    data_name = (
            ["{}-{}-{}".format(each[0], each[1], num_genes) for each in
             itertools.product(["Cortex1", "Cortex2"], ["10xChromium", "Smart_seq2"])] +
            ["{}-{}-{}".format(each[0], each[1], num_genes) for each in
             itertools.product(["pbmc1", "pbmc2"], ["Drop", "inDrops"])]
    )
    for i, d in enumerate(data_name):
        print("[ {} ] Embedding".format(d))
        exp_file_name = "{}/{}.csv".format(pbmc_path if "pbmc" in d else cortex_path, d)
        norta_file_name = "{}/{}-NORTA-data_mat.csv".format(norta_path, d)
        exp_data_mat = loadData(exp_file_name)
        norta_data_mat = loadData(norta_file_name)
        print("Exp data shape: ", exp_data_mat.shape)
        print("Norta data shape: ", norta_data_mat.shape)
        # -----
        exp_adata = scanpy.AnnData(exp_data_mat)
        norta_adata = scanpy.AnnData(norta_data_mat)

        scanpy.pp.neighbors(exp_adata, n_neighbors=20, n_pcs=50)
        scanpy.pp.neighbors(norta_adata, n_neighbors=20, n_pcs=50)

        scanpy.tl.pca(exp_adata, svd_solver='arpack')
        scanpy.tl.pca(norta_adata, svd_solver='arpack')

        scanpy.tl.leiden(exp_adata)
        scanpy.tl.leiden(norta_adata)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
        scanpy.tl.paga(exp_adata)
        scanpy.pl.paga(exp_adata, plot=False)
        scanpy.tl.umap(exp_adata, init_pos='paga')
        exp_s_score = silhouette_score(exp_adata.obsm['X_pca'], exp_adata.obs["leiden"])
        scanpy.pl.umap(exp_adata, color="leiden", legend_loc='on data', add_outline =True, legend_fontsize="x-large", legend_fontoutline=5, title='', frameon=False, show=False, ax=ax1, s=50)
        ax1.set_title("Silhouette={:.3f}".format(exp_s_score))

        scanpy.tl.paga(norta_adata)
        scanpy.pl.paga(norta_adata, plot=False)
        scanpy.tl.umap(norta_adata, init_pos='paga')
        norta_s_score = silhouette_score(norta_adata.obsm['X_pca'], norta_adata.obs["leiden"])
        scanpy.pl.umap(norta_adata, color="leiden", legend_loc='on data', add_outline =True, legend_fontsize="x-large", legend_fontoutline=5, title='', frameon=False, show=False, ax=ax2, s=50)
        ax2.set_title("Silhouette={:.3f}".format(norta_s_score))
        plt.tight_layout()
        plt.savefig("../fig/Fig6/exp_vs_NORTA_{}.pdf".format(i+1))
        plt.show()


def SERGIOLowDim():
    sergio_path = "../data/simulation/SERGIO/"
    data_name = [
        "100gene-9groups-1","100gene-9groups-5","100gene-9groups-10",
        "100gene-9groups-15","100gene-9groups-20"
    ]
    for i, d in enumerate(data_name):
        print("[ {} ] Embedding".format(d))
        sergio_file_name = "{}/{}sparsity.csv".format(sergio_path, d)
        sergio_data_mat = loadData(sergio_file_name)
        print("SERGIO data shape: ", sergio_data_mat.shape)
        # -----
        sergio_adata = scanpy.AnnData(sergio_data_mat)
        scanpy.pp.neighbors(sergio_adata, n_neighbors=20, n_pcs=50)
        scanpy.tl.pca(sergio_adata, svd_solver='arpack')
        scanpy.tl.leiden(sergio_adata)

        fig, ax1 = plt.subplots(1, 1, figsize=(7.5, 7))
        scanpy.tl.paga(sergio_adata)
        scanpy.pl.paga(sergio_adata, plot=False)
        scanpy.tl.umap(sergio_adata, init_pos='paga')
        s_score = silhouette_score(sergio_adata.obsm['X_pca'], sergio_adata.obs["leiden"])
        scanpy.pl.umap(sergio_adata, color="leiden", legend_loc='on data', add_outline =True, legend_fontsize="x-large", legend_fontoutline=5, title='', frameon=False, show=False, ax=ax1, s=50)
        ax1.set_title("Silhouette={:.3f}".format(s_score))
        plt.tight_layout()
        plt.savefig("../fig/Fig6/SERGIO_{}.pdf".format(i+1))
        plt.show()


# ========================================================
#          Fig 7: method performance
# ========================================================

all_models_list = ["pearson", "spearman", "glasso", "GENIE3", "minet", "PIDC", "scDesign2", "scLink", "ZILGM"]
bulk_model_list = ["pearson", "spearman", "glasso", "GENIE3", "minet"]
sc_model_list = ["PIDC", "scDesign2", "scLink", "ZILGM"]

def performanceOnNORTA():
    '''
    Compare model performance on NORTA simulations
    '''
    print("Model Performance on NORTA Simulations")
    metric = np.load("../fig_res/NORTA_sim_metric.npy", allow_pickle=True).item()["100hvg"]
    F1_metric = metric["F1"]
    FDR_metric = metric["FDR"]
    auroc_metric = metric["AUROC"]
    auprc_metric = metric["AUPRC"]
    # -----
    model_F1_list = [
        [F1_metric[x][m] for x in F1_metric if (m in F1_metric[x].index) and (not pd.isna(F1_metric[x][m]))] for m in
        all_models_list]
    model_FDR_list = [
        [FDR_metric[x][m] for x in FDR_metric if (m in FDR_metric[x].index) and (not pd.isna(FDR_metric[x][m]))] for m
        in all_models_list]
    model_PPV_list = [
        [1 - FDR_metric[x][m] for x in FDR_metric if (m in FDR_metric[x].index) and (not pd.isna(FDR_metric[x][m]))] for
        m in all_models_list]
    model_AUROC_list = [
        [auroc_metric[x][m] for x in auroc_metric if (m in auroc_metric[x].index) and (not pd.isna(auroc_metric[x][m]))]
        for m in all_models_list]
    model_AUPRC_list = [
        [auprc_metric[x][m] for x in auprc_metric if (m in auprc_metric[x].index) and (not pd.isna(auprc_metric[x][m]))]
        for m in all_models_list]
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
    plt.vlines(x=[5.5], ymin=0.0, ymax=np.max(np.max(model_F1_list)) + 0.05, linestyles="solid", colors=gray_color)
    plt.ylim(0.0, np.max(np.max(model_F1_list)) + 0.05)
    plt.tight_layout()
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.savefig("../fig/Fig7/NORTA_F1_score.pdf")
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
    plt.vlines(x=[5.5], ymin=0.0, ymax=np.max(np.max(model_AUROC_list)) + 0.05, linestyles="solid", colors=gray_color)
    plt.ylim(0.4, np.max(np.max(model_AUROC_list)) + 0.05)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/Fig7/NORTA_AUROC.pdf")
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
    plt.vlines(x=[5.5], ymin=0.0, ymax=np.max(np.max(model_AUPRC_list)) + 0.05, linestyles="solid", colors=gray_color)
    plt.ylim(0.0, np.max(np.max(model_AUPRC_list)) + 0.05)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/Fig7/NORTA_AUPRC.pdf")
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
    plt.legend(bbox_to_anchor=(0.3, 1.0), ncol=2)
    plt.tight_layout()
    plt.savefig("../fig/Fig7/NORTA_FDR.pdf")
    plt.show()


def performanceOnSERGIO():
    '''
    Compare model performance on SERGIO simulations
    '''
    print("Model Performance on SERGIO Simulations")
    metric = np.load("../fig_res/SERGIO_sim_metric.npy", allow_pickle=True).item()["100hvg"]
    F1_metric = metric["F1"]
    FDR_metric = metric["FDR"]
    auroc_metric = metric["AUROC"]
    auprc_metric = metric["AUPRC"]
    # -----
    model_F1_list = [
        [F1_metric[x][m] for x in F1_metric if (m in F1_metric[x].index) and (not pd.isna(F1_metric[x][m]))] for m in
        all_models_list]
    model_FDR_list = [
        [FDR_metric[x][m] for x in FDR_metric if (m in FDR_metric[x].index) and (not pd.isna(FDR_metric[x][m]))] for m
        in all_models_list]
    model_PPV_list = [
        [1 - FDR_metric[x][m] for x in FDR_metric if (m in FDR_metric[x].index) and (not pd.isna(FDR_metric[x][m]))] for
        m in all_models_list]
    model_AUROC_list = [
        [auroc_metric[x][m] for x in auroc_metric if (m in auroc_metric[x].index) and (not pd.isna(auroc_metric[x][m]))]
        for m in all_models_list]
    model_AUPRC_list = [
        [auprc_metric[x][m] for x in auprc_metric if (m in auprc_metric[x].index) and (not pd.isna(auprc_metric[x][m]))]
        for m in all_models_list]
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
    plt.vlines(x=[5.5], ymin=0.0, ymax=np.max(np.max(model_F1_list)) + 0.05, linestyles="solid", colors=gray_color)
    plt.ylim(0.025, np.max(np.max(model_F1_list)) + 0.05)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/Fig7/SERGIO_F1_score.pdf")
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
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.75, linestyles="solid", colors=gray_color)
    plt.ylim(0.4, 0.75)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=plt.gca().get_ylim()[1] - 0.01, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=plt.gca().get_ylim()[1] - 0.01, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/Fig7/SERGIO_AUROC.pdf")
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
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.14, linestyles="solid", colors=gray_color)
    plt.ylim(0.025, 0.14)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=0.14, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=0.14, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/Fig7/SERGIO_AUPRC.pdf")
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
    plt.legend(bbox_to_anchor=(0.3, 1.0), ncol=2)
    plt.tight_layout()
    plt.savefig("../fig/Fig7/SERGIO_FDR.pdf")
    plt.show()


def performanceOnExperimental():
    '''
    Compare model performance on real data
    '''
    print("Model Performance on Real Data")
    metric = np.load("../fig_res/experimental_metric.npy", allow_pickle=True).item()
    F1_metric = metric["F1"]
    FDR_metric = metric["FDR"]
    # -----
    model_F1_list = [
        [F1_metric[x][m] for x in F1_metric if (m in F1_metric[x].index) and (not pd.isna(F1_metric[x][m]))] for m in
        all_models_list]
    model_FDR_list = [
        [FDR_metric[x][m] for x in FDR_metric if (m in FDR_metric[x].index) and (not pd.isna(FDR_metric[x][m]))] for m
        in all_models_list]
    # -----
    print(all_models_list)
    print("Avg F1: ", [np.mean(each) for each in model_F1_list])
    print("Std F1: ", [np.std(each) for each in model_F1_list])
    print("Avg FDR: ", [np.mean(each) for each in model_FDR_list])
    print("Std FDR: ", [np.std(each) for each in model_FDR_list])
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
    plt.vlines(x=[5.5], ymin=0.0, ymax=0.025, linestyles="solid", colors=gray_color)
    plt.yticks([0.0, 0.005, 0.010, 0.015, 0.020, 0.025], [0.0, 0.005, 0.010, 0.015, 0.020, 0.025])
    plt.ylim(0.0, 0.025)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=0.025, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=0.025, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/Fig7/experimental_F1_score.pdf")
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
    plt.legend(bbox_to_anchor=(0.3, 1.0), ncol=2)
    plt.tight_layout()
    plt.savefig("../fig/Fig7/experimental_FDR.pdf")
    plt.show()


def performanceMouseCortex():
    '''
    Compare model performance on real data
    '''
    print("Model Performance on mouse cortex data ")
    metric = np.load("../fig_res/moust_TF_PCC-cortex-metric.npy", allow_pickle=True).item()
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
    plt.ylim(0.0, 0.4)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=0.41, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=0.41, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/Fig7/mouse_cortex-F1.pdf")
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
    plt.savefig("../fig/Fig7/mouse_cortex-TPR.pdf")
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
    plt.savefig("../fig/Fig7/mouse_cortex-FDR.pdf")
    plt.show()


def performancePBMC():
    '''
    Compare model performance on real data
    '''
    print("Model Performance on PBMC ")
    metric = np.load("../fig_res/human_TF_similarity-cortex-metric.npy", allow_pickle=True).item()
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
    plt.ylim(0.0, 0.15)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.text(x=3, y=0.16, s="bulk RNA-seq/microarray", color="k", fontsize=20,
             horizontalalignment="center")
    plt.text(x=7.5, y=0.16, s="scRNA-seq", color="k",
             fontsize=20, horizontalalignment="center")
    plt.tight_layout()
    plt.savefig("../fig/Fig7/PBMC-F1.pdf")
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
    plt.savefig("../fig/Fig7/PBMC-TPR.pdf")
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
    plt.savefig("../fig/Fig7/PBMC-FDR.pdf")
    plt.show()
    # -----


# ========================================================
#       Fig 9: simulation w/ different settings
# ========================================================

def diffCellNum():
    '''
    Simulations with different number of cells
    '''
    print("Model Performance on NORTA Simulations w/ Different Cell Num")
    # -----
    F1_metric = np.load("../fig_res/diff_cell_metric.npy", allow_pickle=True).item()["F1"]
    cell_num_list = [0.1, 0.5, 1.0, 1.5, 2.0]
    F1_mat = np.zeros((len(cell_num_list), len(all_models_list)))
    F1_mat[F1_mat == 0] = np.nan
    for i, each_c in enumerate(cell_num_list):
        for j, each_m in enumerate(all_models_list):
            tmp_F1 = F1_metric["Cortex1-10xChromium-100hvg-{}cell-metrics.csv".format(each_c)]
            if each_m in tmp_F1.index:
                F1_mat[i, j] = tmp_F1[each_m]
    # -----
    F1_df = pd.DataFrame(data=F1_mat, index=cell_num_list, columns=all_models_list)
    print("-" * 70)
    print("F1 score:")
    print(F1_df.to_string())
    # -----
    x_ticks = np.arange(len(cell_num_list))
    x_tick_labels = cell_num_list
    plt.figure(figsize=(8, 4))
    for i, x in enumerate(all_models_list):
        plt.plot(x_ticks, F1_mat[:, i], "o-." if x in bulk_model_list else "s-", lw=4, ms=10, color=diff_model_color_dict[x], label=model_name_dict[x].replace("\n", ""))
    plt.xticks(x_ticks, x_tick_labels)
    plt.xlabel("# cells / # genes")
    plt.ylabel("F1 score")
    plt.legend(ncol=3, fontsize=14)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig9/cell_num_F1.pdf")
    plt.show()
    # -----
    auroc_metric = np.load("../fig_res/diff_cell_metric.npy", allow_pickle=True).item()["AUROC"]
    auroc_mat = np.zeros((len(cell_num_list), len(all_models_list)))
    auroc_mat[auroc_mat == 0] = np.nan
    for i, each_c in enumerate(cell_num_list):
        for j, each_m in enumerate(all_models_list):
            tmp_auroc = auroc_metric["Cortex1-10xChromium-100hvg-{}cell-metrics.csv".format(each_c)]
            if each_m in tmp_auroc.index:
                auroc_mat[i, j] = tmp_auroc[each_m]
    # -----
    auroc_df = pd.DataFrame(data=auroc_mat, index=cell_num_list, columns=all_models_list)
    print("-" * 70)
    print("AUROC:")
    print(auroc_df.to_string())
    # -----
    x_ticks = np.arange(len(cell_num_list))
    x_tick_labels = cell_num_list
    plt.figure(figsize=(8, 4))
    for i, x in enumerate(all_models_list):
        plt.plot(x_ticks, auroc_mat[:, i], "o-." if x in bulk_model_list else "s-", lw=4, ms=10, color=diff_model_color_dict[x], label=model_name_dict[x].replace("\n", ""))
    plt.xticks(x_ticks, x_tick_labels)
    plt.xlabel("# cells / # genes")
    plt.ylabel("AUROC")
    plt.ylim(0.49, 0.60)
    plt.legend(ncol=3, fontsize=14)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig9/cell_num_AUROC.pdf")
    plt.show()
    # -----
    auprc_metric = np.load("../fig_res/diff_cell_metric.npy", allow_pickle=True).item()["AUPRC"]
    auprc_mat = np.zeros((len(cell_num_list), len(all_models_list)))
    auprc_mat[auprc_mat == 0] = np.nan
    for i, each_c in enumerate(cell_num_list):
        for j, each_m in enumerate(all_models_list):
            tmp_auprc = auprc_metric["Cortex1-10xChromium-100hvg-{}cell-metrics.csv".format(each_c)]
            if each_m in tmp_auprc.index:
                auprc_mat[i, j] = tmp_auprc[each_m]
    # -----
    auprc_df = pd.DataFrame(data=auprc_mat, index=cell_num_list, columns=all_models_list)
    print("-" * 70)
    print("AUPRC:")
    print(auprc_df.to_string())
    # -----
    x_ticks = np.arange(len(cell_num_list))
    x_tick_labels = cell_num_list
    plt.figure(figsize=(8, 4))
    for i, x in enumerate(all_models_list):
        plt.plot(x_ticks, auprc_mat[:, i], "o-." if x in bulk_model_list else "s-", lw=4, ms=10, color=diff_model_color_dict[x], label=model_name_dict[x].replace("\n", ""))
    plt.xticks(x_ticks, x_tick_labels)
    plt.xlabel("# cells / # genes")
    plt.ylabel("AUPRC")
    plt.legend(ncol=3, fontsize=14)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig9/cell_num_AUPRC.pdf")
    plt.show()


def diffGraph():
    '''
    Simulations with different graph structures
    '''
    print("Model Performance on NORTA Simulations w/ Different Graphs")
    # -----
    F1_metric = np.load("../fig_res/diff_graph_metric.npy", allow_pickle=True).item()["F1"]
    graph_list = ["power", "scale", "GRN"] # power is actually for hub-based results
    F1_mat = np.zeros((len(graph_list), len(all_models_list)))
    F1_mat[F1_mat == 0] = np.nan
    for i, each_c in enumerate(graph_list):
        for j, each_m in enumerate(all_models_list):
            tmp_F1 = F1_metric["Cortex1-10xChromium-100hvg-{}-metrics.csv".format(each_c)]
            if each_m in tmp_F1.index:
                F1_mat[i, j] = tmp_F1[each_m]
    # -----
    F1_df = pd.DataFrame(data=F1_mat, index=graph_list, columns=all_models_list)
    print("-" * 70)
    print("F1 score:")
    print(F1_df.to_string())
    # -----
    x_ticks = np.arange(len(graph_list)) + 1
    x_tick_labels = ["hub", "scale-free", "GRN"]
    plt.figure(figsize=(5, 4))
    F1_mat = [each[~np.isnan(each)] for each in F1_mat]
    bplt = plt.boxplot(F1_mat, patch_artist=True)
    for patch in bplt["boxes"]:
        patch.set_facecolor(gray_color)
        patch.set_linewidth(1.5)
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
    plt.xticks(x_ticks, x_tick_labels)
    plt.ylabel("F1 score")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig9/graph_F1.pdf")
    plt.show()
    # -----
    auroc_metric = np.load("../fig_res/diff_graph_metric.npy", allow_pickle=True).item()["AUROC"]
    auroc_mat = np.zeros((len(graph_list), len(all_models_list)))
    auroc_mat[auroc_mat == 0] = np.nan
    for i, each_c in enumerate(graph_list):
        for j, each_m in enumerate(all_models_list):
            tmp_auroc = auroc_metric["Cortex1-10xChromium-100hvg-{}-metrics.csv".format(each_c)]
            if each_m in tmp_auroc.index:
                auroc_mat[i, j] = tmp_auroc[each_m]
    # -----
    auroc_df = pd.DataFrame(data=auroc_mat, index=graph_list, columns=all_models_list)
    print("-" * 70)
    print("AUROC:")
    print(auroc_df.to_string())
    # -----
    plt.figure(figsize=(4, 4))
    auroc_mat = [each[~np.isnan(each)] for each in auroc_mat]
    bplt = plt.boxplot(auroc_mat, patch_artist=True)
    for patch in bplt["boxes"]:
        patch.set_facecolor(gray_color)
        patch.set_linewidth(1.5)
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
    plt.xticks(x_ticks, x_tick_labels)
    plt.ylabel("AUROC")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig9/graph_AUROC.pdf")
    plt.show()
    # -----
    auprc_metric = np.load("../fig_res/diff_graph_metric.npy", allow_pickle=True).item()["AUPRC"]
    auprc_mat = np.zeros((len(graph_list), len(all_models_list)))
    auprc_mat[auprc_mat == 0] = np.nan
    for i, each_c in enumerate(graph_list):
        for j, each_m in enumerate(all_models_list):
            tmp_auprc = auprc_metric["Cortex1-10xChromium-100hvg-{}-metrics.csv".format(each_c)]
            if each_m in tmp_auprc.index:
                auprc_mat[i, j] = tmp_auprc[each_m]
    # -----
    auprc_df = pd.DataFrame(data=auprc_mat, index=graph_list, columns=all_models_list)
    print("-" * 70)
    print("AUPRC:")
    print(auprc_df.to_string())
    # -----
    plt.figure(figsize=(4, 4))
    auprc_mat = [each[~np.isnan(each)] for each in auprc_mat]
    bplt = plt.boxplot(auprc_mat, patch_artist=True)
    for patch in bplt["boxes"]:
        patch.set_facecolor(gray_color)
        patch.set_linewidth(1.5)
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
    plt.xticks(x_ticks, x_tick_labels)
    plt.ylabel("AUPRC")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig9/graph_AUPRC.pdf")
    plt.show()


def diffSparsity():
    '''
    Simulations with different number of cells
    '''
    print("Model Performance on NORTA Simulations w/ Different Sparsity")
    # -----
    F1_metric = np.load("../fig_res/diff_sparsity_metric.npy", allow_pickle=True).item()["F1"]
    sparsity_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    F1_mat = np.zeros((len(sparsity_list), len(all_models_list)))
    F1_mat[F1_mat == 0] = np.nan
    for i, each_c in enumerate(sparsity_list):
        for j, each_m in enumerate(all_models_list):
            tmp_F1 = F1_metric["100gene-hub-{}sparsity-metrics.csv".format(each_c)]
            if each_m in tmp_F1.index:
                F1_mat[i, j] = tmp_F1[each_m]
    # -----
    F1_df = pd.DataFrame(data=F1_mat, index=sparsity_list, columns=all_models_list)
    print("-" * 70)
    print("F1 score:")
    print(F1_df.to_string())
    # -----
    x_ticks = np.arange(len(sparsity_list))
    x_tick_labels = sparsity_list
    plt.figure(figsize=(8, 4))
    for i, x in enumerate(all_models_list):
        plt.plot(x_ticks, F1_mat[:, i][::-1], "o-." if x in bulk_model_list else "s-", lw=4, ms=10, color=diff_model_color_dict[x], label=model_name_dict[x].replace("\n", ""))
    plt.xticks(x_ticks, x_tick_labels[::-1])
    plt.xlabel("Sparsity (ratio of non-zeros)")
    plt.ylabel("F1 score")
    plt.legend(ncol=3, fontsize=14)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig9/sparsity_F1.pdf")
    plt.show()
    # -----
    auroc_metric = np.load("../fig_res/diff_sparsity_metric.npy", allow_pickle=True).item()["AUROC"]
    auroc_mat = np.zeros((len(sparsity_list), len(all_models_list)))
    auroc_mat[auroc_mat == 0] = np.nan
    for i, each_c in enumerate(sparsity_list):
        for j, each_m in enumerate(all_models_list):
            tmp_auroc = auroc_metric["100gene-hub-{}sparsity-metrics.csv".format(each_c)]
            if each_m in tmp_auroc.index:
                auroc_mat[i, j] = tmp_auroc[each_m]
    # -----
    auroc_df = pd.DataFrame(data=auroc_mat, index=sparsity_list, columns=all_models_list)
    print("-" * 70)
    print("AUROC:")
    print(auroc_df.to_string())
    # -----
    plt.figure(figsize=(8, 4))
    for i, x in enumerate(all_models_list):
        plt.plot(x_ticks, auroc_mat[:, i][::-1], "o-." if x in bulk_model_list else "s-", lw=4, ms=10, color=diff_model_color_dict[x], label=model_name_dict[x].replace("\n", ""))
    plt.xticks(x_ticks, x_tick_labels[::-1])
    plt.xlabel("Sparsity (ratio of non-zeros)")
    plt.ylabel("AUROC")
    plt.legend(ncol=3, fontsize=14)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig9/sparsity_AUROC.pdf")
    plt.show()
    # -----
    auprc_metric = np.load("../fig_res/diff_sparsity_metric.npy", allow_pickle=True).item()["AUPRC"]
    auprc_mat = np.zeros((len(sparsity_list), len(all_models_list)))
    auprc_mat[auprc_mat == 0] = np.nan
    for i, each_c in enumerate(sparsity_list):
        for j, each_m in enumerate(all_models_list):
            tmp_auprc = auprc_metric["100gene-hub-{}sparsity-metrics.csv".format(each_c)]
            if each_m in tmp_auprc.index:
                auprc_mat[i, j] = tmp_auprc[each_m]
    # -----
    auprc_df = pd.DataFrame(data=auprc_mat, index=sparsity_list, columns=all_models_list)
    print("-" * 70)
    print("AUPRC:")
    print(auprc_df.to_string())
    # -----
    plt.figure(figsize=(8, 4))
    for i, x in enumerate(all_models_list):
        plt.plot(x_ticks, auprc_mat[:, i][::-1], "o-." if x in bulk_model_list else "s-", lw=4, ms=10, color=diff_model_color_dict[x], label=model_name_dict[x].replace("\n", ""))
    plt.xticks(x_ticks, x_tick_labels[::-1])
    plt.xlabel("Sparsity (ratio of non-zeros)")
    plt.ylabel("AUPRC")
    plt.legend(ncol=3, fontsize=14)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig9/sparsity_AUPRC.pdf")
    plt.show()


def diffSparsityCMD():
    '''
    Simulations with different number of cells
    '''
    print("CMD on NORTA Simulations w/ Different Sparsity")
    # -----
    metric = np.load("../fig_res/diff_sparsity_corr.npy", allow_pickle=True).item()
    CMD_metric = metric["CMD"]
    pos_diff_metric = metric["pos_abs_diff"]
    neg_diff_metric = metric["neg_abs_diff"]
    sparsity_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, "full"][::-1]
    # -----
    print("-" * 70)
    print("CMD: ", CMD_metric)
    plt.figure(figsize=(8, 4))
    ax = plt.subplot()
    sbn.stripplot(data=neg_diff_metric[::-1], ax=ax, color=light_gray_color, edgecolor="k", linewidth=0.5)
    sbn.stripplot(data=pos_diff_metric[::-1], ax=ax, color=Vivid_10.mpl_colors[-1], edgecolor="k", linewidth=0.5)
    ax.scatter([], [], color=Vivid_10.mpl_colors[-1], label="Positive", s=50)
    ax.scatter([], [], color=light_gray_color, label="Negative", s=50)
    plt.xticks(np.arange(len(sparsity_list)), sparsity_list)
    plt.xlabel("Sparsity (ratio of non-zeros)")
    plt.ylabel("abs(corr - reference)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("../fig/Fig9/sparsity_corr_diff.pdf")
    plt.show()

# ========================================================
#          Fig 10: pre-processed simulations
# ========================================================

def normalizeNORTA():
    '''
    Normalization.
    '''
    print("Performance on Normalized Data")
    processed_metric = np.load("../fig_res/NORTA_lognorm_metric.npy", allow_pickle=True).item()["100hvg"]
    sctransform_metric = np.load("../fig_res/NORTA_sctransform_metric.npy", allow_pickle=True).item()["100hvg"]
    psinorm_metric = np.load("../fig_res/NORTA_psiNorm_metric.npy", allow_pickle=True).item()["100hvg"]
    raw_metric = np.load("../fig_res/NORTA_sim_metric.npy", allow_pickle=True).item()["100hvg"]
    # -----
    processed_F1_metric = processed_metric["F1"]
    sctransform_F1_metric = sctransform_metric["F1"]
    psinorm_F1_metric = psinorm_metric["F1"]
    raw_F1_metric = raw_metric["F1"]
    # -----
    processed_F1_list = [
        [
            processed_F1_metric[x][m] for x in processed_F1_metric
            if (m in processed_F1_metric[x].index) and (not pd.isna(processed_F1_metric[x][m]))
        ] + [
            sctransform_F1_metric[x][m] for x in sctransform_F1_metric
            if (m in sctransform_F1_metric[x].index) and (not pd.isna(sctransform_F1_metric[x][m]))
        ] + [
            psinorm_F1_metric[x][m] for x in psinorm_F1_metric
            if (m in psinorm_F1_metric[x].index) and (not pd.isna(psinorm_F1_metric[x][m]))
        ]
        for m in all_models_list
    ]
    raw_F1_list = [
        [
            raw_F1_metric[x][m] for x in raw_F1_metric
            if (m in raw_F1_metric[x].index) and (not pd.isna(raw_F1_metric[x][m]))
        ] for m in all_models_list
    ]
    # -----
    print(all_models_list)
    print("-" * 70)
    print("Raw: ")
    print("Avg F1: ", [np.mean(each) for each in raw_F1_list])
    print("Std F1: ", [np.std(each) for each in raw_F1_list])
    print("-" * 70)
    print("Processed: ")
    print("Avg F1: ", [np.mean(each) for each in processed_F1_list])
    print("Std F1: ", [np.std(each) for each in processed_F1_list])
    # -----
    plt.figure(figsize=(9, 5))
    width = 0.2
    x_ticks = np.arange(1, len(all_models_list) + 1)
    bplt1 = plt.boxplot(raw_F1_list, positions=x_ticks - 0.15, widths=width, patch_artist=True)
    bplt2 = plt.boxplot(processed_F1_list, positions=x_ticks + 0.15, widths=width, patch_artist=True)
    colors = [simulation_color_dict["ref"], Vivid_10.mpl_colors[9]]
    for i, bplt in enumerate((bplt1, bplt2)):
        for patch in bplt["boxes"]:
            patch.set_facecolor(colors[i])
            patch.set_linewidth(1.5)
        for patch in bplt["fliers"]:
            patch.set_markeredgecolor(colors[i])
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
    plt.bar([1], [np.nan], color=colors[0], label="raw")
    plt.bar([1], [np.nan], color=colors[1], label="normalization")
    plt.legend(ncol=2, bbox_to_anchor=(0.1, 1.0), fontsize=18)
    plt.xticks(x_ticks, [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1 score")
    # plt.ylim(0.0, np.max(np.max(processed_F1_list)) + 0.05)
    plt.ylim(0.0, 0.55)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.savefig("../fig/Fig10/normalization-F1_score-NORTA.pdf")
    plt.tight_layout()
    plt.show()


def normalizeSERGIO():
    '''
    Normalization.
    '''
    print("Performance on Pre-Processed SERGIO Data")
    processed_metric = np.load("../fig_res/SERGIO_lognorm_metric.npy", allow_pickle=True).item()["100hvg"]
    sctransform_metric = np.load("../fig_res/SERGIO_sctransform_metric.npy", allow_pickle=True).item()["100hvg"]
    psinorm_metric = np.load("../fig_res/SERGIO_psiNorm_metric.npy", allow_pickle=True).item()["100hvg"]
    raw_metric = np.load("../fig_res/SERGIO_sim_metric.npy", allow_pickle=True).item()["100hvg"]
    # -----
    processed_F1_metric = processed_metric["F1"]
    sctransform_F1_metric = sctransform_metric["F1"]
    psinorm_F1_metric = psinorm_metric["F1"]
    raw_F1_metric = raw_metric["F1"]
    # -----
    processed_F1_list = [
        [
            processed_F1_metric[x][m] for x in processed_F1_metric
            if (m in processed_F1_metric[x].index) and (not pd.isna(processed_F1_metric[x][m]))
        ] + [
            sctransform_F1_metric[x][m] for x in sctransform_F1_metric
            if (m in sctransform_F1_metric[x].index) and (not pd.isna(sctransform_F1_metric[x][m]))
        ] + [
            psinorm_F1_metric[x][m] for x in psinorm_F1_metric
            if (m in psinorm_F1_metric[x].index) and (not pd.isna(psinorm_F1_metric[x][m]))
        ]
        for m in all_models_list
    ]
    raw_F1_list = [
        [
            raw_F1_metric[x][m] for x in raw_F1_metric
            if (m in raw_F1_metric[x].index) and (not pd.isna(raw_F1_metric[x][m]))
        ] for m in all_models_list
    ]
    # -----
    print(all_models_list)
    print("-" * 70)
    print("Raw: ")
    print("Avg F1: ", [np.mean(each) for each in raw_F1_list])
    print("Std F1: ", [np.std(each) for each in raw_F1_list])
    print("-" * 70)
    print("Processed: ")
    print("Avg F1: ", [np.mean(each) for each in processed_F1_list])
    print("Std F1: ", [np.std(each) for each in processed_F1_list])
    # -----
    plt.figure(figsize=(9, 5))
    width = 0.2
    x_ticks = np.arange(1, len(all_models_list) + 1)
    bplt1 = plt.boxplot(raw_F1_list, positions=x_ticks - 0.15, widths=width, patch_artist=True)
    bplt2 = plt.boxplot(processed_F1_list, positions=x_ticks + 0.15, widths=width, patch_artist=True)
    colors = [simulation_color_dict["ref"], Vivid_10.mpl_colors[9]]
    for i, bplt in enumerate((bplt1, bplt2)):
        for patch in bplt["boxes"]:
            patch.set_facecolor(colors[i])
            patch.set_linewidth(1.5)
        for patch in bplt["fliers"]:
            patch.set_markeredgecolor(colors[i])
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
    plt.bar([1], [np.nan], color=colors[0], label="raw")
    plt.bar([1], [np.nan], color=colors[1], label="normalization")
    plt.legend(ncol=2, bbox_to_anchor=(0.1, 1.0), fontsize=18)
    plt.xticks(x_ticks, [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1 score")
    # plt.ylim(0.025, np.max(np.max(processed_F1_list)) + 0.05)
    plt.ylim(0.025, 0.18)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.savefig("../fig/Fig10/normalization-F1_score-SERGIO.pdf")
    plt.tight_layout()
    plt.show()


def pseudobulkNORTA():
    '''
    Pseudo-bulk data (100 HVGs)
    '''
    print("Performance on pseudo-bulks")
    processed_metric = np.load("../fig_res/NORTA_pseudobulk_metric.npy", allow_pickle=True).item()["100hvg"]
    raw_metric = np.load("../fig_res/NORTA_sim_metric.npy", allow_pickle=True).item()["100hvg"]
    # -----
    processed_F1_metric = processed_metric["F1"]
    raw_F1_metric = raw_metric["F1"]
    # -----
    processed_F1_list = [
        [
            processed_F1_metric[x][m] for x in processed_F1_metric
            if (m in processed_F1_metric[x].index) and (not pd.isna(processed_F1_metric[x][m]))
        ] for m in all_models_list
    ]
    raw_F1_list = [
        [
            raw_F1_metric[x][m] for x in raw_F1_metric
            if (m in raw_F1_metric[x].index) and (not pd.isna(raw_F1_metric[x][m]))
        ] for m in all_models_list
    ]
    # -----
    print(all_models_list)
    print("-" * 70)
    print("Raw: ")
    print("Avg F1: ", [np.mean(each) for each in raw_F1_list])
    print("Std F1: ", [np.std(each) for each in raw_F1_list])
    print("-" * 70)
    print("Processed: ")
    print("Avg F1: ", [np.mean(each) for each in processed_F1_list])
    print("Std F1: ", [np.std(each) for each in processed_F1_list])
    # -----
    plt.figure(figsize=(8.5, 5))
    width = 0.2
    x_ticks = np.arange(1, len(all_models_list) + 1)
    bplt1 = plt.boxplot(raw_F1_list, positions=x_ticks - 0.15, widths=width, patch_artist=True)
    bplt2 = plt.boxplot(processed_F1_list, positions=x_ticks + 0.15, widths=width, patch_artist=True)
    colors = [simulation_color_dict["ref"], Vivid_10.mpl_colors[9]]
    for i, bplt in enumerate((bplt1, bplt2)):
        for patch in bplt["boxes"]:
            patch.set_facecolor(colors[i])
            patch.set_linewidth(1.5)
        for patch in bplt["fliers"]:
            patch.set_markeredgecolor(colors[i])
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
    plt.bar([1], [np.nan], color=colors[0], label="raw")
    plt.bar([1], [np.nan], color=colors[1], label="pseudo-bulk")
    plt.legend(ncol=2, bbox_to_anchor=(0.1, 1.0), fontsize=18)
    plt.xticks(x_ticks, [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1 score")
    plt.ylim(0.0, 0.5)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig10/pseudobulk-F1_score-NORTA.pdf")
    plt.show()


def pseudobulkSERGIO():
    '''
    Pseudo-bulk of SERGIO
    '''
    print("Performance on pseudo-bulk")
    processed_metric = np.load("../fig_res/SERGIO_pseudobulk_metric.npy", allow_pickle=True).item()["100hvg"]
    raw_metric = np.load("../fig_res/SERGIO_sim_metric.npy", allow_pickle=True).item()["100hvg"]
    # -----
    processed_F1_metric = processed_metric["F1"]
    raw_F1_metric = raw_metric["F1"]
    # -----
    processed_F1_list = [
        [
            processed_F1_metric[x][m] for x in processed_F1_metric
            if (m in processed_F1_metric[x].index) and (not pd.isna(processed_F1_metric[x][m]))
        ] for m in all_models_list
    ]
    raw_F1_list = [
        [
            raw_F1_metric[x][m] for x in raw_F1_metric
            if (m in raw_F1_metric[x].index) and (not pd.isna(raw_F1_metric[x][m]))
        ] for m in all_models_list
    ]
    # -----
    print(all_models_list)
    print("-" * 70)
    print("Raw: ")
    print("Avg F1: ", [np.mean(each) for each in raw_F1_list])
    print("Std F1: ", [np.std(each) for each in raw_F1_list])
    print("-" * 70)
    print("Processed: ")
    print("Avg F1: ", [np.mean(each) for each in processed_F1_list])
    print("Std F1: ", [np.std(each) for each in processed_F1_list])
    # -----
    plt.figure(figsize=(8.5, 5))
    width = 0.2
    x_ticks = np.arange(1, len(all_models_list) + 1)
    bplt1 = plt.boxplot(raw_F1_list, positions=x_ticks - 0.15, widths=width, patch_artist=True)
    bplt2 = plt.boxplot(processed_F1_list, positions=x_ticks + 0.15, widths=width, patch_artist=True)
    colors = [simulation_color_dict["ref"], Vivid_10.mpl_colors[9]]
    for i, bplt in enumerate((bplt1, bplt2)):
        for patch in bplt["boxes"]:
            patch.set_facecolor(colors[i])
            patch.set_linewidth(1.5)
        for patch in bplt["fliers"]:
            patch.set_markeredgecolor(colors[i])
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
    plt.bar([1], [np.nan], color=colors[0], label="raw")
    plt.bar([1], [np.nan], color=colors[1], label="pseudo-bulk")
    plt.legend(ncol=2, bbox_to_anchor=(0.1, 1.0), fontsize=18)
    plt.xticks(x_ticks, [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1 score")
    # plt.ylim(0.0, 0.25)
    plt.ylim(0.0, 0.22)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig10/pseudobulk-F1_score-SERGIO.pdf")
    plt.show()


def imputationNORTA():
    '''
    Compare performance on imputed NORTA simulations
    '''
    print("Performance on imputations")
    processed_metric = np.load("../fig_res/NORTA_imputation_metric.npy", allow_pickle=True).item()["100hvg"]
    raw_metric = np.load("../fig_res/NORTA_sim_metric.npy", allow_pickle=True).item()["100hvg"]
    # -----
    processed_F1_metric = processed_metric["F1"]
    raw_F1_metric = raw_metric["F1"]
    # -----
    processed_F1_list = [
        [
            processed_F1_metric[x][m] for x in processed_F1_metric
            if (m in processed_F1_metric[x].index) and (not pd.isna(processed_F1_metric[x][m]))
        ] for m in all_models_list
    ]
    raw_F1_list = [
        [
            raw_F1_metric[x][m] for x in raw_F1_metric
            if (m in raw_F1_metric[x].index) and (not pd.isna(raw_F1_metric[x][m]))
        ] for m in all_models_list
    ]
    # -----
    print(all_models_list)
    print("-" * 70)
    print("Raw: ")
    print("Avg F1: ", [np.mean(each) for each in raw_F1_list])
    print("Std F1: ", [np.std(each) for each in raw_F1_list])
    print("-" * 70)
    print("Processed: ")
    print("Avg F1: ", [np.mean(each) for each in processed_F1_list])
    print("Std F1: ", [np.std(each) for each in processed_F1_list])
    # -----
    plt.figure(figsize=(8.5, 5))
    width = 0.2
    x_ticks = np.arange(1, len(all_models_list) + 1)
    bplt1 = plt.boxplot(raw_F1_list, positions=x_ticks - 0.15, widths=width, patch_artist=True)
    bplt2 = plt.boxplot(processed_F1_list, positions=x_ticks + 0.15, widths=width, patch_artist=True)
    colors = [simulation_color_dict["ref"], Vivid_10.mpl_colors[9]]
    for i, bplt in enumerate((bplt1, bplt2)):
        for patch in bplt["boxes"]:
            patch.set_facecolor(colors[i])
            patch.set_linewidth(1.5)
        for patch in bplt["fliers"]:
            patch.set_markeredgecolor(colors[i])
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
    plt.bar([1], [np.nan], color=colors[0], label="raw")
    plt.bar([1], [np.nan], color=colors[1], label="imputation")
    plt.legend(ncol=2, bbox_to_anchor=(0.1, 1.0), fontsize=18)
    plt.xticks(x_ticks, [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1 score")
    plt.ylim(0.0, 0.5)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig10/imputation-F1_score-NORTA.pdf")
    plt.show()


def imputationSERGIO():
    '''
    Compare performance on imputed SERGIO data
    '''
    print("Performance on imputations")
    processed_metric = np.load("../fig_res/SERGIO_imputation_metric.npy", allow_pickle=True).item()["100hvg"]
    raw_metric = np.load("../fig_res/SERGIO_sim_metric.npy", allow_pickle=True).item()["100hvg"]
    # -----
    processed_F1_metric = processed_metric["F1"]
    raw_F1_metric = raw_metric["F1"]
    # -----
    processed_F1_list = [
        [
            processed_F1_metric[x][m] for x in processed_F1_metric
            if (m in processed_F1_metric[x].index) and (not pd.isna(processed_F1_metric[x][m]))
        ] for m in all_models_list
    ]
    raw_F1_list = [
        [
            raw_F1_metric[x][m] for x in raw_F1_metric
            if (m in raw_F1_metric[x].index) and (not pd.isna(raw_F1_metric[x][m]))
        ] for m in all_models_list
    ]
    # -----
    print(all_models_list)
    print("-" * 70)
    print("Raw: ")
    print("Avg F1: ", [np.mean(each) for each in raw_F1_list])
    print("Std F1: ", [np.std(each) for each in raw_F1_list])
    print("-" * 70)
    print("Processed: ")
    print("Avg F1: ", [np.mean(each) for each in processed_F1_list])
    print("Std F1: ", [np.std(each) for each in processed_F1_list])
    # -----
    plt.figure(figsize=(8.5, 5))
    width = 0.2
    x_ticks = np.arange(1, len(all_models_list) + 1)
    bplt1 = plt.boxplot(raw_F1_list, positions=x_ticks - 0.15, widths=width, patch_artist=True)
    bplt2 = plt.boxplot(processed_F1_list, positions=x_ticks + 0.15, widths=width, patch_artist=True)
    colors = [simulation_color_dict["ref"], Vivid_10.mpl_colors[9]]
    for i, bplt in enumerate((bplt1, bplt2)):
        for patch in bplt["boxes"]:
            patch.set_facecolor(colors[i])
            patch.set_linewidth(1.5)
        for patch in bplt["fliers"]:
            patch.set_markeredgecolor(colors[i])
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
    plt.bar([1], [np.nan], color=colors[0], label="raw")
    plt.bar([1], [np.nan], color=colors[1], label="imputation")
    plt.legend(ncol=2, bbox_to_anchor=(0.1, 1.0), fontsize=18)
    plt.xticks(x_ticks, [model_name_dict[x] for x in all_models_list], fontsize=13)
    plt.ylabel("F1 score")
    # plt.ylim(0.025, 0.2)
    plt.ylim(0.025, 0.175)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("../fig/Fig10/imputationF1_score-SERGIO.pdf")
    plt.show()


if __name__ == '__main__':
    # -----
    # Fig4: experimental data vs. old simulations (ZI-Gaussian / ZI-Poisson)
    compareOldSimCellSparsity()
    compareOldSimGeneSparsity()
    compareOldSimDistShape()
    # -----
    # Fig5: experimental data vs. new simulations (NORTA / SERGIO)
    compareNORTAGeneSparsity()
    compareNORTADistShape()
    compareNORTACorr()
    compareSERGIOGeneSparsity()
    compareSERGIODistShape()
    # -----
    # Fig6: clustering of simulations
    experimentalVsNORTALowDim()
    SERGIOLowDim()
    # -----
    # Fig7: method performance on simulated and experimental data
    performanceOnNORTA() # NORTA simulations
    performanceOnSERGIO() # SERGIO simulations
    performanceMouseCortex() # mouse cortex data
    performancePBMC() # PBMC data
    # -----
    # Fig9: simulations w/ different settings
    diffCellNum()
    diffGraph()
    diffSparsity()
    diffSparsityCMD()
    # -----
    # Fig10: pre-processing approaches
    # Normalization
    normalizeNORTA()
    normalizeSERGIO()
    # Pseudo-bulk
    pseudobulkNORTA()
    pseudobulkSERGIO()
    # Imputation
    imputationNORTA()
    imputationSERGIO()
