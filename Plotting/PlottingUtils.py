'''
Description:
    Util functions for plotting.

Author:
    Jiaqi Zhang <jiaqi_zhang2@brown.edu>
'''

# =====================================
#        Network Loading
# =====================================

import pandas as pd
def loadNetMat(filename):
    net_mat = pd.read_csv(filename, header=0, index_col=0)
    return net_mat

# =====================================
#       Plotting Configuration
# =====================================
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sbn

# Params configuration
params = {
    "legend.fontsize": 20,
    "legend.frameon": False,
    "ytick.labelsize": 20,
    "xtick.labelsize": 20,
    "figure.dpi": 600,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    # "pdf.fonttype": 42,
    # "font.sans-serif": "CMU Serif",
    # "font.family": "sans-serif",
    "font.sans-serif": "Arial",
    "font.family": "sans-serif",
    # "font.weight": "bold",
    "axes.unicode_minus": False,
}
plt.rcParams.update(params)

# Color palette
gray_color = (173 / 255, 181 / 255, 189 / 255)
light_gray_color = (192 / 255, 192 / 255, 192 / 255)
white_color = (1.0, 1.0, 1.0)
from palettable.tableau import BlueRed_12
from palettable.cartocolors.qualitative import Vivid_10

# Simulation colors
simulation_color_dict = {
    "ref": gray_color,
    "log_ref": gray_color,

    "NORTA": Vivid_10.mpl_colors[-1],
    "log_NORTA": Vivid_10.mpl_colors[-1],

    "netSmooth": Vivid_10.mpl_colors[0],
    "log_netSmooth": Vivid_10.mpl_colors[0],

    "scLink": Vivid_10.mpl_colors[1],
    "log_scLink": Vivid_10.mpl_colors[1],

    "ZILGM": Vivid_10.mpl_colors[2],
    "log_ZILGM": Vivid_10.mpl_colors[2],

    "SERGIO": Vivid_10.mpl_colors[3],
    "log_SERGIO": Vivid_10.mpl_colors[3],
}

# Model Colors
diff_model_color_dict = {
    "pearson": Vivid_10.mpl_colors[0],
    "spearman": Vivid_10.mpl_colors[1],
    "glasso": Vivid_10.mpl_colors[2],
    "scLink": Vivid_10.mpl_colors[3],
    "GENIE3": Vivid_10.mpl_colors[4],
    "PIDC": Vivid_10.mpl_colors[5],
    "scDesign2": Vivid_10.mpl_colors[6],
    "minet": Vivid_10.mpl_colors[7],
    "ZILGM": Vivid_10.mpl_colors[8],
}
model_color_dict = {
    "pearson": gray_color,
    "spearman": gray_color,
    "glasso": gray_color,
    "scLink": gray_color,
    "GENIE3": gray_color,
    "PIDC": gray_color,
    "scDesign2": gray_color,
    "minet": gray_color,
    "ZILGM": gray_color,
}

# Model names
model_name_dict = {
    "pearson": "Pearson",
    "spearman": "Spearman",
    "glasso": "GLasso\n(Cov)",
    "scLink": "scLink",
    "GENIE3": "GENIE3",
    "PIDC": "PIDC",
    "scDesign2": "GLasso\n(scDesign2)",
    "minet": "MINet",
    "ZILGM": "ZILNBGM",
}

# =====================================
#     Data Properties Computation
# =====================================
import numpy as np


def ISE(list1, list2):
    list1 = np.asarray(list1)
    list2 = np.asarray(list2)
    return np.sum(np.square(list1 - list2))


def IAE(list1, list2):
    list1 = np.asarray(list1)
    list2 = np.asarray(list2)
    return np.sum(np.abs(list1 - list2))


def _CMD(corr1, corr2):
    CMD = 1 - np.trace(corr1 @ corr2) / (np.linalg.norm(corr1, ord="fro") * np.linalg.norm(corr2, ord="fro"))
    return CMD


def _MSE(corr1, corr2):
    MSE = np.linalg.norm(corr1 - corr2, ord="fro") / np.prod(corr1.shape)
    return MSE


def computeCorrSim(ref_cor, sim_cor, true_net, type):
    # Compute similarity score between correlation matrix
    if isinstance(ref_cor, float) or isinstance(sim_cor, float):
        return None, None
    ref_cor = ref_cor[type]
    sim_cor = sim_cor[type]
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


# =====================================
#       Simulation Clustering
# =====================================
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from umap import UMAP
import itertools
import scanpy

def loadData(filename):
    data = pd.read_csv(filename, header=0, index_col=0)
    return data



# =====================================
#          Utils for Plotting
# =====================================

def removeAllBorders():
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)


