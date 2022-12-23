import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../util/SERGIO-master/")
from SERGIO.sergio import sergio


def SERGIOSimulate(sparsity_level):
    print("=" * 70)
    print("***** Sparsity level = {} *****".format(sparsity_level))
    sim = sergio(number_genes=100, number_bins=9, number_sc=300, noise_params=1, decays=0.8, sampling_state=15, noise_type='dpd')
    sim.build_graph(
        input_file_taregts='../util/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1/Interaction_cID_4.txt',
        input_file_regs='../util/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1/Regs_cID_4.txt',
        shared_coop_state=2
    )
    sim.simulate()
    expr = sim.getExpressions()
    # expr_clean = np.concatenate(expr, axis=1)

    """
    Add outlier genes
    """
    expr_O = sim.outlier_effect(expr, outlier_prob = 0.01, mean = 0.8, scale = 1)

    """
    Add Library Size Effect
    """
    libFactor, expr_O_L = sim.lib_size_effect(expr_O, mean = 4.6, scale = 0.4)

    """
    # Add Dropouts
    # """
    # binary_ind = sim.dropout_indicator(expr_O_L, shape = 6.5, percentile = 82)
    # expr_O_L_D = np.multiply(binary_ind, expr_O_L)

    # binary_ind = sim.dropout_indicator(expr, shape=6.5, percentile=100-sparsity_level)
    # expr_O_L_D = np.multiply(binary_ind, expr)
    binary_ind = sim.dropout_indicator(expr_O_L, shape=6.5, percentile=100-sparsity_level)
    expr_O_L_D = np.multiply(binary_ind, expr_O_L)

    """
    Convert to UMI count
    """
    count_matrix = sim.convert_to_UMIcounts(expr_O_L_D)

    """
    Make a 2d gene expression matrix
    """
    count_matrix = np.concatenate(count_matrix, axis=1).T
    num_cells, num_genes = count_matrix.shape
    pd.DataFrame(
        data=count_matrix,
        index=["cell{}".format(i) for i in range(num_cells)],
        columns=["gene{}".format(i) for i in range(num_genes)]
    # ).to_csv("./SERGIO_simulation/100gene-9groups-{}sparsity.csv".format(sparsity_level))
    ).to_csv("./SERGIO_simulation_all/100gene-9groups-{}sparsity.csv".format(sparsity_level))


def checkSERGIO(filename):
    expr_mat = pd.read_csv(filename, header=0, index_col=0)
    adj_expr_mat = np.sign(np.abs(expr_mat))
    num_cell, num_gene = expr_mat.shape
    print("# cells = {} | # genes = {}".format(num_cell, num_gene))
    print("sparsity={}".format(np.mean(np.sign(np.abs(expr_mat.values)))))
    # -----
    zero_per_cell = np.mean(adj_expr_mat, axis=1)
    zero_per_gene = np.mean(adj_expr_mat, axis=0)
    library_size = np.sum(expr_mat, axis=1)
    gene_avg = np.mean(expr_mat, axis=0)
    gene_var = np.var(expr_mat, axis=0)
    # -----
    plt.subplot(2, 3, 1)
    plt.title("zero per cell")
    plt.boxplot([zero_per_cell])
    plt.subplot(2, 3, 2)
    plt.title("zero per gene")
    plt.boxplot([zero_per_gene])
    plt.subplot(2, 3, 3)
    plt.title("library size")
    plt.boxplot([library_size])
    plt.subplot(2, 3, 4)
    plt.title("gene expression mean")
    plt.boxplot([gene_avg])
    plt.subplot(2, 3, 5)
    plt.title("gene expression var")
    plt.boxplot([gene_var])
    plt.subplot(2, 3, 6)
    plt.title("gene expression mean vs. var")
    plt.scatter(x=gene_var, y=gene_avg)
    plt.tight_layout()
    plt.show()
    # -----
    norm_expr_mat = (expr_mat - expr_mat.min()) / (expr_mat.max() - expr_mat.min())
    plt.subplot(1, 2, 1)
    plt.hist(norm_expr_mat.values.reshape(-1))
    plt.subplot(1, 2, 2)
    plt.hist(norm_expr_mat[norm_expr_mat > 0].values.reshape(-1))
    plt.tight_layout()
    plt.show()


def extract500Genes():
    sim = sergio(number_genes=400, number_bins=9, number_sc=300, noise_params=1, decays=0.8, sampling_state=15,
                 noise_type='dpd')
    sim.build_graph(
        input_file_taregts='../util/SERGIO-master/data_sets/De-noised_400G_9T_300cPerT_5_DS2/Interaction_cID_5.txt',
        input_file_regs='../util/SERGIO-master/data_sets/De-noised_400G_9T_300cPerT_5_DS2/Regs_cID_5.txt',
        shared_coop_state=2
    )
    sim.simulate()
    expr = sim.getExpressions()
    # expr_clean = np.concatenate(expr, axis=1)

    """
    Add outlier genes
    """
    expr_O = sim.outlier_effect(expr, outlier_prob=0.01, mean=0.8, scale=1)

    """
    Add Library Size Effect
    """
    libFactor, expr_O_L = sim.lib_size_effect(expr_O, mean=4.6, scale=0.4)
    # =====
    for sparsity_level in [1, 5, 10, 15, 20]:
        print("=" * 70)
        print("***** Sparsity level = {} *****".format(sparsity_level))
        """
        # Add Dropouts
        # """
        binary_ind = sim.dropout_indicator(expr_O_L, shape=6.5, percentile=100 - sparsity_level)
        expr_O_L_D = np.multiply(binary_ind, expr_O_L)
        """
        Convert to UMI count
        """
        count_matrix = sim.convert_to_UMIcounts(expr_O_L_D)
        """
        Make a 2d gene expression matrix
        """
        count_matrix = np.concatenate(count_matrix, axis=1).T
        num_cells, num_genes = count_matrix.shape
        print(count_matrix.shape)
        pd.DataFrame(
            data=count_matrix,
            index=["cell{}".format(i) for i in range(num_cells)],
            columns=["gene{}".format(i) for i in range(num_genes)]
            # ).to_csv("./SERGIO_simulation/100gene-9groups-{}sparsity.csv".format(sparsity_level))
        ).to_csv("./SERGIO_simulation_all/400gene-9groups-{}sparsity.csv".format(sparsity_level))


def construct400GeneNetwork():
    net_edge_idx = pd.read_csv("../util/SERGIO-master/data_sets/De-noised_400G_9T_300cPerT_5_DS2/gt_GRN.csv", header=None, index_col=None).values
    net_mat = np.zeros((400, 400))
    for each in net_edge_idx:
        net_mat[each[0]-1, each[1]-1] = 1
        net_mat[each[1]-1, each[0]-1] = 1
    net_mat[np.diag_indices(400)] = 1
    # -----
    pd.DataFrame(
        data=net_mat,
    ).to_csv("./SERGIO_simulation_all/400gene-true_network.csv")


def construct100GeneNetwork():
    net_edge_idx = pd.read_csv("../util/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1/gt_GRN.csv", header=None, index_col=None).values
    net_mat = np.zeros((100, 100))
    for each in net_edge_idx:
        net_mat[each[0]-1, each[1]-1] = 1
        net_mat[each[1]-1, each[0]-1] = 1
    net_mat[np.diag_indices(100)] = 1
    # -----
    pd.DataFrame(
        data=net_mat,
    ).to_csv("./SERGIO_simulation_all/100gene-true_network.csv")



if __name__ == '__main__':
    # SERGIOSimulate(sparsity_level = 1)
    # SERGIOSimulate(sparsity_level = 5)
    # SERGIOSimulate(sparsity_level = 10)
    # SERGIOSimulate(sparsity_level = 15)
    # SERGIOSimulate(sparsity_level = 20)
    # # -----
    # for sparsity in [1, 5, 10, 15, 20]:
    #     filename = "./SERGIO_simulation_all/100gene-9groups-{}sparsity.csv".format(sparsity)
    #     checkSERGIO(filename)
    # =========================================
    # extract500Genes()
    # construct100GeneNetwork()
    # construct400GeneNetwork()
    # for sparsity in [1, 5, 10, 15, 20]:
    #     filename = "./SERGIO_simulation_all/400gene-9groups-{}sparsity.csv".format(sparsity)
    #     checkSERGIO(filename)
    pass