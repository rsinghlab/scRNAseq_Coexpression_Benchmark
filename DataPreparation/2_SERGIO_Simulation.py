'''
Description:
    SERGIO simulation.

Author:
    Jiaqi Zhang <jiaqi_zhang2@brown.edu>
'''
import numpy as np
import pandas as pd
from SERGIO.sergio import sergio

# ======================================
#TODO: file path

def SERGIOSimulate_100G():
    # Simulate clean data
    sim = sergio(number_genes=100, number_bins=9, number_sc=300, noise_params=1, decays=0.8, sampling_state=15, noise_type='dpd')
    sim.build_graph(
        input_file_taregts='../util/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1/Interaction_cID_4.txt',
        input_file_regs='../util/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1/Regs_cID_4.txt',
        shared_coop_state=2
    )
    sim.simulate()
    expr = sim.getExpressions()
    # Add outlier genes
    expr_O = sim.outlier_effect(expr, outlier_prob = 0.01, mean = 0.8, scale = 1)
    # Add Library Size Effect
    libFactor, expr_O_L = sim.lib_size_effect(expr_O, mean = 4.6, scale = 0.4)
    for sparsity_level in [1, 5, 10, 15, 20]:
        print("=" * 70)
        print("***** Sparsity level = {} *****".format(sparsity_level))
        # Add Dropouts
        binary_ind = sim.dropout_indicator(expr_O_L, shape=6.5, percentile=100 - sparsity_level)
        expr_O_L_D = np.multiply(binary_ind, expr_O_L)
        # Convert to UMI count
        count_matrix = sim.convert_to_UMIcounts(expr_O_L_D)
        # Make a 2d gene expression matrix
        count_matrix = np.concatenate(count_matrix, axis=1).T
        num_cells, num_genes = count_matrix.shape
        pd.DataFrame(
            data=count_matrix,
            index=["cell{}".format(i) for i in range(num_cells)],
            columns=["gene{}".format(i) for i in range(num_genes)]
        ).to_csv("./SERGIO_simulation_all/100gene-9groups-{}sparsity.csv".format(sparsity_level))


def constructGRN_100G():
    net_edge_idx = pd.read_csv(
        "../util/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1/gt_GRN.csv", header=None, index_col=None
    ).values
    net_mat = np.zeros((100, 100))
    for each in net_edge_idx:
        net_mat[each[0]-1, each[1]-1] = 1
        net_mat[each[1]-1, each[0]-1] = 1
    net_mat[np.diag_indices(100)] = 1
    # -----
    pd.DataFrame(
        data=net_mat,
    ).to_csv("./SERGIO_simulation_all/100gene-true_network.csv")

# ======================================

def SERGIOSimulate_400G():
    # Simulate clean data
    sim = sergio(number_genes=400, number_bins=9, number_sc=300, noise_params=1, decays=0.8, sampling_state=15,
                 noise_type='dpd')
    sim.build_graph(
        input_file_taregts='../util/SERGIO-master/data_sets/De-noised_400G_9T_300cPerT_5_DS2/Interaction_cID_5.txt',
        input_file_regs='../util/SERGIO-master/data_sets/De-noised_400G_9T_300cPerT_5_DS2/Regs_cID_5.txt',
        shared_coop_state=2
    )
    sim.simulate()
    expr = sim.getExpressions()
    # Add outlier genes
    expr_O = sim.outlier_effect(expr, outlier_prob=0.01, mean=0.8, scale=1)
    # Add Library Size Effect
    libFactor, expr_O_L = sim.lib_size_effect(expr_O, mean=4.6, scale=0.4)
    # =====
    for sparsity_level in [1, 5, 10, 15, 20]:
        print("=" * 70)
        print("***** Sparsity level = {} *****".format(sparsity_level))
        # Add Dropouts
        binary_ind = sim.dropout_indicator(expr_O_L, shape=6.5, percentile=100 - sparsity_level)
        expr_O_L_D = np.multiply(binary_ind, expr_O_L)
        # Convert to UMI count
        count_matrix = sim.convert_to_UMIcounts(expr_O_L_D)
        # Make a 2d gene expression matrix
        count_matrix = np.concatenate(count_matrix, axis=1).T
        num_cells, num_genes = count_matrix.shape
        print(count_matrix.shape)
        pd.DataFrame(
            data=count_matrix,
            index=["cell{}".format(i) for i in range(num_cells)],
            columns=["gene{}".format(i) for i in range(num_genes)]
        ).to_csv("./SERGIO_simulation_all/400gene-9groups-{}sparsity.csv".format(sparsity_level))


def constructGRN_400G():
    net_edge_idx = pd.read_csv(
        "../util/SERGIO-master/data_sets/De-noised_400G_9T_300cPerT_5_DS2/gt_GRN.csv", header=None, index_col=None
    ).values
    net_mat = np.zeros((400, 400))
    for each in net_edge_idx:
        net_mat[each[0]-1, each[1]-1] = 1
        net_mat[each[1]-1, each[0]-1] = 1
    net_mat[np.diag_indices(400)] = 1
    # -----
    pd.DataFrame(
        data=net_mat,
    ).to_csv("./SERGIO_simulation_all/400gene-true_network.csv")

# ======================================

if __name__ == '__main__':
    SERGIOSimulate_100G()
    constructGRN_100G()
    SERGIOSimulate_400G()
    constructGRN_400G()