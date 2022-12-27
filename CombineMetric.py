import pandas as pd
import numpy as np
import os


def loadMetric(dir_path, filename):
    if filename in os.listdir(dir_path):
        metric = pd.read_csv("{}/{}".format(dir_path, filename), header=0, index_col=0)
        return metric
    else:
        return pd.DataFrame()


NORTA_sim_list = [
    "pbmc1-Drop", "pbmc2-Drop", "pbmc1-inDrops", "pbmc2-inDrops",
    "Cortex1-10xChromium", "Cortex2-10xChromium", "Cortex1-Smart_seq2", "Cortex2-Smart_seq2"
]
SERGIO_sim_list = [
    "100gene-9groups-1", "100gene-9groups-5", "100gene-9groups-10",
    "100gene-9groups-15", "100gene-9groups-20"
]
ZILGM_sim_list = ["100hvg-ZILGM_0.8pi", "100hvg-ZILGM_0.9pi", "100hvg-ZILGM_1.0pi"]
scLink_sim_list = ["100hvg-scLink_0.07rho", "100hvg-scLink_0.10rho", "100hvg-scLink_0.13rho"]

#
# for each in NORTA_sim_list:
#     tmp_name = "{}-100hvg".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/simulation/", "{}-NORTA-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_new_sim/", "{}-NORTA-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/simulation/NORTA/{}-NORTA-metrics.csv".format(tmp_name))


# for each in SERGIO_sim_list:
#     tmp_name = each
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_simulation_all/", "{}-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_simulation_all/", "{}-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_model_SERGIO/", "{}sparsity-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/simulation/SERGIO/{}-metrics.csv".format(tmp_name))


# for each in scLink_sim_list:
#     tmp_name = each
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/scLink_simulation/", "{}-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/scLink_simulation/", "{}-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/scLink_simulation/", "{}-ZILGM-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/simulation/ZI_Gaussian/{}-metrics.csv".format(tmp_name))


# for each in ZILGM_sim_list:
#     tmp_name = each
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_simulation/", "{}-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_simulation/", "{}-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_simulation/", "{}-ZILGM-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/simulation/ZI_Poisson/{}-metrics.csv".format(tmp_name))


# for each in [0.1, 0.5, 1.0, 1.5, 2.0]:
#     tmp_name = "Cortex1-10xChromium-100hvg-{}cell".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/diff_cell/", "{}-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/diff_cell/", "{}-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/diff_cell/", "{}-ZILGM-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/diff_settings/diff_cell/{}-metrics.csv".format(tmp_name))


# for each in ["scale", "power", "GRN"]:
#     tmp_name = "Cortex1-10xChromium-100hvg-{}".format(each)
#     save_name = "Cortex1-10xChromium-100hvg-{}".format(each if each != "scale" else "hub")
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/diff_graph/", "{}-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/diff_graph/", "{}-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/diff_graph/", "{}-ZILGM-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/diff_settings/diff_graph/{}-metrics.csv".format(save_name))


# for each in np.arange(0.1, 1.0, 0.1):
#     tmp_name = "100gene-hub-{:.1f}sparsity".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/diff_sparsity/", "{}-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/diff_sparsity/", "{}-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/diff_sparsity/", "{}-ZILGM-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/diff_settings/diff_sparsity/{}-metrics.csv".format(tmp_name))


# for each in [
#     "pbmc1-inDrops", "pbmc2-inDrops",
#     "Cortex1-Smart_seq2", "Cortex2-Smart_seq2"
# ]:
#     tmp_name = "{}-500hvg".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/high_dim/", "{}-NORTA-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/high_dim/", "{}-NORTA-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_model_highdim/", "{}-NORTA-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/diff_settings/high_dim/NORTA/{}-NORTA-metrics.csv".format(tmp_name))


# for each in [
#     "400gene-9groups-1", "400gene-9groups-5", "400gene-9groups-10",
#     "400gene-9groups-15", "400gene-9groups-20"
# ]:
#     tmp_name = each
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_400genes/", "{}-SERGIO-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_400genes/", "{}-SERGIO-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_400genes/", "{}-SERGIO-ZILGM-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/diff_settings/high_dim/SERGIO/{}-metrics.csv".format(tmp_name))



# for each in NORTA_sim_list:
#     tmp_name = "{}-100hvg".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/simulation/", "{}-NORTA-process-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_new_sim/", "{}-NORTA-process-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/pre-processing/normalization/NORTA/{}-NORTA-lognorm-metrics.csv".format(tmp_name))


# for each in NORTA_sim_list:
#     tmp_name = "{}-100hvg".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/other_norm/", "{}-sctransform-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/other_norm/", "{}-sctransform-ZILGM-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/pre-processing/normalization/NORTA/{}-NORTA-sctransform-metrics.csv".format(tmp_name))


# for each in NORTA_sim_list:
#     tmp_name = "{}-100hvg".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/other_norm/", "{}-psiNorm-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/other_norm/", "{}-psiNorm-ZILGM-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/pre-processing/normalization/NORTA/{}-NORTA-psiNorm-metrics.csv".format(tmp_name))

# for each in SERGIO_sim_list:
#     tmp_name = "{}".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_simulation_all/", "{}-process-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_simulation_all/", "{}-process-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_model_SERGIO/", "{}sparsity-process-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/pre-processing/normalization/SERGIO/{}-lognorm-metrics.csv".format(tmp_name))


# for each in SERGIO_sim_list:
#     tmp_name = "{}".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/other_norm/", "{}-sctransform-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/other_norm/", "{}-sctransform-ZILGM-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/pre-processing/normalization/SERGIO/{}-sctransform-metrics.csv".format(tmp_name))
#
#
# for each in SERGIO_sim_list:
#     tmp_name = "{}".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/other_norm/", "{}-psiNorm-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/other_norm/", "{}-psiNorm-ZILGM-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/pre-processing/normalization/SERGIO/{}-psiNorm-metrics.csv".format(tmp_name))


# for each in NORTA_sim_list:
#     tmp_name = "{}-100hvg".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/pseudo_bulk/", "{}-NORTA--metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/pseudo_bulk/", "{}-NORTA-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_model_pseudobulk/", "{}-NORTA-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/pre-processing/pseudo_bulk/NORTA/{}-NORTA-metrics.csv".format(tmp_name))

# for each in SERGIO_sim_list:
#     tmp_name = "{}".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_pseudobulk/", "{}-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_pseudobulk/", "{}-PIDC-metrics.csv".format(tmp_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_model_SERGIO_pseudobulk/", "{}-metrics.csv".format(tmp_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/pre-processing/pseudo_bulk/SERGIO/{}-metrics.csv".format(tmp_name))


# for each in NORTA_sim_list:
#     impute_name = "ALRA"
#     tmp_name = "{}-100hvg".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/imputed_simulation/", "{}-NORTA-{}-metrics.csv".format(tmp_name, impute_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/imputed_simulation/", "{}-NORTA-{}-PIDC-metrics.csv".format(tmp_name, impute_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_model_imputation/", "{}-NORTA-{}-metrics.csv".format(tmp_name, impute_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/pre-processing/imputation/NORTA/{}-NORTA-{}-metrics.csv".format(tmp_name, impute_name))


# for each in SERGIO_sim_list:
#     impute_name = "SAVER"
#     tmp_name = "{}".format(each)
#     metric = pd.concat([
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_imputation/", "{}-{}-metrics.csv".format(tmp_name, impute_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/SERGIO_imputation/", "{}-{}-PIDC-metrics.csv".format(tmp_name, impute_name)),
#         loadMetric("../Gene_Coexpression_Benchmark/res/logs/evaluation/ZILGM_model_SERGIO_imputation/", "{}-{}-metrics.csv".format(tmp_name, impute_name)),
#     ], axis=0)
#     model_name = metric.Model.unique()
#     metric.to_csv("./res/pre-processing/imputation/SERGIO/{}-{}-metrics.csv".format(tmp_name, impute_name))