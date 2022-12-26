'''
Description:
    Model running for PBMC data.

Author:
    Jiaqi Zhang
'''

import sys
sys.path.append("../")
sys.path.append("/")
from Py_Utils import readH5ad, getTimeStr, addDefaultArg, getDefualtLayerSize
from Running import normalAugmentation, clusterAugmentation

# from Models.Py_Utils import readH5ad, getTimeStr, addDefaultArg, getDefualtLayerSize
# from Models.StructureModel.Running import clusterAugmentation, normalAugmentation

# ------------------------------------------------------------------------------------
# Augmentation for different sampling ratios
# def clusterAugmentExp():
#     # parameters configuration
#     config = {"need_save": True}
#     # ------------------------------------------------
#     cluster_ls = [1]
#     cluster_size = {1:5008}
#     for c in cluster_ls:
#         print("=" * 70)
#         print("Cluster {}".format(c))
#         config["train_data"] = "../../Data/splat_simulation/wo_preprocess/cluster{}_data-{}.mtx".format(c, cluster_size[c])
#         for s in [0.5, 0.25, 0.1, 0.05, 0.03, 0.01]: # different sampling ratio of training data
#             print("*" * 70)
#             print("Train data size : {}".format(s))
#             config["train_size"] = s
#             for t in range(5): # 5 trials
#                 print("#" * 70)
#                 print("TRIAL {}".format(t))
#                 config["prediction_save_path"] = \
#                     "../../Prediction/DCA/splat_simulation_cluster{}-augmented-trial{}-DCA_estimation-{}.npy".format(
#                         c, t, s)
#                 clusterAugmentation(config)

# ------------------------------------------------------------------------------------
# Train, validate, SERGIO_simulation_all

def normalAugmentExp():
    config = {
        "train_data": "../../Data/WOT/wo_preprocess/training_all.mtx",
        "need_save": True,

        "model_save_path": "/gpfs/data/rsingh47/jzhan322/Prediction/WOT-DCA_model.h5",
        "prediction_save_path": "/gpfs/data/rsingh47/jzhan322/Prediction/WOT-DCA_estimation.mtx",
        "train_data_save_path": "/gpfs/data/rsingh47/jzhan322/Prediction/WOT-DCA_train_data.mtx",
        "filter_save_path": "/gpfs/data/rsingh47/jzhan322/Prediction/WOT-DCA_filters.npy",
    }
    normalAugmentation(config)




if __name__ == '__main__':
    # clusterAugmentExp()
    normalAugmentExp()
    pass