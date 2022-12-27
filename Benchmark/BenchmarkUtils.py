'''
Description:
    Util functions for benchmarks.

Author:
    Jiaqi Zhang <jiaqi_zhang2@brown.edu>
'''

import pandas as pd


def loadExprMat(filename):
    data_mat = pd.read_csv(filename, header=0, index_col=0)
    return data_mat


def loadNetMat(filename):
    net_mat = pd.read_csv(filename, header=0, index_col=0)
    return net_mat


def loadMetric(filename):
    metric_df = pd.read_csv(filename, header=0, index_col=0)
    return metric_df