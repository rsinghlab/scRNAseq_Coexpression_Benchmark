'''
Description:
    Util functions for imputation.

Author:
    Jiaqi Zhang <jiaqi_zhang2@brown.edu>
'''

# ======================================

def loadData(filename):
    expr_mat = pd.read_csv(filename, header=0, index_col=0)
    return expr_mat
