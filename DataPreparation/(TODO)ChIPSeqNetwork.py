import pandas as pd
import numpy as np

# # =======================================
# # Pathway Common gene similarity matrix
# # (https://maayanlab.cloud/Harmonizome/dataset/Pathway+Commons+Protein-Protein+Interactions)
# # Format the gene similarity matrix
# with open("./ChIP_seq_data/Pathway_gene_similarity_matrix_cosine.txt") as f:
#     first_line = f.readline().strip('\n')
# gene_list = first_line.split("\t")
# num_cols = len(gene_list)
# gene_list = gene_list[3:]
#
# net = pd.read_csv("./ChIP_seq_data/Pathway_gene_similarity_matrix_cosine.txt", sep="\t", skiprows=2, usecols=np.arange(3, num_cols))
# net.index = gene_list
# net.columns = gene_list
# net.to_csv("./ChIP_seq_data/formatted_Pathway_gene_similarity_matrix_cosine.csv")
# # =======================================


# # =======================================
# # BioGPS gene similarity matrix
# # (https://maayanlab.cloud/Harmonizome/dataset/BioGPS+Mouse+Cell+Type+and+Tissue+Gene+Expression+Profiles)
# # Format the gene similarity matrix
# with open("./ChIP_seq_data/BioGPS_gene_similarity_matrix_cosine.txt") as f:
#     first_line = f.readline().strip('\n')
# gene_list = first_line.split("\t")
# num_cols = len(gene_list)
# gene_list = gene_list[3:]
#
# net = pd.read_csv("./ChIP_seq_data/BioGPS_gene_similarity_matrix_cosine.txt", sep="\t", skiprows=2, usecols=np.arange(3, num_cols))
# net.index = gene_list
# net.columns = gene_list
# net.to_csv("./ChIP_seq_data/formatted_BioGPS_gene_similarity_matrix_cosine.csv")
# # =======================================


# =======================================
def loadTFList():
    # http://humantfs.ccbr.utoronto.ca/download.php
    tf_list = pd.read_csv("./ChIP_seq_data/human_TF_names_v_1.01.txt", header=None, index_col=None)
    return tf_list.values


def loadPathwaySimilarityNet():
    net = pd.read_csv("ChIP_seq_data/formatted_Pathway_gene_similarity_matrix_cosine.csv", header=0, index_col=0)
    return net


def loadBioGPSSimilarityNet():
    net = pd.read_csv("ChIP_seq_data/formatted_BioGPS_gene_similarity_matrix_cosine.csv", header=0, index_col=0)
    return net


def loadMouseTFNet():
    net = pd.read_excel("./ChIP_seq_data/mouse_TF_atlas_gene_coexpression.xlsx", sheet_name="PPC", header=0, index_col=0)
    return net


def getGeneList(data_name):
    if "pbmc" in data_name:
        data_dir_path = "./experimental/PBMC/processed/expr/"
        meta_filename = "./experimental/PBMC/SCP424/metadata/meta.txt"
    elif "Cortex" in data_name:
        data_dir_path = "./experimental/mouse_cortex/processed/expr/"
        meta_filename = "./experimental/mouse_cortex/SCP425/metadata/meta_combined.txt"
    else:
        raise ValueError("Unknown  data name {}!".format(data_name))
    with open("{}/{}-1000hvg.csv".format(data_dir_path, data_name)) as f:
        first_line = f.readline().strip('\n').strip(",")
    gene_list = first_line.split(",")
    gene_list = [each.split("_")[-1] for each in gene_list]
    return gene_list


def getTabulaGeneList(data_name):
    data_dir_path = "./Tabula_Muris/500hvg/"
    with open("{}/{}.csv".format(data_dir_path, data_name)) as f:
        first_line = f.readline().strip('\n').strip(",")
    gene_list = first_line.split(",")
    # gene_list = [each.split("_")[-1] for each in gene_list]
    return gene_list


def extractSubNet(net, gene_list, TF_list = None):
    num_genes = len(gene_list)
    sub_net = np.zeros((num_genes, num_genes))
    sub_net[np.diag_indices(num_genes)] = 1.0
    sub_net = pd.DataFrame(sub_net, index=gene_list, columns=gene_list)
    # -----
    if TF_list is None:
        exist_gene_list = [each for each in gene_list if each in net.index.values]
    else:
        exist_gene_list = [each for each in gene_list if each in net.index.values and each in TF_list]
    print("# total genes={} | # existing genes={}".format(num_genes, len(exist_gene_list)))
    ref_sub_net = net.loc[exist_gene_list, exist_gene_list]
    sub_net.loc[exist_gene_list, exist_gene_list] = ref_sub_net
    return sub_net, ref_sub_net


if __name__ == '__main__':
    pbmc_sim_list = [
        "pbmc1-Drop", "pbmc2-Drop", "pbmc1-inDrops", "pbmc2-inDrops"
    ]
    cortex_sim_list = [
        "Cortex1-10xChromium", "Cortex2-10xChromium", "Cortex1-Smart_seq2", "Cortex2-Smart_seq2"
    ]
    tabula_data_list = ["skeletal_muscle_satellite_stem_cell", "T_cell", "type_B_pancreatic_cell"]
    # -----
    # Pathway Common gene similarity network
    # ref_net = loadPathwaySimilarityNet()
    # for d in pbmc_sim_list:
    #     gene_list = getGeneList(d)
    #     sub_net, _ = extractSubNet(ref_net, gene_list)
    #     sub_net.to_csv("./ChIP_seq_data/{}-Pathway_similarity-net_mat.csv".format(d))
    # # BioGPS gene similarity network
    # ref_net = loadBioGPSSimilarityNet()
    # for d in cortex_sim_list:
    #     gene_list = getGeneList(d)
    #     sub_net, _ = extractSubNet(ref_net, gene_list)
    #     sub_net.to_csv("./ChIP_seq_data/{}-BioGPS_similarity-net_mat.csv".format(d))
    # # ----
    # # Mouse TF ATLAS net
    # ref_net = loadMouseTFNet()
    # for d in cortex_sim_list:
    #     print("-" * 70)
    #     print(d)
    #     gene_list = getGeneList(d)
    #     _, ref_sub_net = extractSubNet(ref_net, gene_list)
    #     ref_sub_net.to_csv("./ChIP_seq_data/{}-mouse_TF_PCC-sub_net_mat.csv".format(d))
    # # ----
    # # Human TF similarity network
    # tf_list = loadTFList()
    # ref_net = loadPathwaySimilarityNet()
    # for d in pbmc_sim_list:
    #     print("-" * 70)
    #     print(d)
    #     gene_list = getGeneList(d)
    #     _, ref_sub_net = extractSubNet(ref_net, gene_list, tf_list)
    #     ref_sub_net.to_csv("./ChIP_seq_data/{}-human_TF_similarity-sub_net_mat.csv".format(d))
    # # ----
    # # Mouse TF ATLAS net (for Tabula Muris three cell types)
    # ref_net = loadMouseTFNet()
    # for d in tabula_data_list:
    #     print("-" * 70)
    #     print(d)
    #     gene_list = getTabulaGeneList(d)
    #     _, ref_sub_net = extractSubNet(ref_net, gene_list)
    #     ref_sub_net.to_csv("./ChIP_seq_data/{}-mouse_TF_PCC-sub_net_mat.csv".format(d))
