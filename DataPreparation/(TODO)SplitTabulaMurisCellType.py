import scanpy
import numpy as np
import pandas as pd
import os
import anndata


def readTabulaData(
        path : str,
        must_include : str = "counts",
        num_tasks : int = None):
    """
    Read the Tabula Muris dataset.

    Args:
        path (str): Path of the directory containing original count data.
        must_include (str): The string that must be included in the file name.
                            This is used to exclude useless files. (default: "counts")
        num_tasks (int): The number of tasks to be read.
                         None means read data of all tasks. (default: None)

    Return:
        Dict[str, anndata.AnnData]: Data of each task.

    Exceptions:
        ValueError("The number of tasks to be read should not exceed the total task number!")
    """
    dir_list = os.listdir(path)
    dir_list = [each for each in dir_list if (must_include in each) and (".csv" in each)]
    print("The number of tasks = {}".format(len(dir_list)))
    if num_tasks is not None and isinstance(num_tasks, int):
        if num_tasks > len(dir_list):
            raise ValueError("The number of tasks to be read should not exceed the total task number!")
        dir_list = dir_list[:num_tasks]
        print("Read only the first {} tasks.".format(num_tasks))
    task_data_dict = {}
    for each in dir_list:
        task_data = scanpy.read_csv("{}/{}".format(path, each), first_column_names=True)
        task_data = task_data.T
        task_name = each.split("-")[0]
        task_data_dict[task_name] = task_data.copy()
        print("| Task : {} | Num of genes = {}, \t Num of cells = {}".format(
            task_name, task_data.shape[1], task_data.shape[0]
        ))
    return anndata.concat(list(task_data_dict.values()))


def readAnnotation(file_path):
    annot_data = pd.read_csv(file_path)
    annot_data = annot_data[["cell", "cell_ontology_class", "tissue"]].reset_index(drop=True)
    return annot_data


def getCell(annot_data, cell_type):
    if cell_type == "T cell":
        tmp_annot_data = annot_data[annot_data.cell_ontology_class == cell_type]
        tmp_tissue = tmp_annot_data.tissue.unique()
    elif cell_type == "skeletal muscle satellite stem cell":
        tmp_annot_data = annot_data[annot_data.cell_ontology_class == cell_type]
        tmp_tissue = tmp_annot_data.tissue.unique()
    elif cell_type == "type B pancreatic cell":
        tmp_annot_data = annot_data[annot_data.cell_ontology_class == cell_type]
        tmp_tissue = tmp_annot_data.tissue.unique()
    else:
        raise ValueError("Unused cell type {}!".format(cell_type))
    # -----
    data_path = "../../GeneNetwork-GGM/data/TabulaMuris/FACS/"
    relevant_task_data = anndata.concat([
        scanpy.read_csv("{}/{}-counts.csv".format(data_path, t), first_column_names=True).T
        for t in tmp_tissue
    ])
    # -----
    needed_cell = np.intersect1d(relevant_task_data.obs_names.values, tmp_annot_data.cell.values)
    cell_type_data = relevant_task_data.to_df().loc[needed_cell]
    return cell_type_data



def preprocessing(adata):
    scanpy.pp.normalize_per_cell(  # normalize with total UMI count per cell
        adata, key_n_counts='n_counts_all'
    )
    scanpy.pp.log1p(adata)
    scanpy.pp.highly_variable_genes(adata, n_top_genes=500)
    gene_idx = adata.var.highly_variable[adata.var.highly_variable.values == True].index.values
    adata = adata.to_df().loc[:, gene_idx]
    return adata


if __name__ == '__main__':
    # # Split by cell type
    annot_data = readAnnotation("../../GeneNetwork-GGM/data/TabulaMuris/description/annotations_facs.csv")
    cell_type_label = ["T cell", "skeletal muscle satellite stem cell", "type B pancreatic cell"]
    # cell_type_label = ["skeletal muscle satellite stem cell", "type B pancreatic cell"]
    for t in cell_type_label:
        tmp_t = t.replace(" ", "_")
        print("=" * 70)
        print("[ {} ]".format(t))
        cell_type_data = getCell(annot_data, cell_type=t)
        cell_type_data.to_csv("./Tabula_Muris/raw/{}.csv".format(tmp_t))
        print("Data shape (before): ", cell_type_data.shape)
        cell_type_data = preprocessing(anndata.AnnData(cell_type_data))
        cell_type_data.to_csv("./Tabula_Muris/500hvg/{}.csv".format(tmp_t))
        print("Data shape (after): ", cell_type_data.shape)
