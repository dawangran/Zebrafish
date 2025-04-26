# -*- coding: utf-8 -*-
"""
@File    :   MetaCell.py
@Time    :   2024/05/21 
@Author  :   Dawn
@Version :   1.0
@Desc    :   Metacell
"""


import scanpy as sc
import omicverse as ov
import pickle
import os



def metacell_auto(input_file, output_dir, cluster_key="stage_cluster", min_cell_num=10, scaling_factor=4000):
    """
    Metacell auto process for single-cell RNA-seq data.

    Parameters:
        input_file (str): Path to the input h5ad file.
        output_dir (str): Path to the output directory.
        cluster_key (str): Column name used for grouping, default is "stage_cluster".
        min_cell_num (int): Minimum number of cells, default is 10.
        scaling_factor (int): Scaling factor used to calculate cell_num, default is 4000.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    with open(input_file, 'rb') as f:
        adata = pickle.load(f)
    
    # Add a cell_id column
    adata.obs['cell_id'] = adata.obs.index.to_list()

    # Calculate the number of cells and average gene counts for each cluster_key
    data = adata.obs.groupby(cluster_key).agg({"n_genes_by_counts": "mean", 'cell_id': 'count'})
    data['sum'] = data['n_genes_by_counts'] * data['cell_id']
    data['cell_num'] = data['sum'] / scaling_factor  # Use scaling_factor for scaling
    data['cell_num'] = data['cell_num'].astype(int)
    data.loc[data['cell_num'] < min_cell_num, "cell_num"] = min_cell_num

    # Create a dictionary mapping cluster_key to cell_num
    cluster_dict = {k: v for k, v in zip(data.index, data['cell_num'])}

    # Initialize variables
    names = locals()
    adata_list = []
    raw_adata_list = []

    # Iterate over each cluster, generate MetaCell, and save the model
    for i in cluster_dict:
        sub_adata = adata[adata.obs[cluster_key] == i]
        names['meta_obj_' + str(i)] = ov.single.MetaCell(sub_adata, use_rep='X_pca',
                                                         n_metacells=cluster_dict[i],
                                                         use_gpu=False)
        names['meta_obj_' + str(i)].initialize_archetypes()
        names['meta_obj_' + str(i)].train(min_iter=10, max_iter=30)
        names['meta_obj_' + str(i)].save(os.path.join(output_dir, f'model.log.auto.mean.2k.{i}.pkl'))

        # Generate the predicted adata object
        names['ad_' + str(i)] = names['meta_obj_' + str(i)].predicted(method='hard', summarize_layer='counts')
        names['ad_' + str(i)].obs[cluster_key] = str(i)
        adata_list.append(names['ad_' + str(i)])
        raw_adata_list.append(sub_adata)

    # Concatenate all adata objects
    adata_meta = sc.concat(adata_list)
    adata_meta.var = adata.var

    # Save the concatenated adata object
    adata_meta.write(os.path.join(output_dir, 'merged_adata_meta.h5ad'))

    print("Processing complete. Results have been saved to the specified directory.")
    return adata_meta