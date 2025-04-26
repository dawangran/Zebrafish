# -*- coding: utf-8 -*-
"""
@File    :   SpliceActivate.py
@Time    :   2024/07/01 
@Author  :   Dawn
@Version :   1.0
@Desc    :   Splice activate for gene events
"""


import scanpy as sc
import pandas as pd

def process_gene_events_for_stage(adata_psi, stage_column, stage_value, min_event_count=10):
    """
    Process gene events for a specific stage in the adata object.

    Parameters:
        adata_psi (AnnData): Input AnnData object containing PSI (Percent Spliced In) data.
        stage_column (str): The column name in adata.obs that contains the stage information.
        stage_value (str): The specific stage value to process (e.g., '128C', '512C', '3h', etc.).
        min_event_count (int): Minimum number of events required per gene to be considered, default is 10.

    Returns:
        pd.DataFrame: A DataFrame containing gene information for the specified stage.
    """
    # Subset the adata object by the specified stage value
    sub_adata = adata_psi[adata_psi.obs[stage_column] == stage_value]
    
    # Convert to DataFrame
    df = sub_adata.to_df()
    
    # Count events with values between 0 and 1
    count_per_column = (df > 0) & (df < 1)
    counts = count_per_column.sum()
    
    # Select events with more than min_event_count occurrences
    psi_event = list(counts[counts > min_event_count].index)
    
    # Subset the adata object by selected events
    sub_adata = sub_adata[:, psi_event]
    
    # Aggregate gene information
    gene_sub = sub_adata.var.groupby("gene_name").agg({"event": "count"})
    gene_sub['stage'] = stage_value
    
    return gene_sub
