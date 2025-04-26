# Source Code Directory Explanation for Zebrafish Paper

This directory contains the core source code files of the Zebrafish project. Each file is responsible for different functional modules. Below is a detailed introduction to each file.

## File List and Functions

### 1. `GenerateMatrix.py`
- **Function**: Used to generate different types of `AnnData` objects, process Cyclone-related isoform count data, PSI data, aggregate isoform expression data to the gene level, and generate isoform score data.
- **Main Functions**:
  - `generate_Iso_adata`: Process Cyclone isoform data and save it as an `AnnData` object.
  - `generate_PSI_adata`: Read Cyclone PSI data and save it as an `AnnData` object.
  - `generate_Gene_adata`: Aggregate isoform expression data to the gene level and save it as an `AnnData` object.
  - `generate_IF_adata`: Generate an `AnnData` object containing isoform score data.

### 2. `MetaCell.py`
- **Function**: Automatically perform metacell processing on single-cell RNA-seq data, generate a metacell model based on clustering results, and save it.
- **Main Functions**:
  - `metacell_auto`: Automatically process single-cell RNA-seq data, generate a metacell model, and save the results.

### 3. `SpliceActivate.py`
- **Function**: Process gene events at specific stages, filter gene events that meet the criteria, and aggregate gene information.
- **Main Functions**:
  - `process_gene_events_for_stage`: Process gene events at a specific stage in the `AnnData` object and return a `DataFrame` containing gene information.

### 4. `TrajectoryKNN.py`
- **Function**: Use the K-Nearest Neighbors (KNN) algorithm for trajectory analysis and generate Sankey diagram data and charts.
- **Main Functions**:
  - `generate_labels`: Generate labels based on the replaced indices and distances.
  - `most_frequent_with_positions`: Return a dictionary of the most frequent elements and their positions in the input `Series`.
  - `traject_knn`: Perform trajectory analysis using the KNN algorithm.
  - `generate_sankey_data`: Generate links and nodes for the Sankey diagram based on clustering data.
  - `plot_sankey_chart`: Draw a Sankey diagram based on the provided nodes and links.

### 5. `TranscriptEntropy.py`
- **Function**: Calculate the isoform entropy score for a given gene list.
- **Main Functions**:
  - `_compute_score`: Calculate the transcript score based on the Kullback-Leibler divergence.
  - `isoform_entropy_score`: Calculate the isoform entropy score for a given gene list.

## Dependencies
The project uses several Python libraries, mainly including:
- `numpy`
- `pandas`
- `anndata`
- `scanpy`
- `omicverse`


## Usage
You can call the functions in each file according to your specific needs. For example, to generate an `AnnData` object for isoform data, you can call the `generate_Iso_adata` function in `GenerateMatrix.py` as follows:

```python
from GenerateMatrix import generate_Iso_adata

data_path = "your_isoform_count_data.csv"
adata = generate_Iso_adata(data_path)