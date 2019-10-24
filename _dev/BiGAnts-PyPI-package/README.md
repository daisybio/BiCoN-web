# BiGAnts: network-constrained biclustering of patients and multi-omics data
PyPI package for conjoint clustering of networks and omics data

An application example is given in the file script_main.py in the project's [GitHub](https://github.com/biomedbigdata/BiGAnts-PyPI-package).

To install the package please run:
`pip install bigants` 

## Data input

The algorithm needs as an input one CSV matrix with gene expression/methylation/any other numerical data and one TSV file with a network.

### Numerical data

Numerical data is accepted in the following format:
- genes as rows.
- patients as columns.
- first column - genes IDs (can be any IDs).

For instance:

|   | Unnamed: 0 | GSM748056 | GSM748059 | ... | GSM748278 | GSM748279 | GSM1465989 |
|---|------------|-----------|-----------|-----|-----------|-----------|------------|
| 0 | 1454       | 0.053769  | 0.117412  | ... | -0.392363 | -1.870838 | -1.432554  |
| 1 | 201931     | -0.618279 | 0.278637  | ... | 0.803541  | -0.514947 | 2.361925   |
| 2 | 8761       | 0.215820  | -0.343865 | ... | 0.700430  | 0.073281  | -0.977656  |
| 3 | 2703       | -0.504701 | 1.295049  | ... | 1.861972  | 0.601808  | 0.191013   |
| 4 | 26207      | -0.626415 | -0.646977 | ... | 2.331724  | 2.339122  | -0.100924  |

There are 2 examples of gene expression datasets that can be placed in the "data" folder
- GSE30219 - a Non-Small Cell Lung Cancer dataset from GEO for patients with either adenocarcinoma or squamous cell carcinoma. 
- TCGA pan-cancer dataset with patients that have luminal or basal breast cancer.
Both can be found [here](https://drive.google.com/drive/folders/1J0XRrklwcV_Cgy_9Ay_6yJrN_x28Cosk?usp=sharing)

### Network

An interaction network should be present as a TSV table with two columns that represent two interacting genes. **Without a header!**

For instance:

|   | 6416 | 2318 |
|---|------|------|
| 0 | 6416 | 5371 |
| 1 | 6416 | 351  |
| 2 | 6416 | 409  |
| 3 | 6416 | 5932 |
| 4 | 6416 | 1956 |

In the data folder (on the [GitHub](https://github.com/biomedbigdata/BiGAnts-PyPI-package) page of the project) there is an example of a PPI network from Bioigrid with experimentally validated interactions.

## Functions

1. bigants.**data_preprocessing**(path_expr, path_net, log2 = False, size = 2000)

Parameters:

- path_to_expr: *string*, path to the numerical data 
- path_to_net: *string*, path to the network file
- log2: *bool, (default = False)*, indicates if log2 transformation should be applied to the data 
- size: *int, optional (default = 2000)* determines the number of genes that should be pre-selected by variance for the analysis. Shouldn't be higher than 5000.

Returns:

- GE: *pandas data frame*, processed expression data
- G: *networkX graph*, processed network data
- labels: *dict*, for mapping between real genes/patients IDs and the internal ones
- rev_labels: *dict*, additional dictionary for mapping between real genes/patients IDs and the internal ones

2. bigants.**BiGAnts**(GE,G,L_g_min,L_g_max) creates a model for the given data:

Parameters:

- GE: *pandas dataframe*, processed expression data
- G: *networkX graph*, processed network data
- L_g_min: *int*, minimal solution subnetwork size
- L_g_max: *int*, maximal solution subnetwork size

Methods:

bigants.BiGAnts.**run**(self, n_proc = 1, K = 20, evaporation = 0.5, show_plot = False)

- K: *int, default = 20*, number of ants. Fewer ants - less space exploration. Usually set between 20 and 50      
- n_proc: *int, default = 1*, number of processes that should be used
- evaporation, *float, default = 0.5*, the rate at which pheromone evaporates
- show_plot: *bool, default = False*, set true if convergence plots should be during the analysis
