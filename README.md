<!-- README.md is generated from README.Rmd. Please edit that file -->

# lionmap

Hello! Welcome to the Lionmap package :)

This package is meant to provide efficient and robust cell type classification while being easy to use.

Optional features to increase accuracy include the use of a custom hierarchy that describes the relationship between cell types, and user-defined confidence thresholds at each stage of the classification.

An example using a sample dataset where we run the entire pipeline (finding marker genes, creating models, and classifying cells) in a few lines of code is shown below.

## Example

``` r
library(lionmap)
#load reference and query datasets
data("hierarchy","train_ex_data_bpcells","test_ex_data_bpcells","train_ex_metadata")
#find marker genes 
marker_genes = FindMarkerGenes(train_ex_data_bpcells, train_ex_metadata, tree = hierarchy, metadata_cluster_column = "seurat_annotations", metadata_cell_label_column = "cell_label")
#create models
models = GetModels(marker_genes, train_ex_data_bpcells, train_ex_metadata, tree = hierarchy, metadata_cluster_column = "seurat_annotations", metadata_cell_label_column = "cell_label")
#classify cells using models
classifications = Classify(test_ex_data_bpcells, models, hierarchy)
```

## Installation

You can install the development version of lionmap from [GitHub](https://github.com/jonathan-columbiau/lionmap) with:

``` r
# install.packages("devtools")
devtools::install_github("jonathan-columbiau/lionmap")
```
