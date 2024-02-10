Hello! Welcome to the Lionmap package :)

This package is meant to provide efficient and robust cell type classification while being easy to use.

Optional features that might increase mapping accuracy include using a custom hierarchy that describes the relationship between cell types, and setting different confidence thresholds for a classification to proceed to the next level of the hierarchy.

An example using a sample dataset where we run the entire pipeline (finding marker genes, creating models, and classifying cells) in a few lines of code is shown below.

```{r, eval=FALSE}
#load reference and query datasets
data("hierarchy","train_ex_data_bpcells","test_ex_data_bpcells","train_ex_metadata")
#find marker genes 
marker_genes = FindMarkerGenes(train_ex_data_bpcells, train_ex_metadata, tree = hierarchy, metadata_cluster_column = "seurat_annotations", metadata_cell_label_column = "cell_label")
#create models
models = GetModels(marker_genes, train_ex_data_bpcells, train_ex_metadata, tree = hierarchy, metadata_cluster_column = "seurat_annotations", metadata_cell_label_column = "cell_label")
#classify cells using models
classifications = Classify(test_ex_data_bpcells, models, hierarchy)

```

The package only works using BPCells matrices as gene expression input. BPCells is a package that enables fast gene expression transformations and IO operations while requiring little memory on your computing environment. For more info on BPCells, check out the package [here](https://bnprks.github.io/BPCells/).
