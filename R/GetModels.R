#' Function that gives models trained to differentiate all pairwise matchups.
#'
#' @param marker_genes List with marker genes returned by the FindMarkerGenes function.
#' @param ref_bpcells BPCells obj with reference dataset GE values
#' @param ref_metadata Dataframe with metadata on each cell in reference dataset.
#' @param tree Tree structure in treedata format.
#' @param metadata_cluster_column Metadata celltype label column.
#' @param metadata_cell_label_column Metadata cell ID column
#' @param n_cells_sampled Number of cells used in pairwise model determination for each class.
#' @param models_to_include Optional vector which provides the names of models to include. If using
#' this parameter, include a subset of the following (make sure the names match or it won't work):
#' "linear_svm", "polynomial_svm", "naive_bayes", "ridge", "lasso", "elastic_net",
#'  "linear_da", "knn", "rf", "quadratic_da"
#' @param npcs Optional parameter giving number of PCs to use in model creation.
#'
#' @return List of models that differentiates each pairwise matchup.
#' @export
#'
#' @examples
#' GetModels(marker_genes, ref_bpcells, ref_metadata, tree,
#' metadata_cluster_column = "cluster_label",
#' metadata_cell_label_column = "cell_label",
#' n_cells_sampled = 500, models_to_include = NULL, npcs = 5)
GetModels <- function(marker_genes, ref_bpcells, ref_metadata, tree, metadata_cluster_column = "cluster_label", metadata_cell_label_column = "cell_label", n_cells_sampled = 500, models_to_include = NULL, npcs = 5) {
  #unit test: ref_metadata is a dataframe
  testthat::test_that("ref_metadata is a dataframe (not a tibble)", {
    # Assuming ref_metadata is defined in your environment
    testthat::expect_true(is.data.frame(ref_metadata) & ! dplyr::is.tbl(ref_metadata))
  })

  #1) Normalize reference atlas.
  # Normalize by reads-per-cell
  ref_bpcells <- BPCells::multiply_cols(ref_bpcells, 1/Matrix::colSums(ref_bpcells))
  # Log normalization
  ref_bpcells <- log1p(ref_bpcells * 10000) # Log normalization
  #0) Create list with same overall structure in terms of names and matchups as marker genes. Just set value as NA for now.
  model_list <- marker_genes
  for (i in 1:length(model_list)) {#iterate over parent nodes
    for (j in 1:length(model_list[[i]])) {#iterate over pairwise matchup
      model_list[[i]][[j]] <- NA
    }
  }

  #get list of tipnodes
  tipnodes <- tidytree::tip.label(tree)
  ref_bpcells %<>% BPCells::t()
  for (i in 1:length(marker_genes)) {  #Iterate through parent nodes
    for (j in 1:length(marker_genes[[i]])) { ##2) Iterate through matchups
      print(paste0("Model: ", as.character(i * length(marker_genes) + j - 1)))
      #3) Subset seurat object to only have cells pertaining to this matchup
      node1 <- marker_genes[[i]][[j]]$compared_nodes[1]
      node2 <- marker_genes[[i]][[j]]$compared_nodes[2]
      if(node1 %in% tipnodes) {
        node1_tip_nodes <- node1
      } else{
        node1_tip_nodes <- tree %>% tidytree::offspring(node1, type = "tips") %>% tidytree::nodelab(tree,.)
      }
      if(node2 %in% tipnodes) {
        node2_tip_nodes <- node2
      } else {
        node2_tip_nodes <- tree %>%  tidytree::offspring(node2, type = "tips") %>%  tidytree::nodelab(tree,.)
      }

      #if on tip node
      cells_node_1 <- ref_metadata[ref_metadata[,metadata_cluster_column] %in% c(node1_tip_nodes),metadata_cell_label_column]
      cells_node_2 <- ref_metadata[ref_metadata[,metadata_cluster_column] %in% c(node2_tip_nodes),metadata_cell_label_column]

      matchup_marker_genes <- marker_genes[[i]][[j]]$marker_genes
      #subset to only have particular genes and cells. Cells are ordered on whether they're from cell 1 or cell 2.
      subset_dataset <- ref_bpcells[c(cells_node_1, cells_node_2), matchup_marker_genes]
      #get average expression and variance of each gene in log normalized space
      n_threads <- 1
      gene_level_stats <- BPCells::matrix_stats(subset_dataset, col_stats = "variance", threads = n_threads)$col_stats
      avg_log_exp <-  gene_level_stats["mean",]
      #get stdev of each gene
      stdev <- sqrt(gene_level_stats["variance",])
      #z-score dataset
      subset_dataset <- subset_dataset %>% BPCells::add_cols(-avg_log_exp) %>% BPCells::multiply_cols(1/stdev)
      #sample n_cells_sampled # of cells from each node
      if(length(cells_node_1) > n_cells_sampled) {
        cells_node_1 <- cells_node_1 %>% sample(n_cells_sampled)
      }
      if(length(cells_node_2) > n_cells_sampled) {
        cells_node_2 <- cells_node_2 %>% sample(n_cells_sampled)
      }
      subset_dataset <- subset_dataset[c(cells_node_1,cells_node_2),]


      #write to memory to make it easier to load datasets
      #perform pca
      marker_ge_pca <- prcomp(subset_dataset, center = F, rank = npcs)

      #Unit test 1:
      if (i ==1 & j ==1) {
        testthat::test_that("Same Number rows in pca-transformed matrix as number of cells in each class for 1st matchup in 1st class", {
          testthat::expect_equal(marker_ge_pca$x %>% nrow(), length(cells_node_1) + length(cells_node_2))
        })
      }

      #save pca loadings (contribution of each variable to each pc) to variable, to use for later prediction
      #get labels for each cell in the training dataset in order of how it appears
      #reorder to get all node1 cells then node2 cells
      celltype_labels <- c(rep(node1, length(cells_node_1)), rep(node2, length(cells_node_2))) %>% as.factor()
      if(!is.null(models_to_include)) {
        classification_models <- CreateAllModels(marker_ge_pca$x, celltype_labels, models_to_include)
      } else {
        classification_models <- CreateAllModels(marker_ge_pca$x, celltype_labels)
      }


      model_list[[i]][[j]] <- model_list[[i]][[j]] %>% as.list()
      model_list[[i]][[j]][["Models"]] <- classification_models
      model_list[[i]][[j]][["avg_log_exp"]] <- avg_log_exp
      model_list[[i]][[j]][["stdev"]] <- stdev
      model_list[[i]][[j]][["pc_loadings"]] <- marker_ge_pca$rotation %>% t()
      model_list[[i]][[j]][["tip_labels"]] <- marker_genes[[i]][[j]][["tip_labels"]]
      model_list[[i]][[j]][["compared_nodes"]] <- marker_genes[[i]][[j]][["compared_nodes"]]

      #remove model list element without name in model list, just an artifact of how I created the list in the beginning of the function.
      model_list[[i]][[j]] <- model_list[[i]][[j]][model_list[[i]][[j]] %>% names() != ""]
    }
  }

  #7) Return list of pairwise models.
  model_list

}
