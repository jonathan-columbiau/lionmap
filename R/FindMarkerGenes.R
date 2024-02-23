#' Finds marker genes at each hierarchical level specified by the tree, using the GE matrix
#' provided in the ref_bpcells parameter, and the celltype labels provided in the input to ref_metadata. It
#' identifies marker genes by using a function provided in the BPCells package,
#' marker_features, which finds genes using the Wilcoxon test.
#'
#' @param ref_bpcells A GE reference dataset in BPCells format.
#' @param ref_metadata A dataframe with metadata that includes a column providing
#' celltype labels used for classification and a column providing cell ids.
#' @param tree A tree structure (in treedata format) to find marker genes for. Will find marker genes
#' that distinguish pairs of classes at each level of the hierarchy.
#' @param n_genes Number of marker genes, per pairwise class, you want to find.
#' @param metadata_cluster_column The name of the column in the metadata giving the celltype labels.
#' @param metadata_cell_label_column The name of the column in the metadata giving the cell IDs.
#' @param n_cells_sampled Number of cells per class used to find marker genes.
#'
#' @return A list providing marker genes that distinguish each pairwise combination of celltypes, at each
#' level of the hierarchy in the tree you provided.
#' @export
#'
#' @examples
#' FindMarkerGenes(ref_bpcells, ref_metadata)
FindMarkerGenes = function(ref_bpcells, ref_metadata, tree, n_genes = 50, metadata_cluster_column = "cell_type", metadata_cell_label_column = "cellid",n_cells_sampled = 500) {

  #unit test: ref_metadata is a dataframe
  testthat::test_that("ref_metadata is a dataframe (not a tibble)", {
    # Assuming ref_metadata is defined in your environment
    testthat::expect_true(is.data.frame(ref_metadata) & ! dplyr::is.tbl(ref_metadata))
  })
  testthat::test_that("The cluster column is in the metadata", {
    # Assuming ref_metadata is defined in your environment
    testthat::expect_true(metadata_cluster_column %in% colnames(ref_metadata))
  })

  testthat::test_that("The cell label column is in the metadata", {
    # Assuming ref_metadata is defined in your environment
    testthat::expect_true(metadata_cell_label_column %in% colnames(ref_metadata))
  })

  testthat::test_that("None of the cell types have any spaces in their names (one of the requirements of using Lionmap) ", {
    # Assuming ref_metadata is defined in your environment
    testthat::expect_true(all(!grepl(" ",unique(ref_metadata[,metadata_cell_label_column]))))
  })



  #Unit test 1: ref_bpcells is a bpcells object - else throw error
  testthat::test_that("ref_bpcells param is a bpcells object", {
    testthat::expect_equal(class(ref_bpcells) %>% attr("package"),"BPCells")
  })
  #1) Normalize reference atlas.
  # Normalize by reads-per-cell
  ref_bpcells <- BPCells::multiply_cols(ref_bpcells, 1/Matrix::colSums(ref_bpcells))
  # Log normalization
  ref_bpcells <-log1p(ref_bpcells * 10000) # Log normalization
  marker_genes <- vector(mode = "list")
  internal_nodes <- tree@phylo$node.label
  direct_child_nodes <- vector(mode = "list", length = length(internal_nodes))
  for (i in 1:length(internal_nodes)) {
    child_node_number_ids <- treeio::child(tree, internal_nodes[i])
    child_node_labels <- treeio::nodelab(tree, id = child_node_number_ids)
    direct_child_nodes[[i]] <- child_node_labels
  }
  names(direct_child_nodes) <- internal_nodes
  #remove internal nodes with only one child node (nothing to compare against)
  direct_child_nodes <- direct_child_nodes[lapply(direct_child_nodes,length)!=1]


  child_node_labels <- direct_child_nodes %>% unlist(use.names = F)
  descendant_tip_nodes <- vector(mode = "list", length = length(child_node_labels))
  names(descendant_tip_nodes) <- child_node_labels

  #add functionality for same-level classification only (tree structure is useless in this case)
  if(length(internal_nodes) == 1) {
    for (i in 1:length(child_node_labels)) {
      descendant_tip_nodes[[i]] <- child_node_labels[i]
    }

  } else {
    for (i in 1:length(child_node_labels)) {
      descendant_tip_nodes[[i]] <- treeio::offspring(tree,child_node_labels[i], type = "tips") %>% nodelab(tree,.)
    }
    #remove tip nodes from the above list of lists (they don't have any children nodes so their positions will have a length of 0)
    descendant_tip_nodes <- descendant_tip_nodes[lapply(descendant_tip_nodes,length)>0]
    specified_tip_nodes <- descendant_tip_nodes %>% unlist(use.names = F) %>% unique()

  }

  for (i in 1:length(direct_child_nodes)) { #iterate over each parent node
    specified_ancestor_node <- names(direct_child_nodes)[i]
    direct_children_of_specified_ancestor_nodes_vector <- direct_child_nodes[[specified_ancestor_node]]
    child_node_round_robin_matchups <- pairwise_combinations(direct_children_of_specified_ancestor_nodes_vector)
    list_with_matchups <- vector(mode = "list")
    for (j in 1:length(child_node_round_robin_matchups)) {
      node1 <- child_node_round_robin_matchups[[j]]$cluster1
      node2 <- child_node_round_robin_matchups[[j]]$cluster2
      compared_nodes <- c(node1, node2)
      if (node1 %in% names(descendant_tip_nodes)) {
        node1_tip_nodes <- descendant_tip_nodes[[node1]]
      } else {
        node1_tip_nodes <- node1
      }
      if (node2 %in% names(descendant_tip_nodes)) {
        node2_tip_nodes <- descendant_tip_nodes[[node2]]
      } else {
        node2_tip_nodes <- node2
      }

      #name the matchup of interest
      list_with_matchups[[j]] <- list("compared_nodes" = compared_nodes)
      names(list_with_matchups)[j] <- paste0(node1, " vs ", node2)
      #match colnames with class
      cells_node_1 <- ref_metadata[ref_metadata[,metadata_cluster_column] %in% c(node1_tip_nodes),metadata_cell_label_column]
      cells_node_2 <- ref_metadata[ref_metadata[,metadata_cluster_column] %in% c(node2_tip_nodes),metadata_cell_label_column]
      #sample cells
      if(length(cells_node_1) > n_cells_sampled) {
        cells_node_1 <- cells_node_1 %>% sample(n_cells_sampled)
      }
      if(length(cells_node_2) > n_cells_sampled) {
        cells_node_2 <- cells_node_2 %>% sample(n_cells_sampled)
      }

      subset_atlas <-ref_bpcells[, c(cells_node_1, cells_node_2)]
      celltype_labels <- c(rep(node1, length(cells_node_1)), rep(node2, length(cells_node_2))) %>% as.factor()
      pairwise_markers <- BPCells::marker_features(subset_atlas, celltype_labels, method = "wilcoxon")
      #remove genes with less than .5 logCPM in either class
      pairwise_markers %<>% dplyr::filter(foreground_mean > .5 |background_mean > .5) %>% dplyr::select(-background) %>% dplyr::distinct(feature, .keep_all = TRUE) %>% dplyr::mutate(log2_fc = log2(foreground_mean/background_mean))

      #get log2fc, and select the top marker genes with the highest abs value log2fc
      pairwise_markers %<>% dplyr::mutate(abs_log2_fc = log2(foreground_mean/background_mean) %>% abs()) %>% dplyr::arrange(abs_log2_fc) %>% dplyr::slice_max(abs_log2_fc, n = n_genes) %>% dplyr::pull(feature)
      list_with_matchups[[j]] <- c(list_with_matchups[[j]],list(marker_genes = pairwise_markers))

    }
    marker_genes[[i]] <-  list_with_matchups

  }

  #set naming and return list
  names(marker_genes) <- names(direct_child_nodes)
  marker_genes

}
