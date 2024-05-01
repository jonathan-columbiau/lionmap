#' @param ref_bpcells A GE reference dataset in BPCells format.
#' @param ref_metadata A dataframe with metadata that includes a column providing
#' celltype labels used for classification and a column providing cell ids.
#' @param tree A tree structure (in treedata format) to find marker genes for. Will find marker genes
#' that distinguish pairs of classes at each level of the hierarchy.
#' @param metadata_cluster_column The name of the column in the metadata giving the celltype labels.
#' @param metadata_cell_id_column The name of the column in the metadata giving the cell IDs.
#' @param min_n_cells_sampled Minimum  number of cells to find marker genes per class.
#' @param max_n_cells_sampled Max number of cells to find marker genes per class if greater than threshold.
#'
#'
#' @return A list providing marker genes that distinguish each pairwise
#' combination of celltypes, whether p-value significant or not, the top certain
#' number of marker genes arranged by log fold change. Given
#' at each level of the hierarchy in the tree you provided.
#'
#' @importFrom tidytree child
#' @importFrom tidytree nodelab
#' @importFrom testthat test_that
#' @importFrom testthat expect_true
#' @importFrom testthat expect_equal
#' @importFrom tidytree offspring
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr slice_max
#' @importFrom dplyr arrange
#' @importFrom dplyr is.tbl
#' @export
#'
#' @examples
#' data("train_ex_data_bpcells")
#' data("train_ex_metadata")
#' data("test_ex_data_bpcells")
#' data("test_ex_metadata")
#' possible_cell_classes = train_ex_metadata$seurat_annotations %>% unique()
#' equal_tree = CreateEqualTree(cell_labels = possible_cell_classes)
#' marker_genes = FindMarkerGenes(ref_bpcells = train_ex_data_bpcells, ref_metadata = train_ex_metadata, tree = equal_tree, metadata_cluster_column = "seurat_annotations", metadata_cell_id_column = "cell_label")
FindMarkerGenes = function(ref_bpcells, ref_metadata, tree, metadata_cluster_column = "cell_type", metadata_cell_id_column = "cellid",min_n_cells_sampled = 100, max_n_cells_sampled = 5000, n_marker_genes_per_comparison = 50) {

  #unit test: ref_metadata is a dataframe
  test_that("ref_metadata is a dataframe (not a tibble)", {
    # Assuming ref_metadata is defined in your environment
    expect_true(is.data.frame(ref_metadata) & ! is.tbl(ref_metadata))
  })
  test_that("The cluster column is in the metadata", {
    # Assuming ref_metadata is defined in your environment
    expect_true(metadata_cluster_column %in% colnames(ref_metadata))
  })

  test_that("The cell label column is in the metadata", {
    # Assuming ref_metadata is defined in your environment
    expect_true(metadata_cell_id_column %in% colnames(ref_metadata))
  })

  test_that("None of the cell types have any spaces in their names (one of the requirements of using Lionmap) ", {
    # Assuming ref_metadata is defined in your environment
    expect_true(all(!grepl(" ",unique(ref_metadata[,metadata_cell_id_column]))))
  })



  #Unit test 1: ref_bpcells is a bpcells object - else throw error
  test_that("ref_bpcells param is a bpcells object", {
    expect_equal(class(ref_bpcells) %>% attr("package"),"BPCells")
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
    child_node_number_ids <- child(tree, internal_nodes[i])
    child_node_labels <- nodelab(tree, id = child_node_number_ids)
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
      descendant_tip_nodes[[i]] <- offspring(tree,child_node_labels[i], type = "tips") %>% nodelab(tree,.)
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
      cells_node_1 <- ref_metadata[ref_metadata[,metadata_cluster_column] %in% c(node1_tip_nodes),metadata_cell_id_column]
      cells_node_2 <- ref_metadata[ref_metadata[,metadata_cluster_column] %in% c(node2_tip_nodes),metadata_cell_id_column]
      #amended for metric 2 - if not enough mapped clusters (at the most granular level) are present, skip this pairwise combo of finding marker genes
      if(length(cells_node_1) > min_n_cells_sampled) {
        cells_node_1 <- cells_node_1 %>% sample(min_n_cells_sampled)
      } else {
        list_with_matchups[[j]] <- c(list_with_matchups[[j]],list(marker_genes = "NA - Not enough cells of one group"))
        next
      }
      if(length(cells_node_2) > min_n_cells_sampled) {
        cells_node_2 <- cells_node_2 %>% sample(min_n_cells_sampled)
      } else {
        stop(paste0("Error: Not enough cells of one group in a matchup - either ",node1, "or ", node2))
        list_with_matchups[[j]] <- c(list_with_matchups[[j]],list(marker_genes = "NA - Not enough cells of one group"))
        next
      }

      #now check if > max_n_cells_sampled cells in each group, in which case sample up to max_n_cells_sampled
      if(length(cells_node_1) < max_n_cells_sampled | length(cells_node_2) < max_n_cells_sampled) {
        num_cells = min(length(cells_node_1),length(cells_node_2))
        cells_node_1 <- cells_node_1 %>% sample(num_cells)
        cells_node_2 <- cells_node_2 %>% sample(num_cells)
      } else {
        #sample max_n_cells_sampled cells if both have >= max_n_cells_sampled
        cells_node_1 <- cells_node_1 %>% sample(max_n_cells_sampled)
        cells_node_2 <- cells_node_2 %>% sample(max_n_cells_sampled)
      }

      #if greater than min (n_cells_sampled) and less than max (max_n_cells_sampled), then get the number of cells equal to the lower of the two node cells

      subset_atlas <-ref_bpcells[, c(cells_node_1, cells_node_2)]
      celltype_labels <- c(rep(node1, length(cells_node_1)), rep(node2, length(cells_node_2))) %>% as.factor()
      pairwise_markers <- BPCells::marker_features(subset_atlas, celltype_labels, method = "wilcoxon")

      if(typeof(pairwise_markers$feature) == "integer") {
        pairwise_markers$feature <- rownames(subset_atlas)[pairwise_markers$feature]
      }
      #for some reason this isn't returning gene names at high levels of rownames so changed if integer to character

      #add functionality for no markers - do nothing, write no marker genes
      if(nrow(pairwise_markers) == 0) {
        list_with_matchups[[j]] <- c(list_with_matchups[[j]],list(marker_genes = NA))
        next
      }

      #remove genes with less than .5 logCPM in either class
      pairwise_markers %<>% filter(foreground_mean > .5 |background_mean > .5) %>% dplyr::select(-background) %>% distinct(feature, .keep_all = TRUE) %>% mutate(log2_fc = log2(foreground_mean/background_mean))

      #get log2fc, and dplyr::select marker genes arranged by abs value log2fc
      pairwise_markers %<>% mutate(abs_log2_fc = log2(foreground_mean/background_mean) %>% abs()) %>% slice_max(n = n_marker_genes_per_comparison, order_by = abs_log2_fc) %>%  pull(feature)
      list_with_matchups[[j]] <- c(list_with_matchups[[j]],list(marker_genes = pairwise_markers %>% as.character()))
    }
    marker_genes[[i]] <-  list_with_matchups

  }

  #set naming and return list

  names(marker_genes) <- names(direct_child_nodes)
  marker_genes

}

