#' Classification function that takes in input dataset (in BPCells format), alongside
#' created models and the tree structure and provides a vector of classifications.
#'
#' @param bpcells_query Query dataset already aligned to reference dataset used to find marker genes and create models.
#' @param models Models created through GetMarkerGenes function.
#' @param tree Tree used in model creation.
#' @param prop_max_threshold Proportion of evidence required for a test
#'
#' @return A vector providing classifications of cells in bpcells_query in the same order.
#' Performs specific unit testing on each run through.
#' Unit test 1: all remaining cells assigned to internal nodes
#' Unit test 2: expected number of elements returned
#' @export
#'
#' @examples
#' data("train_ex_data_bpcells")
#' data("train_ex_metadata")
#' data("test_ex_data_bpcells")
#' data("test_ex_metadata")
#' possible_cell_classes = train_ex_metadata$seurat_annotations %>% unique()
#' equal_tree = CreateEqualTree(cell_labels = possible_cell_classes)
#' marker_genes = FindMarkerGenes(ref_bpcells = train_ex_data_bpcells, ref_metadata = train_ex_metadata, tree = equal_tree, metadata_cluster_column = "seurat_annotations", metadata_cell_label_column = "cell_label")
#' models <- GetModels(marker_genes = marker_genes, ref_bpcells = train_ex_data_bpcells, ref_metadata = train_ex_metadata, tree = equal_tree, metadata_cluster_column = "seurat_annotations", metadata_cell_label_column = "cell_label")
#' query_classifications = Classify(bpcells_query = test_ex_data_bpcells,models = models,tree_struc = equal_tree)
#'
#'
#'
Classify <- function(bpcells_query, models, tree_struc, prop_max_threshold = .66) {
  # initial rootnode level - do all at once. All tree_struc have rootnodes, even if all clusters are at one level
  rootnode <- tree_struc %>%
    tidytree::rootnode() %>%
    tidytree::nodelab(tree_struc, .)
  tipnodes <- tidytree::nodelab(tree_struc, tidytree::offspring(tree_struc, rootnode, type = "tips"))
  # Normalize by reads-per-cell
  bpcells_query <- BPCells::multiply_cols(bpcells_query, 1 / Matrix::colSums(bpcells_query))
  # Log normalization
  bpcells_query <- log1p(bpcells_query * 10000) # Log normalization
  # get internal nodes in hierarchical order
  internal_nodes <- c(rootnode, tidytree::offspring(tree_struc, rootnode) %>% tidytree::nodelab(tree_struc, .) %>% .[-which(. %in% tidytree::nodelab(tree_struc, tidytree::offspring(tree_struc, rootnode, type = "tips")))])
  internal_node_assignment <- vector(mode = "list", length = length(internal_nodes)) %>% purrr::set_names(internal_nodes)
  query_cells <- bpcells_query
  internal_node_assignment[[rootnode]] <- colnames(query_cells)

  #add support for trees with redundant internal nodes
  redundant_internal_nodes = tidytree::child(tree_struc, internal_nodes) %>% lapply(function(x) {length(x) ==1}) %>% unlist() %>% .[. == T] %>% names()

  final_classifications <- vector(mode = "character")
  for (j in 1:length(internal_nodes)) { # iterate node
    print(j)
    ## track cells that don't go past internal nodes
    node <- internal_nodes[j]

    cells <- internal_node_assignment[[node]]

    #check if redundant internal node. If so, assign to offspring
    if (node %in% redundant_internal_nodes) {
      child_node = tidytree::child(tree_struc, node) %>% nodelab(tree_struc,.)
      #check if offspring is internal node or tip - since redundant, only has one element
      #by definition
      #if internal node assign to that in internal_node_assignment list
      if(!child_node %in% tipnodes) {
        internal_node_assignment[[child_node]] <- cells
      } else {#if tip assign to child in final_classifications vector
        these_classifications = rep(child_node, length(cells))
        names(these_classifications) = cells
        final_classifications <- c(final_classifications, these_classifications)
      }
      next
    }

    res_list <- vector(mode = "list", length = length(models[[node]]))

    ## final classification for cells that don't go past internal node = internal node
    for (i in 1:length(models[[node]])) { # iterate over models, classify
      first_lev_avg_counts <- models[[node]][[i]]$avg_log_exp # scale
      first_lev_std_counts <- models[[node]][[i]]$stdev
      first_lev_markers <- models[[node]][[i]]$pc_loadings %>% colnames()
      first_lev_bpcells <- query_cells[first_lev_markers, cells] %>% BPCells::t() # select markers and cells of this internal node level
      first_lev_bpcells <- first_lev_bpcells %>%
        BPCells::add_cols(-first_lev_avg_counts) %>%
        BPCells::multiply_cols(1 / first_lev_std_counts)
      first_lev_pc_loadings <- models[[node]][[i]]$pc_loadings[, first_lev_markers] %>% t()
      # transform data using matrix multiplication operator %*%
      first_lev_bpcells <- first_lev_bpcells %*% first_lev_pc_loadings
      nonsparse_mat <- first_lev_bpcells %>% as.matrix()
      # use models
      first_lev_models <- models[[node]][[i]]$Models
      res_list[[i]] <- purrr::map2(first_lev_models, first_lev_models %>% names(), predict_models, nonsparse_mat) %>%
        as.data.frame() %>%
        magrittr::set_colnames(paste0(colnames(.), "_", i)) %>%
        t()
    }

    count_threshold <- (tidytree::child(tree_struc,node) %>% length() - 1) * length(models[[1]][[1]][["Models"]]) * prop_max_threshold #2/3 of max score (num matchups for each class - 1)*num models per pairwise comparison

    obs_above_threshold <- res_list %>%
      plyr::rbind.fill.matrix(res_list) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(everything(), names_to = "obs", values_to = "class") %>%
      dplyr::group_by(obs) %>%
      dplyr::count(class, name = "count") %>%
      dplyr::group_by(obs) %>%
      dplyr::filter(count == max(count))
    #filter obs with mult max classes
    tied_obs <- obs_above_threshold %>% dplyr::group_by(obs) %>% dplyr::summarise(n = dplyr::n()) %>% dplyr::filter(n > 1) %>% dplyr::pull(obs)
    obs_above_threshold <- obs_above_threshold %>%
      dplyr::filter(!obs %in% tied_obs) %>%
      dplyr::filter(count >= count_threshold) %>%
      dplyr::pull(class, name = obs)
    ## assign "tip" cells to final classification
    tip_cells <- obs_above_threshold[obs_above_threshold %in% tipnodes]
    final_classifications <- final_classifications %>% append(tip_cells)
    ## assign "stuck" cells to final classification - ties/not threshold
    stuck_cells <- rownames(first_lev_bpcells)[!rownames(first_lev_bpcells) %in% names(final_classifications) & !rownames(first_lev_bpcells) %in% names(obs_above_threshold)] %>% purrr::set_names(., .)
    stuck_cells <- rep(node, length(stuck_cells))
    names(stuck_cells) <- rownames(first_lev_bpcells)[!rownames(first_lev_bpcells) %in% names(final_classifications) & !rownames(first_lev_bpcells) %in% names(obs_above_threshold)]
    final_classifications <- final_classifications %>% append(stuck_cells)
    ## update internal node assignment list
    if (obs_above_threshold[!obs_above_threshold %in% tipnodes & !obs_above_threshold %in% names(final_classifications)] %>% length() > 0) {
      obs_above_threshold <- obs_above_threshold[obs_above_threshold %in% internal_nodes]
      obs_above_threshold <- split(obs_above_threshold, obs_above_threshold)
      # unit test 1: all remaining cells assigned to internal nodes
      testthat::test_that("all remaining cells assigned to internal nodes", {
        testthat::expect_contains(internal_nodes, names(obs_above_threshold))
      })
      for (name in names(obs_above_threshold)) {
        internal_node_assignment[[name]] <- obs_above_threshold[[name]] %>% names()
      }
    }
  }
  testthat::test_that("expected number of elements returned", {
    testthat::expect_equal(length(final_classifications), ncol(bpcells_query))
  })
  #order match bpcells_query
  final_classifications[match(colnames(bpcells_query),names(final_classifications))]
}
