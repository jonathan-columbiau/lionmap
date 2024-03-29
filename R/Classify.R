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
#' @importFrom tidytree rootnode
#' @importFrom tidytree nodelab
#' @importFrom tidytree offspring
#' @importFrom tidytree rootnode
#' @importFrom tidytree child
#' @importFrom purrr set_names
#' @importFrom purrr map2
#' @importFrom magrittr set_colnames
#' @importFrom plyr rbind.fill.matrix
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by
#' @importFrom dplyr count
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @importFrom dplyr n
#' @importFrom dplyr pull
#' @importFrom testthat test_that
#' @importFrom testthat expect_contains
#' @importFrom testthat expect_equal

#'
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
Classify <- function (bpcells_query, models, tree_struc, prop_max_threshold = 0.66)
{
  rootnode <- tree_struc@phylo %>% rootnode() %>% nodelab(tree_struc,
                                                          .)
  all_nodes <- offspring(tree_struc, rootnode) %>% nodelab(tree_struc,
                                                           .) %>% c(rootnode, .)
  tip_nodes = tidytree::offspring(tree_struc, rootnode, type = "tips") %>%
    nodelab(tree_struc, .)
  tipnodes = tip_nodes
  internal_nodes = all_nodes[!all_nodes %in% tip_nodes]
  bpcells_query <- BPCells::multiply_cols(bpcells_query, 1/Matrix::colSums(bpcells_query))
  bpcells_query = bpcells_query * 10000
  bpcells_query = log1p(bpcells_query)
  internal_node_assignment <- vector(mode = "list", length = length(internal_nodes)) %>%
    purrr::set_names(internal_nodes)
  query_cells <- bpcells_query
  internal_node_assignment[[rootnode]] <- colnames(query_cells)
  redundant_internal_nodes = child(tree_struc, internal_nodes) %>%
    lapply(function(x) {
      length(x) == 1
    }) %>% unlist() %>% .[. == T] %>% names()
  final_classifications <- vector(mode = "character")
  for (j in 1:length(internal_nodes)) {
    print(j)
    node <- internal_nodes[j]
    cells <- internal_node_assignment[[node]]
    if (length(cells) == 0) {
      next
    }
    if (node %in% redundant_internal_nodes) {
      child_node = child(tree_struc, node) %>% nodelab(tree_struc,
                                                       .)
      if (!child_node %in% tipnodes) {
        internal_node_assignment[[child_node]] <- cells
      }
      else {
        these_classifications = rep(child_node, length(cells))
        names(these_classifications) = cells
        final_classifications <- c(final_classifications,
                                   these_classifications)
      }
      next
    }
    res_list <- vector(mode = "list", length = length(models[[node]]))
    for (i in 1:length(models[[node]])) {
      first_lev_avg_counts <- models[[node]][[i]]$avg_log_exp
      first_lev_std_counts <- models[[node]][[i]]$stdev
      first_lev_markers <- models[[node]][[i]]$pc_loadings %>%
        colnames()
      first_lev_bpcells <- query_cells[first_lev_markers,
                                       cells] %>% BPCells::t()
      first_lev_bpcells <- first_lev_bpcells %>% BPCells::add_cols(-first_lev_avg_counts) %>%
        BPCells::multiply_cols(1/first_lev_std_counts)
      first_lev_pc_loadings <- models[[node]][[i]]$pc_loadings[,
                                                               first_lev_markers] %>% t()
      first_lev_bpcells <- first_lev_bpcells %*% first_lev_pc_loadings
      nonsparse_mat <- first_lev_bpcells %>% as.matrix()
      first_lev_models <- models[[node]][[i]]$Models
      res_list[[i]] <- map2(first_lev_models, first_lev_models %>%
                              names(), predict_models, nonsparse_mat) %>% as.data.frame() %>%
        set_colnames(paste0(colnames(.), "_", i)) %>%
        t()
    }
    count_threshold <- (child(tree_struc, node) %>% length() -
                          1) * length(models[[1]][[1]][["Models"]]) * prop_max_threshold
    obs_above_threshold <- res_list %>% rbind.fill.matrix(res_list) %>%
      as.data.frame() %>% pivot_longer(everything(), names_to = "obs",
                                       values_to = "class") %>% group_by(obs) %>% count(class,
                                                                                        name = "count") %>% group_by(obs) %>% filter(count ==
                                                                                                                                       max(count))
    tied_obs <- obs_above_threshold %>% group_by(obs) %>%
      summarise(n = n()) %>% filter(n > 1) %>% pull(obs)
    obs_above_threshold <- obs_above_threshold %>% filter(!obs %in%
                                                            tied_obs) %>% filter(count >= count_threshold) %>%
      pull(class, name = obs)
    tip_cells <- obs_above_threshold[obs_above_threshold %in%
                                       tip_nodes]
    final_classifications <- final_classifications %>% append(tip_cells)
    stuck_cells <- rownames(first_lev_bpcells)[!rownames(first_lev_bpcells) %in%
                                                 names(final_classifications) & !rownames(first_lev_bpcells) %in%
                                                 names(obs_above_threshold)] %>% set_names(., .)
    stuck_cells <- rep(node, length(stuck_cells))
    names(stuck_cells) <- rownames(first_lev_bpcells)[!rownames(first_lev_bpcells) %in%
                                                        names(final_classifications) & !rownames(first_lev_bpcells) %in%
                                                        names(obs_above_threshold)]
    final_classifications <- final_classifications %>% append(stuck_cells)
    if (obs_above_threshold[!obs_above_threshold %in% tipnodes &
                            !obs_above_threshold %in% names(final_classifications)] %>%
        length() > 0) {
      obs_above_threshold <- obs_above_threshold[obs_above_threshold %in%
                                                   internal_nodes]
      obs_above_threshold <- split(obs_above_threshold,
                                   obs_above_threshold)
      test_that("all remaining cells assigned to internal nodes",
                {
                  expect_contains(internal_nodes, names(obs_above_threshold))
                })
      for (name in names(obs_above_threshold)) {
        internal_node_assignment[[name]] <- obs_above_threshold[[name]] %>%
          names()
      }
    }
  }
  test_that("expected number of elements returned", {
    expect_equal(length(final_classifications), ncol(bpcells_query))
  })
  final_classifications[match(colnames(bpcells_query), names(final_classifications))]
}
