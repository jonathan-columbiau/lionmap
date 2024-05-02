#' Classification function that takes in input dataset (in BPCells format), alongside
#' created models and the tree structure and provides a vector of classifications.
#'
#' @param bpcells_query Query dataset already aligned to reference dataset used to find marker genes and create models.
#' @param models Models created through GetMarkerGenes function.
#' @param tree Tree used in model creation.
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
#' marker_genes = FindMarkerGenes(ref_bpcells = train_ex_data_bpcells, ref_metadata = train_ex_metadata, tree = equal_tree, metadata_cluster_column = "seurat_annotations", metadata_cell_id_column = "cell_label")
#' models <- GetModels(marker_genes = marker_genes, ref_bpcells = train_ex_data_bpcells, ref_metadata = train_ex_metadata, tree = equal_tree, metadata_cluster_column = "seurat_annotations", metadata_cell_label_column = "cell_label")
#' query_classifications = Classify_Return_All(bpcells_query = test_ex_data_bpcells,models = models,tree_struc = equal_tree)
#'
#'
#'change classify so it gives onfidence score at each level of the hierarchy
#'a cell is classified at.
#'Returns a dataframe with the following columns for each classification step
#'performed
#'
#'
#
#

Classify_Return_All <- function (bpcells_query, models, tree_struc)
{
  returned_df = matrix(nrow = 0, ncol = 6) %>% as.data.frame()
  colnames(returned_df) = c("cell_id","best_classification",
                            "count",
                            "max_possible_score",
                            "cur_internal_node",
                            "is_final_classification_for_cell")
  class(returned_df$cell_id) = "character"
  class(returned_df$best_classification) = "character"
  class(returned_df$count) = "numeric"
  class(returned_df$max_possible_score) = "numeric"
  class(returned_df$cur_internal_node) = "character"
  class(returned_df$is_final_classification_for_cell) = "logical"

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
    for (i in 1:length(models[[node]])) { #returns list of dfs with results of each pairwise comparison at an internal node
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
    max_possible_count_cur_internal_node <- (child(tree_struc, node) %>% length() -
                                               1) * length(models[[1]][[1]][["Models"]])
    best_classification_per_obs_with_counts <- res_list %>%
      rbind.fill.matrix(res_list) %>%
      as.data.frame() %>%
      pivot_longer(everything(), names_to = "obs",
                   values_to = "class") %>%
      group_by(obs) %>% count(class, name = "count") %>%
      group_by(obs) %>%
      dplyr::filter(count == max(count))


    best_classification_per_obs_with_counts %<>% rename(cell_id = obs, best_classification= class)
    best_classification_per_obs_with_counts$max_possible_score = max_possible_count_cur_internal_node
    best_classification_per_obs_with_counts$cur_internal_node = node

    tied_obs <- best_classification_per_obs_with_counts %>%
      group_by(cell_id) %>%
      summarise(n = n()) %>%
      filter(n > 1) %>%
      pull(cell_id)

    best_classification_per_obs_with_counts <- best_classification_per_obs_with_counts %>%
      filter(!cell_id %in% tied_obs) %>% mutate(confidence_score = count/max_possible_score)

    best_classification_per_obs_with_counts$is_final_classification_for_cell = F


    #write NA for and best_classification + confidence_score in case of ties
    #add these columns c("cell_id", "cur_internal_node",
    #"best_classification", "count", "max_possible_score",
    #"confidence_score","is_final_classification_for_cell")
    cur_internal_node = node

    #special to return_full_df = T
    if(length(tied_obs) > 0) {
      tied_obs_returned = data.frame("cell_id" = tied_obs,
                                     cur_internal_node = node,
                                     count = NA,
                                     best_classification = "NA/Tie",
                                     confidence_score = NA,
                                     max_possible_score = max_possible_count_cur_internal_node,
                                     is_final_classification_for_cell = TRUE)
        returned_df = bind_rows(returned_df, tied_obs_returned)
        stuck_cells = rep(node,length(tied_obs_returned$cell_id %>% unique())) #cells stuck at current internal node
        names(stuck_cells) = tied_obs_returned$cell_id %>% unique()
        final_classifications <- final_classifications %>% append(stuck_cells)
    }


    tip_cell_ids <- best_classification_per_obs_with_counts$cell_id[best_classification_per_obs_with_counts$best_classification %in%
                                                                   tip_nodes]
    tip_cell_classifications <- best_classification_per_obs_with_counts$best_classification[match(tip_cell_ids,best_classification_per_obs_with_counts$cell_id)]

    names(tip_cell_classifications) <- tip_cell_ids
    #special for returned_df
    if(length(tip_cell_ids) != 0) {
      tip_cell_obs_returned = best_classification_per_obs_with_counts %>% filter(cell_id %in% tip_cell_ids)
      tip_cell_obs_returned$is_final_classification_for_cell = T
      returned_df = bind_rows(returned_df, tip_cell_obs_returned)
      final_classifications <- final_classifications %>% append(tip_cell_classifications)
    }








    remaining_cells = best_classification_per_obs_with_counts %>% filter(!cell_id %in% tip_cell_ids)
    if (nrow(remaining_cells) > 0) {
      test_that("all remaining cells assigned to internal nodes",
                {
                  expect_contains(internal_nodes, remaining_cells$best_classification)
                })
      for (name in remaining_cells$best_classification %>% unique()) {
        #add to internal_node_assignment list
        internal_node_assignment[[name]] <- remaining_cells$cell_id
      }
      #special to returned_df = T
      returned_df = bind_rows(returned_df, remaining_cells)
    }
  }
  test_that("expected number of elements (all elements) in final_classifications vector", {
    expect_equal(length(final_classifications), ncol(bpcells_query))
  })
  #special to returned_df = T
  return(returned_df)
}

