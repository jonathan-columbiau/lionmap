Classify <- function(bpcells_query, models, tree_struc, prop_max_threshold = .66) {
  # initial rootnode level - do all at once. All tree_struc have rootnodes, even if all clusters are at one level
  rootnode <- tree_struc %>%
    rootnode() %>%
    nodelab(tree_struc, .)
  tipnodes <- nodelab(tree_struc, offspring(tree_struc, rootnode, type = "tips"))
  # Normalize by reads-per-cell
  bpcells_query <- multiply_cols(bpcells_query, 1 / Matrix::colSums(bpcells_query))
  # Log normalization
  bpcells_query <- log1p(bpcells_query * 10000) # Log normalization
  # get internal nodes in hierarchical order
  internal_nodes <- c(rootnode, offspring(tree_struc, rootnode) %>% nodelab(tree_struc, .) %>% .[-which(. %in% nodelab(tree_struc, offspring(tree_struc, rootnode, type = "tips")))])
  internal_node_assignment <- vector(mode = "list", length = length(internal_nodes)) %>% set_names(internal_nodes)
  query_cells <- bpcells_query
  internal_node_assignment[[rootnode]] <- colnames(query_cells)
  final_classifications <- vector(mode = "character")
  for (j in 1:length(internal_nodes)) { # iterate node
    print(j)
    ## track cells that don't go past internal nodes
    node <- internal_nodes[j]
    res_list <- vector(mode = "list", length = length(models[[node]]))
    cells <- internal_node_assignment[[node]]
    ## final classification for cells that don't go past internal node = internal node
    for (i in 1:length(models[[node]])) { # iterate over models, classify
      first_lev_avg_counts <- models[[node]][[i]]$avg_log_exp # scale
      first_lev_std_counts <- models[[node]][[i]]$stdev
      first_lev_markers <- models[[node]][[i]]$pc_loadings %>% colnames()
      first_lev_bpcells <- query_cells[first_lev_markers, cells] %>% t() # select markers and cells of this internal node level
      first_lev_bpcells <- first_lev_bpcells %>%
        add_cols(-first_lev_avg_counts) %>%
        multiply_cols(1 / first_lev_std_counts)
      first_lev_pc_loadings <- models[[node]][[i]]$pc_loadings[, first_lev_markers] %>% t()
      # transform data using matrix multiplication operator %*%
      first_lev_bpcells <- first_lev_bpcells %*% first_lev_pc_loadings
      nonsparse_mat <- first_lev_bpcells %>% as.matrix()
      # use models
      first_lev_models <- models[[node]][[i]]$Models
      res_list[[i]] <- purrr::map2(first_lev_models, first_lev_models %>% names(), predict_models, nonsparse_mat) %>%
        as.data.frame() %>%
        set_colnames(paste0(colnames(.), "_", i)) %>%
        t()
    }

    count_threshold <- (child(tree_struc,node) %>% length() - 1) * length(models[[1]][[1]][["Models"]]) * prop_max_threshold #2/3 of max score (num matchups for each class - 1)*num models per pairwise comparison

    obs_above_threshold <- res_list %>%
      plyr::rbind.fill.matrix(res_list) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(everything(), names_to = "obs", values_to = "class") %>%
      dplyr::group_by(obs) %>%
      dplyr::count(class, name = "count") %>%
      dplyr::group_by(obs) %>%
      dplyr::filter(count == max(count))
    #filter obs with mult max classes
    tied_obs <- obs_above_threshold %>% group_by(obs) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(obs)
    obs_above_threshold <- obs_above_threshold %>%
      dplyr::filter(!obs %in% tied_obs) %>%
      dplyr::filter(count >= count_threshold) %>%
      pull(class, name = obs)
    ## assign "tip" cells to final classification
    tip_cells <- obs_above_threshold[obs_above_threshold %in% tipnodes]
    final_classifications <- final_classifications %>% append(tip_cells)
    ## assign "stuck" cells to final classification - ties/not threshold
    stuck_cells <- rownames(first_lev_bpcells)[!rownames(first_lev_bpcells) %in% names(final_classifications) & !rownames(first_lev_bpcells) %in% names(obs_above_threshold)] %>% set_names(., .)
    stuck_cells <- rep(node, length(stuck_cells))
    names(stuck_cells) <- rownames(first_lev_bpcells)[!rownames(first_lev_bpcells) %in% names(final_classifications) & !rownames(first_lev_bpcells) %in% names(obs_above_threshold)]
    final_classifications <- final_classifications %>% append(stuck_cells)
    ## update internal node assignment list
    if (obs_above_threshold[!obs_above_threshold %in% tipnodes & !obs_above_threshold %in% names(final_classifications)] %>% length() > 0) {
      obs_above_threshold <- obs_above_threshold[obs_above_threshold %in% internal_nodes]
      obs_above_threshold <- split(obs_above_threshold, obs_above_threshold)
      # unit test 1: all remaining cells assigned to internal nodes
      test_that("all remaining cells assigned to internal nodes", {
        expect_contains(internal_nodes, names(obs_above_threshold))
      })
      for (name in names(obs_above_threshold)) {
        internal_node_assignment[[name]] <- obs_above_threshold[[name]] %>% names()
      }
    }
  }
  test_that("expected number of elements returned", {
    expect_equal(length(final_classifications), ncol(bpcells_query))
  })
  #order match bpcells_query
  final_classifications[match(colnames(bpcells_query),names(final_classifications))]
}
