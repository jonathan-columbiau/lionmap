#' Helper function that takes in a vector of names, and provides a list with each pairwise combination
#' of names.
#'
#' @param cluster.names Vector of names/celltype labels
#'
#' @return List of pairwise combinations
#' @export
#'
#' @examples
#' pairwise_combinations(names_vector)
pairwise_combinations <- function(cluster.names) {
  name_list <- list()
  cluster.names <- unique(cluster.names)
  list_index_to_add = 1
  for (i in 1:(length(cluster.names))) {
    for(j in (i + 1): length(cluster.names)) {
      if(i != length(cluster.names)) {
        name_list[[list_index_to_add]] <- list(cluster1 = as.character(cluster.names[i]), cluster2 = as.character(cluster.names[j]))
        list_index_to_add <- list_index_to_add + 1
      }
    }
  }
  name_list
}
