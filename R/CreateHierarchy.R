
#' Create Hierarchical Tree Structure from CSV with parent and child nodes
#' named in a csv file
#'
#' @param csv_file_path A CSV file path with a column named 'parent' and a column
#' named 'node' giving the relationships between each node and its parent node.
#' This info will be used to create the tree.
#'
#' @return a treedata object
#' @export
#'
#' @examples
CreateHierarchy <- function(csv_file_path) {
  # Read the CSV file
  tree_tbl = read.csv(csv_file_path)

  # Replace spaces with underscores in 'node' and 'parent' columns
  tree_tbl$node = gsub(" ", "_", tree_tbl$node)
  tree_tbl$parent = gsub(" ", "_", tree_tbl$parent)

  # Add dummy 'branch_length' and 'trait' columns
  tree_tbl$branch_length = NA
  tree_tbl$trait = NA

  # Convert the data frame to a tibble
  tree_tbl <- tree_tbl %>% tibble::as_tibble()

  # Convert the tibble to a treedata object
  tree_tidy = as.treedata(tree_tbl)

  return(tree_tidy)
}
