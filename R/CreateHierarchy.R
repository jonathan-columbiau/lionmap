
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
#' @examples tree_df = data.frame(
#' parent = c("p1","p1","p1","p1","p1",
#'            "p2","p2","p2","p2","p2","p3","p3","p3","p3","p3",
#'            "p4","p4","p4","p4","p4","gp1","gp1","gp2",
#'            "gp2","rootnode","rootnode","rootnode"),
#' node = c("n1","n2","n3","n4","n5",
#'          "n6","n7","n8","n9","n10","n11","n12","n13","n14",
#'          "n15","n16","n17","n18","n19","n20","p1","p2",
#'          "p3","p4","gp1","gp2","rootnode")
#' )
#'
#' new_hierarchy = CreateHierarchy(csv_file_path = NULL, df_hierarchy = tree_df)
CreateHierarchy <- function(csv_file_path = NULL, df_hierarchy = NULL) {
  if(!is.null(csv_file_path)) {
    # Read the CSV file
    tree_tbl = read.csv(csv_file_path)

    # Replace spaces with underscores in 'node' and 'parent' columns
    tree_tbl$node = gsub(" ", "_", tree_tbl$node)
    tree_tbl$parent = gsub(" ", "_", tree_tbl$parent)

    # Add dummy 'branch_length' and 'trait' columns
    tree_tbl$branch_length = NA
    tree_tbl$trait = NA

    # Convert the data frame to a tibble
    tree_tbl <- tree_tbl %>% tibble::as_tibble() %>% dplyr::select(node, parent, branch_length, trait)

    # Convert the tibble to a treedata object
    tree_tidy = as.treedata(tree_tbl)
    return(tree_tidy)
  }
  else {
    tree_tbl <- df_hierarchy %>% as_tibble()
    tree_tbl$branch_length = NA
    tree_tbl$trait = NA
    tree_tbl <- tree_tbl %>% dplyr::select(node, parent, branch_length, trait)
    tree_tidy = as.treedata(tree_tbl)
    return(tree_tidy)
  }

}
