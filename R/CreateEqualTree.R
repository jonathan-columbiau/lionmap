#' This function is meant for users who don't want to use the tree-based/hierarchical functionality of the package.
#' It creates a tree that assumes that all the celltype labels are at the same level of the tree.
#'
#'
#' @param cell_labels Vector of celltype labels we want to create a tree with. Don't include duplicate celltype labels.
#' @param rootnode_name Optional parameter giving the name of the rootnode of the tree.
#' Recommended to keep as "Unmapped".
#'
#' @return Tree in treedata format.
#' @export
#'
#' @examples
#' equal_tree = lionmap::CreateEqualTree(cell_labels = paste0(1:20,"n"),rootnode_name = "Root")
CreateEqualTree <- function(cell_labels, rootnode_name = "Rootnode") {
  CreateHierarchy(df_hierarchy = data.frame(parent = rootnode_name,
                                            node = cell_labels))
}
