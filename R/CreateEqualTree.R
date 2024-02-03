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
#' CreateEqualTree(celltype_labels)
CreateEqualTree <- function(cell_labels, rootnode_name = "Rootnode") {
  tree_newick_format <- paste0("(",stringr::str_c(unique(na.omit(cell_labels)), collapse = ","),")",rootnode_name,";") %>% .[!is.na(.)]
  treeio::read.newick(textConnection(tree_newick_format), node.label = "label") %>% treeio::as.treedata()
}
