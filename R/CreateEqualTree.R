CreateEqualTree <- function(cell_labels, rootnode_name = "Unmapped") {
  tree_newick_format <- paste0("(",str_c(unique(na.omit(cell_labels)), collapse = ","),")",rootnode_name,";") %>% .[!is.na(.)]
  read.newick(textConnection(tree_newick_format), node.label = "label") %>% as.treedata()
}
