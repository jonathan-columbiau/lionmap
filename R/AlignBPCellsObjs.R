#' Align BPCells Objects so they both have the same genes (useful for getting reference and
#' query datasets in the proper format before finding marker genes).
#'
#' @param bpcells_obj1 First BPCells Object
#' @param bpcells_obj2 Second BPCells Object
#'
#' @return The first BPCells object put in the bpcells_obj1 parameter.
#' @export
#'
#' @examples
#' data("train_ex_data_bpcells")
#' data("test_ex_data_bpcells")
#' train_ex_data_bpcells = AlignBPCellsObjs(train_ex_data_bpcells, test_ex_data_bpcells)
#' test_ex_data_bpcells = AlignBPCellsObjs(test_ex_data_bpcells, train_ex_data_bpcells)
AlignBPCellsObjs <- function(bpcells_obj1, bpcells_obj2) {
  #only keep genes with at least 1 cell expressed and found in both datasets
  genes_expressed_obj1 <- BPCells::rowSums(bpcells_obj1) %>% .[.!=0] %>% names()
  genes_expressed_obj2 <- BPCells::rowSums(bpcells_obj2) %>% .[.!=0] %>% names()
  keep_genes = intersect(genes_expressed_obj1, genes_expressed_obj2)
  bpcells_obj1 = bpcells_obj1[keep_genes,]
  bpcells_obj1
}
