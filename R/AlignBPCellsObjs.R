AlignBPCellsObjs <- function(bpcells_obj1, bpcells_obj2) {
  #only keep genes with at least 1 cell expressed and found in both datasets
  genes_expressed_obj1 <- rowSums(bpcells_obj1) %>% .[.!=0] %>% names()
  genes_expressed_obj2 <- rowSums(bpcells_obj2) %>% .[.!=0] %>% names()
  keep_genes = intersect(genes_expressed_obj1, genes_expressed_obj2)
  bpcells_obj1 = bpcells_obj1[keep_genes,]
  bpcells_obj1
}
