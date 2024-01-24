#' Helper function which uses the fread function from the data.table library to
#' read in a GE dataset in BPCells format.
#'
#' @param filename Name of file
#' @param dir_name Name of BPCells obj
#'
#' @return BPCells-formatted GE matrix
#' @keywords internal
#'
#' @examples
#' LoadBPCellsFromTable(filename, dir_name = "bpcells_dataset")
LoadBPCellsFromTable <- function(filename, dir_name = "new_dataset") {
  dataset <- data.table::fread(filename) %>% .[,2:ncol(.)] %>% as.data.frame()
  removed_gene_names <- dataset$gene_name[duplicated(dataset$gene_name)]
  dataset <- dataset[!duplicated(dataset$gene_name), ]
  rownames(dataset) <- dataset$gene_name
  dataset <- dataset %>% .[,2:ncol(.)]
  #write bpcells
  dataset <- dataset %>% as("Matrix") %>% as("dgCMatrix") %>% as("IterableMatrix") %>% write_matrix_dir(paste0(dir_name), overwrite = T)
  open_matrix_dir(dir_name)
}
