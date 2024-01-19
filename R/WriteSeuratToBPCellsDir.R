WriteSeuratToBPCellsDir <- function(seurat_obj,dir_name,cell_label_column = "cellid") {
  write_matrix_dir(
    mat = seurat_obj@assays$RNA@counts,
    dir = dir_name
  )
}
