LoadBPCellsObjToSeurat <- function(directory_name, metadata) {
  mat.bpcells = open_matrix_dir(dir = directory_name)
  new_seurat_obj = CreateSeuratObject(counts = mat.bpcells, meta.data = metadata)
  new_seurat_obj
}
