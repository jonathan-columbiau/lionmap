#' Helper function that takes in BPCells object and provides Seurat object
#' with this dataset (since Seurat also supports BPCells).
#'
#' @param directory_name Name of directory holding bpcells obj
#' @param metadata Dataframe of metadata for BPCells cells
#'
#' @return Seurat object with BPCells obj as GE matrix
#' @keywords internal
#'
#' @examples
#' LoadBPCellsObjToSeurat(bpcells_dir, metadata_df)
LoadBPCellsObjToSeurat <- function(directory_name, metadata) {
  mat.bpcells = open_matrix_dir(dir = directory_name)
  new_seurat_obj = CreateSeuratObject(counts = mat.bpcells, meta.data = metadata)
  new_seurat_obj
}
