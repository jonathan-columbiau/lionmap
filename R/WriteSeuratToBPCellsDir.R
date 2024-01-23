#' Helper function to write counts stored in Seurat object in the field
#' assays/RNA/counts to bpcells directory
#'
#' @param seurat_obj Seurat object
#' @param dir_name BPCells directory name to write to
#' @param cell_label_column Cell label column name in Seurat object metadata
#'
#' @return BPCells GE matrix
#' @export
#'
#' @examples WriteSeuratToBPCellsDir(seurat_obj, bpcells_dir, "cellid")
WriteSeuratToBPCellsDir <- function(seurat_obj,dir_name,cell_label_column = "cellid") {
  write_matrix_dir(
    mat = seurat_obj@assays$RNA@counts,
    dir = dir_name
  )
}
