#' Helper function that loads BPCells obj already created
#'
#' @param directory_name Name of directory with BPCells Object.
#'
#' @return BPCells obj
#' @keywords internal
#'
#' @examples
#' LoadBPCellsObj(bpcells_dir_name)
LoadBPCellsObj <- function(directory_name) {
  open_matrix_dir(dir = directory_name)
}
