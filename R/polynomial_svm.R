#' polynomial_svm  helper function
#'
#' @param reference_dataset PC-transformed dataset in matrix form
#' @param celltype_labels Celltype labels in vector
#'
#' @return polynomial_svm model trained on evenly split dataset (upsamples if classes aren't evenly split)
#' @keywords internal
#'
#' @examples
#' ex_model = polynomial_svm(reference_dataset, celltype_labels)
polynomial_svm <- function(reference_dataset, celltype_labels) {
  reference_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname =  "celltype_labels")
  pairwise_model <- svm(celltype_labels ~ ., data = reference_dataset, kernel = "polynomial", scale = F)
  pairwise_model
}
