#' linear_da  helper function
#'
#' @param reference_dataset PC-transformed dataset in matrix form
#' @param celltype_labels Celltype labels in vector
#'
#' @return linear_da model trained on evenly split dataset (upsamples if classes aren't evenly split)
#' @export
#'
#' @examples
#' ex_model = linear_da(reference_dataset, celltype_labels)
linear_da <- function(reference_dataset, celltype_labels) {
  #upsample minority class to make class frequencies equal
  reference_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname = "celltype_labels")
  pairwise_model <- MASS::lda(celltype_labels ~ ., data = reference_dataset)
  pairwise_model
}
