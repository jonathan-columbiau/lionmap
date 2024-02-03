#' naive_bayes  helper function
#'
#' @param reference_dataset PC-transformed dataset in matrix form
#' @param celltype_labels Celltype labels in vector
#'
#' @return naive_bayes model trained on evenly split dataset (upsamples if classes aren't evenly split)
#' @keywords internal
#'
#' @examples
#' ex_model = naive_bayes(reference_dataset, celltype_labels)
naive_bayes <- function(reference_dataset, celltype_labels) {
  #upsample minority class to make class frequencies equal
  reference_dataset <- caret::upSample(x = reference_dataset, y = celltype_labels, yname = "celltype_labels")
  pairwise_model <- e1071::naiveBayes(celltype_labels ~ ., data = reference_dataset)
  pairwise_model
}
