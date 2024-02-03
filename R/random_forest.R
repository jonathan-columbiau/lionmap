#' Random Forest  helper function (uses ranger backend)
#'
#' @param reference_dataset PC-transformed dataset in matrix form
#' @param celltype_labels Celltype labels in vector
#'
#' @return rf model trained on evenly split dataset (upsamples if classes aren't evenly split)
#' @keywords internal
#'
#' @examples
#' ex_model = random_forest(reference_dataset, celltype_labels)
random_forest <- function(reference_dataset, celltype_labels) {
  #upsample minority class to make class frequencies equal
  reference_dataset <- caret::upSample(x = reference_dataset, y = celltype_labels, yname = "celltype_labels")
  pairwise_model <- ranger::ranger(celltype_labels ~ ., data = reference_dataset, num.trees = 500, classification = T, replace = T)
  pairwise_model
}
