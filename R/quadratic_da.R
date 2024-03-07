#' Quadratic discriminant analysis  helper function
#'
#' @param reference_dataset PC-transformed dataset in matrix form
#' @param celltype_labels Celltype labels in vector
#'
#' @return quadratic_da model trained on evenly split dataset (upsamples if classes aren't evenly split)
#' @keywords internal
#'
#' @examples
#' data("iris")
#' dataset = iris[iris$Species %in% c("setosa","versicolor"),]
#' labels = as.factor(dataset$Species)
#' labels = droplevels(labels)
#' dataset = dataset[,1:4]
#' ex_model = quadratic_da(dataset, labels)
quadratic_da <- function(reference_dataset, celltype_labels) {
  #upsample minority class to make class frequencies equal
  reference_dataset <- caret::upSample(x = reference_dataset, y = celltype_labels, yname = "celltype_labels")
  pairwise_model <- MASS::qda(celltype_labels ~ ., data = reference_dataset)
  pairwise_model
}
