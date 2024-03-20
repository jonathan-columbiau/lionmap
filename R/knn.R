#' knn  helper function
#'
#' @param reference_dataset PC-transformed dataset in matrix form
#' @param celltype_labels Celltype labels in vector
#'
#' @return knn model trained on evenly split dataset (upsamples if classes aren't evenly split)
#' @importFrom caret upSample
#' @importFrom caret knn3
#' @keywords internal
#'
#' @examples
#' data("iris")
#' dataset = iris[iris$Species %in% c("setosa","versicolor"),]
#' labels = as.factor(dataset$Species)
#' labels = droplevels(labels)
#' dataset = dataset[,1:4]
#' ex_model = knn(dataset, labels)
knn <- function(reference_dataset, celltype_labels) {
  #upsample minority class to make class frequencies equal
  reference_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname = "celltype_labels")
  pairwise_model <- knn3(celltype_labels ~ ., data = reference_dataset, k = 5)
  pairwise_model
}
