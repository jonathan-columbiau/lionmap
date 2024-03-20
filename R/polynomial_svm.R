#' polynomial_svm  helper function
#'
#' @param reference_dataset PC-transformed dataset in matrix form
#' @param celltype_labels Celltype labels in vector
#'
#' @return polynomial_svm model trained on evenly split dataset (upsamples if classes aren't evenly split)
#' @importFrom caret upSample
#' @importFrom e1071 svm
#' @keywords internal
#'
#' @examples
#' data("iris")
#' dataset = iris[iris$Species %in% c("setosa","versicolor"),]
#' labels = as.factor(dataset$Species)
#' labels = droplevels(labels)
#' dataset = dataset[,1:4]
#' ex_model = polynomial_svm(dataset, labels)
polynomial_svm <- function(reference_dataset, celltype_labels) {
  reference_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname =  "celltype_labels")
  pairwise_model <- svm(celltype_labels ~ ., data = reference_dataset, kernel = "polynomial", scale = F)
  pairwise_model
}
