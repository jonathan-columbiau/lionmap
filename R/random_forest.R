#' Random Forest  helper function (uses ranger backend)
#'
#' @param reference_dataset PC-transformed dataset in matrix form
#' @param celltype_labels Celltype labels in vector
#'
#' @return rf model trained on evenly split dataset (upsamples if classes aren't evenly split)
#' @importFrom caret upSample
#' @importFrom ranger ranger
#' @keywords internal
#'
#' @examples
#' data("iris")
#' dataset = iris[iris$Species %in% c("setosa","versicolor"),]
#' labels = as.factor(dataset$Species)
#' labels = droplevels(labels)
#' dataset = dataset[,1:4]
#' ex_model = random_forest(dataset, labels)
random_forest <- function(reference_dataset, celltype_labels) {
  #upsample minority class to make class frequencies equal
  reference_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname = "celltype_labels")
  pairwise_model <- ranger(celltype_labels ~ ., data = reference_dataset, num.trees = 500, classification = T, replace = T)
  pairwise_model
}
