#' Binomial elastic net helper function
#'
#' @param reference_dataset PC-transformed dataset in matrix form
#' @param celltype_labels Celltype labels in vector
#'
#' @return Binomial Elastic Net model trained on evenly split dataset (upsamples if classes aren't evenly split)
#' @importFrom caret upSample
#' @importFrom glmnet cv.glmnet
#' @keywords internal
#'
#' @examples
#' data("iris")
#' dataset = iris[iris$Species %in% c("setosa","versicolor"),]
#' labels = as.factor(dataset$Species)
#' labels = droplevels(labels)
#' dataset = dataset[,1:4]
#' ex_model = binomial_elastic_net(dataset, labels)
binomial_elastic_net <- function(reference_dataset, celltype_labels) {
  upsampled_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname =  "celltype_labels", list = T)
  celltype_labels <- upsampled_dataset[["y"]]
  reference_dataset <- upsampled_dataset[["x"]] %>% as.matrix()
  pairwise_model <- cv.glmnet(x = reference_dataset, y = celltype_labels, family = "binomial", alpha = .5)
  pairwise_model
}
