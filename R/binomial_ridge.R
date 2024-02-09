#' Binomial ridge helper function
#'
#' @param reference_dataset PC-transformed dataset in matrix form
#' @param celltype_labels Celltype labels in vector
#'
#' @return Binomial ridge model trained on evenly split dataset (upsamples if classes aren't evenly split)
#' @keywords internal
#'
#' @examples
#' ex_model = binomial_ridge(reference_dataset, celltype_labels)
binomial_ridge <- function(reference_dataset, celltype_labels) {
  upsampled_dataset <- caret::upSample(x = reference_dataset, y = celltype_labels, yname =  "celltype_labels", list = T)
  celltype_labels <- upsampled_dataset[["y"]]
  reference_dataset <- upsampled_dataset[["x"]] %>% as.matrix()
  pairwise_model <- glmnet::cv.glmnet(x = reference_dataset, y = celltype_labels, family = "binomial", alpha = 0)
  pairwise_model
}