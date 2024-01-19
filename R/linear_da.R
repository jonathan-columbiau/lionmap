linear_da <- function(reference_dataset, celltype_labels) {
  #upsample minority class to make class frequencies equal
  reference_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname = "celltype_labels")
  pairwise_model <- MASS::lda(celltype_labels ~ ., data = reference_dataset)
  pairwise_model
}
