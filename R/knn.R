knn <- function(reference_dataset, celltype_labels) {
  #upsample minority class to make class frequencies equal
  reference_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname = "celltype_labels")
  pairwise_model <- caret::knn3(celltype_labels ~ ., data = reference_dataset, k = 5)
  pairwise_model
}
