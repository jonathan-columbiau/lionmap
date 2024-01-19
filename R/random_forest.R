random_forest <- function(reference_dataset, celltype_labels) {
  #upsample minority class to make class frequencies equal
  reference_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname = "celltype_labels")
  pairwise_model <- ranger(celltype_labels ~ ., data = reference_dataset, num.trees = 500, classification = T, replace = T)
  pairwise_model
}
