polynomial_svm <- function(reference_dataset, celltype_labels) {
  reference_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname =  "celltype_labels")
  pairwise_model <- svm(celltype_labels ~ ., data = reference_dataset, kernel = "polynomial", scale = F)
  pairwise_model
}
