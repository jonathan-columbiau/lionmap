binomial_ridge <- function(reference_dataset, celltype_labels) {
  upsampled_dataset <- upSample(x = reference_dataset, y = celltype_labels, yname =  "celltype_labels", list = T)
  celltype_labels <- upsampled_dataset[["y"]]
  reference_dataset <- upsampled_dataset[["x"]] %>% as.matrix()
  pairwise_model <- cv.glmnet(x = reference_dataset, y = celltype_labels, family = "binomial", alpha = 0)
  pairwise_model
}
