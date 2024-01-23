#' A helper function that takes in a ge reference dataset that has been transformed
#' using PCA, and a set of celltype labels in factor format, and produces a set of
#' models used to differentiate celltypes.
#'
#' @param reference_dataset PC-transformed dataset
#' @param celltype_labels Vector of celltype labels
#' @param models_to_include Optional vector giving model names to include. If include
#' this parameter, include a subset of the following (make sure the names match or it won't work):
#' "linear_svm", "polynomial_svm", "naive_bayes", "ridge", "lasso", "elastic_net",
#'  "linear_da", "knn", "rf", "quadratic_da"
#'
#' @return List of Models
#' @export
#'
#' @examples
#' CreateAllModels(reference_dataset, celltype_labels, models_to_include = NULL)
CreateAllModels <- function(reference_dataset, celltype_labels, models_to_include = NULL) {
  output_model_list <- vector(mode = "list", length = 10)
  names(output_model_list) <- c("linear_svm", "polynomial_svm", "naive_bayes", "ridge", "lasso", "elastic_net", "linear_da", "knn", "rf", "quadratic_da")
  function_list <- c(linear_svm,polynomial_svm,naive_bayes,binomial_ridge,binomial_lasso,binomial_elastic_net,linear_da,knn,random_forest,quadratic_da)
  if(!is.null(models_to_include)) {
    model_indices_to_keep <- which(names(output_model_list) %in% models_to_include)
    output_model_list <- output_model_list[model_indices_to_keep]
    function_list <- function_list[model_indices_to_keep]
  }
  for (i in 1:length(function_list)) {
    model <- function_list[[i]](reference_dataset, celltype_labels)
    output_model_list[[i]] <- model
  }
  #return the list of models
  output_model_list
}
