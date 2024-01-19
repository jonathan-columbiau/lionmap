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
