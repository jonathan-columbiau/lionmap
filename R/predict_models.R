predict_models <- function(model, model_name, nonsparse_mat) {
  if (model_name %in% c("ridge", "lasso", "elastic_net")) {
    predict(model, nonsparse_mat, s = "lambda.1se", type = "class") %>% as.character() %>% set_names(rownames(nonsparse_mat)) %>%  return()
  } else if (model_name %in% c("linear_da", "quadratic_da")){
    predict(model,nonsparse_mat %>% as.data.frame())$class %>% as.character() %>% set_names(rownames(nonsparse_mat)) %>% return()
  } else if (model_name %in% c("knn")){
    model %>% predict(nonsparse_mat %>% as.data.frame(), type = "class") %>% as.character() %>% set_names(rownames(nonsparse_mat)) %>% return()
  }  else if (model_name %in% c("rf")){
    model %>% predict(nonsparse_mat, type = "response") %>% .$predictions %>%  as.character() %>% set_names(rownames(nonsparse_mat)) %>% return()
  } else {
    predict(model, nonsparse_mat) %>% as.character() %>% set_names(rownames(nonsparse_mat)) %>% return()
  }
}
