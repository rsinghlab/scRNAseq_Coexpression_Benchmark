loadPred <- function(filename, model_name){
  obj <- readRDS(filename);
  pred_list <- obj[[model_name]]$est_network_list;
  return (pred_list)
}