aggregate_p_f_e = function(p,f,e){
  
  p_label = p$label
  p$label = NULL
  p$label = p_label
  
  p_plus_e = rbind(cbind(colnames(p),t(p)), cbind(f$label,e))
  
  # add some NA to the top of the f.
  f_label_index = which(p_plus_e[,1]=='label')
  NA_plus_f = rbind(matrix(NA, nrow = f_label_index-1, ncol = ncol(f)),rbind(matrix(colnames(f), nrow = 1),sapply(f, as.character)))
  NA_plus_f = NA_plus_f[,-which(colnames(NA_plus_f) %in% 'label')]
  if(is.null(dim(NA_plus_f))){
    data = cbind(data.table(NA_plus_f), p_plus_e)
  }else{
    if(identical(colnames(NA_plus_f),as.character(NA_plus_f[1,]))){
      NA_plus_f = NA_plus_f[-1,]
    }
    data = cbind(NA_plus_f, p_plus_e)
  }
  
  
  
  
  
  return(data)
}