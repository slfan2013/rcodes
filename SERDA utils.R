remove_outlier = function(v){
  out = boxplot.stats(v)$out
  return(list(value = v[!v%in%out],index = which(v%in%out)))
}
RSD = function(data, robust = FALSE){
  return(apply(data,1,function(x){
    if(robust){x = remove_outlier(x)[[1]]}
    return(sd(x,na.rm=T)/mean(x,na.rm=T))
  }))
}
robust_sd = function(x){
  x = remove_outlier(x)[[1]]
  sd(x,na.rm = TRUE)
}
rm_batch = function(d, batch){ # d=e
  return(t(apply(t(d),2,function(x){
    x/(rep(by(x,batch,mean), table(batch))/mean(x))
  })))
  
}
scale_data = function(d){# d = e. 
  sds = apply(d, 1, sd, na.rm = TRUE)
  # !!! how to deal with zero sds.
  sds[sds == 0] = 1
  
  means = apply(d, 1, mean, na.rm = TRUE)
  data_scale = d
  for(i in 1:nrow(d)){
    data_scale[i,] = (d[i,] - means[i])/sds[i]
  }
  return(list(
    data_scale = data_scale,
    means = means,
    sds = sds
  ))
}
transform = function(e, forward = TRUE, y0 = NULL, lambda = NULL, regular_log = FALSE){
  
  
  if(is.null(y0)){y0=0}
  
  if('numeric' %in% class(e)){
    
    if(forward){
      
      
    }else{
      if(is.null(lambda)){
        stop("You forgot to input lambda.")
      }
      e_after = 0.5 * (2*y0 - lambda * exp(-e) + exp(e))
      
    }
    
  }else{
    
    if(forward){
      if(regular_log){
        e_after = log(e)
      }else{
        e_after = e
        lambda = rowMeans(e)^2
        for(i in 1:nrow(e)){
          e_after[i,] = log(e[i,] - y0 + sqrt((e[i,]-y0)^2 + lambda[i]))
        }
      }
    }else{
      if(regular_log){
        e_after = exp(e)
      }else{
        if(is.null(lambda)){
          stop("You forgot to input lambda.")
        }
        e_after = e
        for(i in 1:nrow(e)){
          e_after[i,] = 0.5 * (2*y0 - lambda[i] * exp(-e[i,]) + exp(e[i,]))
        }
      }
    }
    
  }
  
  
  return(list(e_after, lambda = lambda))
}

ae_model = function(x_train_input, x_test_input,x_train_output, x_test_output, patience = 50, layer_units = c(1400), batch_size = 128, epochs = 2000, verbose = 0,activations = c("elu"), drop_out_rate = 0.05, optimizer = "adam",s_index = NULL, t_index = NULL
                    # , customize_performance_evaluation_set
){
  # seed = seed
  # reticulate::py_config()
  # reticulate::py_set_seed(seed)
  # set.seed(seed)
  
  
  # middles = list()
  # tensorflow::tf$random$set_seed(seed)
  # middles[[length(layer_units)+1]] = keras_model_sequential() # here is the final layer.
  # middles[[length(layer_units)+1]] %>% layer_dense(input_shape = layer_units[length(layer_units)], units = nrow(f),use_bias = TRUE)
  # 
  # for(layer_index in length(layer_units):1){
  #   tensorflow::tf$random$set_seed(seed)
  #   middles[[layer_index]] = keras_model_sequential(name = paste0('middle',layer_index))
  #   if(layer_index==1){
  #     middles[[layer_index]] %>% layer_dense(input_shape = nrow(f), units = layer_units[layer_index], activation = activations[layer_index],use_bias = TRUE) # this is the first mid-layer
  #   }else{
  #     middles[[layer_index]] %>% layer_dense(input_shape = layer_units[layer_index-1], units = layer_units[layer_index], activation = activations[layer_index],use_bias = TRUE) # this is other mid-layer
  #   }
  # }
  # model = middles[[1]]%>%layer_dropout(drop_out_rate)
  # for(i in 2:length(middles)){
  #   model_temp = middles[[i]]
  #   model %>% model_temp
  # }
  
  early_stopping <- callback_early_stopping(patience = patience)
  tensorflow::tf$random$set_seed(seed)
  model <- keras_model_sequential() 
  # model <- keras_model_sequential() 
  model %>% 
    layer_dense(units = layer_units, activation =  activations, input_shape = nrow(f)) %>% 
    layer_dropout(rate = drop_out_rate) %>% 
    layer_dense(units = nrow(f)
                # , activation = 'softmax'
    )
  
  model %>% compile(
    loss = "mean_absolute_error",#1997171
    # loss = "mean_squared_error",#0.2111578
    # loss = "mean_absolute_percentage_error",#0.3923684
    # loss = "mean_squared_logarithmic_error",#0.2652711
    # loss = "squared_hinge",#1.961605
    # loss = "hinge",#3.870686
    # loss = "logcosh",#0.2002801
    # loss = "huber_loss",#0.2005046
    # loss = "categorical_crossentropy",#2.752495
    # loss = "sparse_categorical_crossentropy",#error
    # loss = "binary_crossentropy",#1.201967
    # loss = "kullback_leibler_divergence",#2.403364
    # loss = "poisson",#error
    # loss = "cosine_proximity",#error
    optimizer = optimizer
  )
  model %>% save_model_hdf5("model.h5")
  
  # Winsorizing
  # for(i in 1:ncol(x_train)){
  #   
  #   quan = quantile(x_train[,i], probs = c(0.01,0.99))
  #   
  # }
  
  
  model %>% fit(
    x = x_train_input,
    y = x_train_output,
    epochs = epochs,
    verbose = verbose,
    batch_size = batch_size,
    validation_data = list(x_test_input, x_test_output), 
    callbacks = list(early_stopping),
    view_metrics = FALSE
  )
  
  epoch_n = length(model$history$epoch)#-patience
  # epoch_n  
  
  
  
  x_all_input = rbind(x_train_input,x_test_input)
  x_all_output = rbind(x_train_output,x_test_output)
  
  final_model <- load_model_hdf5("model.h5")
  final_model %>% fit(
    x = x_all_input,
    y = x_all_output,
    epochs = epoch_n,
    verbose = verbose,
    batch_size = batch_size,
    view_metrics = FALSE
  ) 
  
  
  
  
  
  
  if(is.null(s_index)){
    return(list(
      final_model = final_model,
      model = model
    ))
  }else{
    x_all = x_all_input
    # pred_train = x_all_output
    # 
    # model2 <- load_model_hdf5("model.h5")
    # model2 %>% fit(
    #   x = x_all_input,
    #   y = pred_train,
    #   epochs = epoch_n,
    #   verbose = verbose,
    #   batch_size = batch_size
    # ) 
    
    
    pred = predict(final_model,  x_all)
    
    x_all[s_index,] = pred[1:nrow(x_current_train_output),]
    x_all[t_index,] = pred[(nrow(x_current_train_output)+1):nrow(pred),]
    
    
    
    return(list(
      final_model = final_model,
      pred_train = x_all,
      model = model
    ))
    
  }
  
  
}