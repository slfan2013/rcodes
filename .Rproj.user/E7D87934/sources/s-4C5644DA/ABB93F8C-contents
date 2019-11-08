read_data = function(path = "C:\\Users\\Sili\\Desktop\\projects\\mx 399706 Chris Morrissey\\mx 399706 Chris Morrissey, mouse serum, March 2018.xlsx", sheet  = 1){
  library(data.table)
  
  
  # path = "C:\\Users\\Sili\\Documents\\Github\\Bajaj_2_5_2019\\Serum\\Bajaj complete urine and serum data transposed 12.31.18pm 90day 6month.csv"
  
  if(grepl("xlsx",path)){
    data = readxl::read_excel(path, col_names = FALSE, sheet  = sheet )
    data = data.table(data)
  }else{
    data = fread(path)
    data[data=='']=NA
  }
  
  
  
  sample_col_range = min(which(!is.na(data[1,]))):ncol(data)
  sample_row_range = 1:min(which(!is.na(data[[1]])))
  compound_col_range = 1:(min(which(!is.na(data[1,]))))
  compound_row_range = (min(which(!is.na(data[[1]])))):nrow(data)
  
  p = t(data[sample_row_range,sample_col_range,with=F])
  colnames(p) = p[1,]
  p = p[-1,]
  p = p[,c(ncol(p),1:(ncol(p)-1))]
  p = data.table(p)
  
  p = as.data.table(p)
  
  colnames(p) = make.unique(colnames(p), sep = "_")
  if(!"label"%in%colnames(p)){
    stop("Cannot find 'label' in your data. Please check the data format requirement.")
  }
  if(sum(is.na(p$label))>0){
    p$label[is.na(p$label)] = "na"
  }
  
  
  
  
  f = data[compound_row_range,compound_col_range,with=F]
  colnames(f) = as.character(f[1,])
  f = f[-1,]
  f = f[,c(ncol(f),1:(ncol(f)-1)),with=F]
  
  f = as.data.table(f)
  colnames(f) = make.unique(colnames(f), sep = "_")
  if(sum(is.na(f$label))>0){
    f$label[is.na(f$label)] = "na"
  }
  
  
  e = data[compound_row_range, sample_col_range, with = F]
  colnames(e) = as.character(e[1,])
  colnames(e)[is.na(colnames(e))] = "na"
  e = e[-1,]
  
  e = data.table(label = e$label, sapply(e[,-1,with=F], function(x){
    as.numeric(x)
  }))
  
  colnames(e) = make.unique(colnames(e), sep = "_")
  e$label[is.na(e$label)] = "na"
  e$label = f$label
  colnames(e) = c("label",p$label)
  
  
  e_matrix = data.matrix(e[,-1,with=F])
  
  return(list(p = p, f = f, e = e, e_matrix = e_matrix))
}