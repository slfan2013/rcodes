read_data = function (input = "C:\\Users\\FanslyFTW\\Documents\\GitHub\\app\\Datasets\\mx 426933_Morrissey_HILIC posESI_negESI_human serum_11-2018 submit_injorder.xlsx")
{
  pacman::p_load(openxlsx, data.table)
  
  if (grepl("xlsx", input)) {
    d <- openxlsx::read.xlsx(input, sheet = 1, colNames = FALSE)
  }else if (grepl("csv", input)) {
    d <- data.table::fread(input)
  }
  data = d
  start_col_index = min(which(diff(which(is.na(data[1,])))>1)+1)
  if(start_col_index==Inf){
    start_col_index = max(which(is.na(data[1,])))+1
  }
  
  start_row_index = min(which(diff(which(is.na(data[,1])))>1)+1)
  if(start_row_index==Inf){
    start_row_index = max(which(is.na(data[,1])))+1
  }
  
  if(start_col_index == -Inf){
    start_col_index = max(which(data[1,]==""))+1
  }
  if(start_row_index == -Inf){
    start_row_index = max(which(data[,1]==""))+1
  }
  
  
  p = data.table(t(data[1:start_row_index,start_col_index:ncol(data)]))
  colnames(p) = as.character(p[1,])
  p = p[-1,]
  if(is.null(dim(p))){
    p = data.table(label = p)
  }
  p$label = make.unique(p$label)
  
  f = data.table(data[start_row_index:nrow(data),1:start_col_index])
  colnames(f) = as.character(f[1,])
  f = f[-1,]
  if(is.null(dim(f))){
    f= data.table(label = f)
  }
  f$label = make.unique(f$label)
  
  e = data.matrix(data[(start_row_index+1):nrow(data),(start_col_index+1):ncol(data)])
  rownames(e) = f$label
  colnames(e) = p$label
  
  result =list(p = p, f = f, e = e)
  
  
  data_matrix = as.matrix(e)
  sum_NA = sum(is.na(data_matrix))
  if (sum_NA > 0) {
    message_from_R = paste0("<em>Please note, your data has ",
                            sum_NA, " missing values. These values will be replace by half-minimum for each compound before normalization.</em>")
  }else {
    message_from_R = ""
  }
  num_NAs = c()
  for(i in 1:nrow(data_matrix)){
    num_NAs[i] = sum(is.na(data_matrix[i,]))
  }
  too_many_na = num_NAs > (0.05 * ncol(data_matrix))
  # data_matrix = data_matrix[!too_many_na,]
  # f = f[!too_many_na,]
  for (i in 1:nrow(data_matrix)) {
    data_matrix[i, is.na(data_matrix[i, ])] = 1/2 * runif(1,min(data_matrix[i,!is.na(data_matrix[i, ])])*0.9,min(data_matrix[i,!is.na(data_matrix[i, ])])*1.1)
  }
  constant_compound = apply(data_matrix,1,sd,na.rm = TRUE) == 0
  data_matrix = data_matrix[!constant_compound,]
  f = f[!constant_compound,]
  if (!"label" %in% names(p)) {
    stop("cannot find 'label' in your dataset. Please download the example file and read the second sheet: 'Explanation'. ")
  }
  if (!"batch" %in% names(p)) {
    stop("cannot find 'batch' in your dataset. Please download the example file and read the second sheet: 'Explanation'.")
  }
  if (!"sampleType" %in% names(p)) {
    stop("cannot find 'sampleType' in your dataset. Please download the example file and read the second sheet: 'Explanation'.")
  }
  if (!"time" %in% names(p)) {
    stop("cannot find 'time' in your dataset. Please download the example file and read the second sheet: 'Explanation'.")
  }
  if (!"qc" %in% unique(p$sampleType)) {
    stop("cannot find 'qc' in your dataset. Please download the example file and read the second sheet: 'Explanation'.")
  }
  if (!"sample" %in% unique(p$sampleType)) {
    stop("cannot find 'sample' in your dataset. Please download the example file and read the second sheet: 'Explanation'.")
  }
  p$time = as.numeric(p$time)
  if (sum(is.na(p$time))>0) {
    stop(paste0("'time' has NA "))
  }
  if (any(table(p$batch[p$sampleType == "qc"]) < 5)) {
    stop("Some batches have too little QCs. At least 5 qc needed for EACH batch.")
  }
  
  
  
  good_index = which(p$sampleType %in% c("qc", "validate", "sample"))
  bad_index = which(!p$sampleType %in% c("qc", "validate", "sample"))
  
  p_good = p[good_index, ]
  data_matrix_good = data_matrix[,good_index]
  
  
  p_bad = p[bad_index, ]
  data_matrix_bad = data_matrix[,bad_index]
  
  sample_rank = rank(p_good$time, ties.method = "first")
  sample_order = order(p_good$time, decreasing = FALSE)
  
  sample_time = p_good$time
  data_matrix = data_matrix_good[, sample_order]
  # p = p[sample_order, ]
  # p$time = order(p$time)
  
  print(length(sample_rank))
  
  # p$sampleType[!p$sampleType %in% c("qc", "validate")] = "sample"
  return(list(e = data_matrix_good, f = f, p = p_good, data_matrix = data_matrix,
              sample_size_summary = as.data.frame.matrix(table(p$batch,
                                                               p$sampleType)), sample_order = sample_order, sample_time = sample_time,
              bad_p = p_bad,bad_data_matrix = data_matrix_bad,
              sample_rank = sample_rank, message_from_R = message_from_R, good_index = good_index, bad_index = bad_index))
}
