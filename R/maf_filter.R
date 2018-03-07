#' Filter from data set SNPs with MAF bellow a threshold
#'
#' \code{maf_filter()} Filter SNPs with MAF bellow a threshold
#' 
#' @import data.table
#' @export 
#'
#' @param freq_data data.table. Data with one allele per line. 
#' @param maf_th float. MAF threshold. Default maf_th = NULL to consider 
#' every polymorphic site
#'
#' @return data.table  

maf_filter <-
  function(freq_data,
           maf_th = NULL,
           inplace = FALSE){

  if(inplace){
    f_data <- copy(freq_data)
  }else{
    f_data <- freq_data 
  }
   
  f_data[, N0 := .N, keyby = .(CHR, POS)]
  if(is.null(maf_th)){
    f_data[, maf_th := 1 / NCHROBS]
    f_data[AF > maf_th, N1 := .N, by = .(CHR, POS)]
  }else{
    f_data[, maf_th := maf_th]
    f_data[AF > maf_th, N1 := .N, by = .(CHR, POS)]
  }
  f_data <- f_data[N0 == N1]
  f_data[, c('N0', 'N1', 'maf_th') := NULL]

  setkey(f_data, CHR, POS)
  return(f_data[])
  }
