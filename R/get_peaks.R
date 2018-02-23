#' Return the windows and the SNPs with the higest scores
#'
#'\code{get_peaks()} Return the windows and the SNPs with the higest scores
#' 
#' @import data.table
#' @export
#' 
#' @param data_raw: data.table with every SNP position 
#' @param data_sum: data.table with the following columns: CHR, SNP, CM, POS, col_name
#' @param col_name: Name of the score column 
#' @param score_th: double with the score thereshold 
#' @param greater: BOOLEAN deafault TRUE 
#' 
#' @return data.table 

get_peaks <-
  function(data_raw, data_sum, col_name, score_th, window_s, greater = TRUE){

    score_th <- 0.01
    col_name <- 'p.value'

    setnames(data_sum, col_name, 'score_col')

    setkey(data_raw, CHR, POS, CM, SNP)
    data_raw[, idx := .I]

    setkey(data_sum, CHR, POS, CM, SNP)

    DT <- data_sum[score_col < score_th][data_raw]
    list_idx <- DT[!is.na(score_col)][, idx]

    n <- length(list_idx)

    i <- 1
    
    id_aux <- list()
    
    for(i in 1:n){
      id_aux[[i]] <- list_idx[(1 + i):n] - list_idx[1:(n - i)] <= window_s
      if (!any(id_aux[[i]])) break
    }
    id_aux
    i

    list_idx[2:(n - 1)][id_aux]

  }

