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
  function(data_raw, data_sum, col_name, score_th, step_s, greater = TRUE){


    data_raw <- pbs_data
    data_sum <- pbs_mean

    
    score_th <- 0.01
    col_name <- 'p.value'
    window_s = 20

    setnames(data_sum, col_name, 'score_col')

    setkey(data_raw, CHR, POS, CM, SNP)
    data_raw[, idx := .I]

    setkey(data_sum, CHR, POS, CM, SNP)

    DT <- data_raw[data_sum[score_col < score_th]]

    conc_i <- 
      unlist(lapply(DT[, idx], function(i) i:(i + window_s - 1)))

    idi <- DT[rep(1:.N, each=window_s), idx]
    
    setkey(data_raw, idx)
    data_top <- data_raw[conc_i]
    data_top[, first_idx := DT[rep(1:.N, each=window_s), idx]]

    data_top[data_top[, .I[PBS == max(PBS)], by = first_idx]$V1]

    
    data_raw[res$seq, max(PBS), by = res$idx]
  }

