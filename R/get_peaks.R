#' Return the windows and the SNPs with the higest sginificance scores
#'
#'\code{get_peaks()} Return the windows and the SNPs with the higest scores
#' 
#' @import data.table
#' @export
#' 
#' @param data_raw: data.table with every SNP position 
#' @param data_sum: data.table with the following columns: CHR, SNP, CM, POS, col_name
#' @param stat_col: Name of the stat column 
#' @param score_col: Name of the score column 
#' @param score_th: double with the score thereshold 
#' @param window_s: integer with the size of the window used in the rolling
#' statistics 
#' @param greater: BOOLEAN deafault TRUE 
#' 
#' @return data.table 

get_peaks <-
  function(data_raw, data_sum, stat_col, score_col, score_th, window_s, greater = FALSE){

    data_r <- copy(data_raw)
    data_s <- copy(data_sum)

    setnames(data_r, stat_col, 'stat_col')
    setnames(data_s, score_col, 'score_col')

    setkey(data_r, CHR, POS, CM, SNP)
    setkey(data_s, CHR, POS, CM, SNP)
    data_r[, idx := .I]

    if(greater){
      DT <- data_r[data_s[score_col >= score_th]]
    }else{
      DT <- data_r[data_s[score_col <= score_th]]
    }

    conc_i <- 
      unlist(lapply(DT[, idx], function(i) i:(i + window_s - 1)))
    
    setkey(data_r, idx)
    data_top <- data_r[conc_i]
    data_top[, first_idx := DT[rep(1:.N, each=window_s), idx]]

    data_out <- data_top[data_top[, .I[stat_col == max(stat_col)], by = first_idx]$V1]
    
    DW <- DT[,.(idx, W_ID)]
    data_out <- data_out[DW, on = c(first_idx = 'idx')]
    data_out[, c('idx', 'first_idx') := NULL]

    setnames(data_out, 'stat_col', stat_col)
    return(data_out)
  }

