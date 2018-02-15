#' Estimate the population branch size of for every SNP 
#'
#'\code{pbs_fst()} estimates the branch size for every SNP 
#' 
#' @import data.table
#' @export
#' 
#' @param  focal_pop: data.table with the following columns: CHR, SNP, CM, POS,
#' NCHROBS, POP, VAR, AF.
#' @param close_pop: data.table
#' @param out_pop: data.table
#' 
#' @return data.table with the estimated branch size for every SNP


pbs_fst <-
    function(focal_pop, close_pop, out_pop){
  
    t_fc <- wc_fst(list(focal_pop, close_pop))
    t_fo <- wc_fst(list(focal_pop, out_pop))
    t_co <- wc_fst(list(close_pop, out_pop))

    t_fc[, T := - log(1 - FST)]


    }
