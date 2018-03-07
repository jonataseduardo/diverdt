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

    t_fc[, BS := - log(1 - FST)]
    t_fo[, BS := - log(1 - FST)]
    t_co[, BS := - log(1 - FST)]

    setkeyv(t_fc, c("CHR", "CM", "POS", "SNP"))
    setkeyv(t_fo, c("CHR", "CM", "POS", "SNP"))
    setkeyv(t_co, c("CHR", "CM", "POS", "SNP"))

    aux <- t_fc[t_fo][t_co][, PBS := 0.5 * (BS + i.BS - i.BS.1)]
    aux <- aux[!is.na(PBS)]
    #aux[PBS <0, PBS := 0]

    return(aux[, .(CHR, CM, POS, SNP, PBS)])
    }
