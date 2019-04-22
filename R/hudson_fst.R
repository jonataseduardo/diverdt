#' Estimate the FST using the Hudson formula  
#'
#'\code{wc_fst()} estimates the Fst for every variant 
#' 
#' @import data.table
#' @export
#' 
#' @param  pop_list: list of where each element is a data.table with the
#' following columns: CHR, SNP, CM, POS, NCHROBS, POP, VAR, AF
#' 
#' @return data.table with Fst estimated for every SNP 



hudson_fst <-
    function(pop_list){

    r <- length(pop_list)
    if(r == 2){
      pop_dt <- rbindlist(pop_list)

      setkeyv(pop_dt, c("CHR", "CM", "POS", "VAR", "SNP"))

      pop_dt <- pop_dt[pop_dt[, .N , by = .(CHR, CM, POS, VAR, SNP)][N == r]]
      pop_dt[, N := NULL]
      pop_dt[, G := .GRP, by = .(CHR, CM, POS, VAR, SNP)]

      # Assuming biallelic data select one variant to estimate the Fst 
      pop_dt <- pop_dt[G %% 2 != 0]
      pop_dt[, `:=`(G = NULL, VAR = NULL)]
      setkeyv(pop_dt, c("CHR", "CM", "POS", "SNP"))

      pops <- pop_dt[, .GRP, POP]

      fst_dt <- 
        pop_dt[POP == pops[1, POP]
               ][pop_dt[POP == pops[2, POP]], 
                 on = c("CHR", "CM", "POS", "SNP")]

      fst_dt[, `:=`(T1 = (AF - i.AF) ^ 2 + 
                         AF * (1 - AF) / NCHROBS + 
                         i.AF * (1 - i.AF) / i.NCHROBS,
                    T2 = AF * (1 - i.AF) + i.AF * (1 - AF)
                    )]

      fst_dt[, c('NCHROBS', 'POP', 'AF', 'i.NCHROBS', 'i.POP', 'i.AF') := NULL]
      fst_dt[, FST := T1 / T2]
      fst_dt[ (FST < 0) | is.na(FST) , FST := 0] 

      return(fst_dt[])
    }else{
      print('Hudson Fst is evalalueter
    }
