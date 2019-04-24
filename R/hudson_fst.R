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

      pop1_dt <- pop_list[[1]]
      pop2_dt <- pop_list[[2]]


      on_cols = c("CHR", "CM", "POS", "VAR", "SNP")
      scols = c("AF", "NCHROBS", "POP", on_cols)

      fst_dt <- merge(pop1_dt[, ..scols], 
                      pop2_dt[, ..scols], 
                      by = on_cols,
                      all = TRUE)

      fst_dt[, `:=`(T1 = (AF.y - AF.x) ^ 2  
                        - AF.y * (1 - AF.y) / (NCHROBS.y - 1)
                        - AF.x * (1 - AF.x) / (NCHROBS.x - 1),
                    T2 = AF.y * (1 - AF.x) + AF.x * (1 - AF.y)
                    )]

      fst_dt[, FST := T1 / T2]
      fst_dt[ (FST < 0) | is.na(FST) , FST := 0] 

      return(fst_dt[, c("CHR", "CM", "POS", "SNP", "T1", "T2", "FST"), with = FALSE])
    }else{
      print('Hudson Fst is evalueted for 2 populations')
    }
    }
