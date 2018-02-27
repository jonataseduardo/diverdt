
allele_fst <-
    function(ms_af){

    r <- length(ms_af[, .GRP, pop_id])
    aux <-
      merge(
        ms_af[VAR == 0], 
        ms_af[VAR == 0, .(nbar = mean(Ns),
                  nc = (sum(Ns) - sum(Ns ^ 2) / sum(Ns)) / (r - 1),
                  #Csqr =  (r / (r - 1)) * ( 1 - mean(Ns ^ 2) / mean(Ns) ^ 2),
                  pbar = mean(Ns * AF) / mean(Ns)),
                  by = .(rep_id, POS)],
        by = c("rep_id", "POS"))

    aux2 <-
      aux[, .(ssqr = sum(Ns * (pbar - AF) ^ 2) / ((r - 1) * nbar),
              pqbar = pbar * (1 - pbar), 
              nbar, 
              nc
               ),
          by = .(rep_id, POS)][, 
          `:=`(T1 = ssqr - (pqbar - (r - 1) * ssqr / r) / (nbar - 1), 
               T2 = (nc - 1) * pqbar / (nbar - 1) + 
                    (1 + (r -1) * (nbar - nc)/ (nbar - 1)) * ssqr / r)][,
           `:=`(al_FST = T1 / T2)]

    return(aux2[seq(1,nrow(aux2),r)])
    }
