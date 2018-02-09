library(data.table)

devtools::document('~/diverdt')
devtools::load_all('~/diverdt')

pops <- c('Xavante', 'Kayapo', 'Maya')

pop_list <- 
    lapply(
      pops, 
      function(pop_name){
        fname <- paste0('~/diverdt/inst/extdata/', pop_name)
        load_bim_frq(fname)
      }
    )

wc_fst <-
    function(pop_list){

    r <- length(pop_list)
    pop_dt <- rbindlist(pop_list)

    setkeyv(pop_dt, c("CHR", "CM", "POS", "VAR", "SNP"))

    pop_dt <- pop_dt[pop_dt[, .N , by = .(CHR, CM, POS, VAR, SNP)][N == r]]
    pop_dt[, N := NULL]
    pop_dt[, G := .GRP, by = .(CHR, CM, POS, VAR, SNP)]

    # Assuming biallelic data
    pop_dt <- pop_dt[G %% 2 != 0]
    pop_dt[, `:=`(G = NULL, VAR = NULL)]

    setkeyv(pop_dt, c("CHR", "CM", "POS", "SNP"))

    pop_dt

    nbar <- pop_dt[seq(1:r), mean(NCHROBS)]

    nc <- pop_dt[seq(1:r), 
                 (sum(NCHROBS) - sum(NCHROBS ^ 2) / sum(NCHROBS)) / (r - 1)]
    

    aux1 <- pop_dt[pop_dt[, .(nbar = nbar, 
                             nc = nc, 
                             pbar = mean(NCHROBS * AF) / nbar),
                  by = .(CHR, CM, POS, SNP)]]
    
    aux1[, G := .GRP, POP]

    aux2 <- 
      aux1[, .(ssqr = sum(NCHROBS * (pbar - AF) ^ 2) / ((r - 1) * nbar),
              pqbar = pbar * (1 - pbar), 
              nbar, 
              nc, 
              G
               ),
          by = .(CHR, CM, POS, SNP)
          ][G == 1]

    aux2[, G := NULL]


    fst_dt <- 
      aux2[,`:=`(T1 = ssqr - (pqbar - (r - 1) * ssqr / r) / (nbar - 1), 
                T2 = (nc - 1) * pqbar / (nbar - 1) + 
                     (1 + (r -1) * (nbar - nc)/ (nbar - 1)) * ssqr / r)
          ][,
           `:=`(FST = T1 / T2)]

    return(fst_dt)
    }

system.time( x <- wc_fst(pop_list))
