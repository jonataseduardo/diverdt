library(data.table)
library(RcppRoll)
library(GenABEL)

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

system.time( x <- wc_fst(pop_list))
system.time( z <- pbs_fst(pop_list[[1]], pop_list[[2]], pop_list[[3]]))
