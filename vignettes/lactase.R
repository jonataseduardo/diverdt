library(data.table)

devtools::document('~/diverdt')
devtools::load_all('~/diverdt')

pops <- c('EUR', 'AFR', 'EAS')

pop_list <- 
    lapply(
      pops, 
      function(pop_name){
        fname <- paste0('~/diverdt/inst/extdata/', pop_name, '_Illumina')
        load_bim_frq(fname)
      }
    )

pop_list <- lapply(pop_list, maf_filter)


system.time(
fst_eur_afr <- 
  wc_fst(pop_list[c(1,2)])
)



system.time( 
  pbs_data <-  
    pbs_fst(pop_list[[1]], pop_list[[2]], pop_list[[3]])
  )


pbs_mean <- 
  rolling_mean_n(pbs_data, col_name = 'PBS', window_s = 20, step_s = 5)

pbs_mean[, p_rank := frankv(PBS_MEAN) / (.N + 1)]
pbs_mean[, p.value := -log(p_rank)]

manhattan_plot(pbs_mean, fig_name = 'mean_illumina.png')

pbs_median <- 
  rolling_median_n(pbs_data, col_name = 'PBS', window_s = 20, step_s = 5)

pbs_median[, p_rank := frankv(PBS_MEDIAN) / (.N + 1)]
pbs_median[, p.value := -log(p_rank)]
manhattan_plot(pbs_median, fig_name = 'median_illumina.png')

pbs_moran <- 
  rolling_moran_n(pbs_data, col_name = 'PBS', window_s = 20, step_s = 5)

pbs_moran
manhattan_plot(pbs_moran, fig_name = 'moran_illumina.png')

