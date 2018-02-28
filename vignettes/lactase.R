library(data.table)
library(ggplot2)
library(diverdt)

pops <- c('EUR', 'EAS', 'AFR')

pop_list1 <- 
    lapply(
      pops, 
      function(pop_name){
        fname <- paste0(system.file("extdata", package = "diverdt"), '/', pop_name, '_Illumina')
        load_bim_frq(fname)
      }
    )

pop_list <- lapply(pop_list1, maf_filter)

fst_eur_afr <- 
  wc_fst(pop_list[c(1,2)])

pbs_data <-  
  pbs_fst(pop_list[[1]], pop_list[[2]], pop_list[[3]])

pbs_mean <- 
  rolling_mean_n(pbs_data, col_name = 'PBS', window_s = 20, step_s = 5)

pbs_mean[, P_RANK := frankv(PBS_MEAN) / (.N + 1)]
pbs_mean[, LOG_PVAL := -log(P_RANK)]

peaks <- get_peaks(pbs_data, pbs_mean, stat_col = 'PBS', 
                   score_col = 'LOG_PVAL', score_th = 0.01, 
                   window_s = 20, greater = FALSE)

manhattan_plot(pbs_mean, col_name = 'LOG_PVAL', fig_name = 'mean_illumina.png')

pbs_median <- 
  rolling_median_n(pbs_data, col_name = 'PBS', window_s = 20, step_s = 5)

pbs_median[, P_RANK := frankv(PBS_MEDIAN) / (.N + 1)]
pbs_median[, LOG_PVAL := -log(P_RANK)]

peaks_median <-  
  get_peaks(pbs_data, pbs_median, stat_col = 'PBS', 
            score_col = 'LOG_PVAL', score_th = 0.01, 
            window_s = 20, greater = FALSE)

manhattan_plot(pbs_median, col_name = 'LOG_PVAL', fig_name = 'median_illumina.png')
