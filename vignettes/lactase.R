library(data.table)
library(ggplot2)

#library(diverdt)
devtools::load_all('~/diverdt')

pops <- c('EUR', 'EAS', 'AFR')

pop_list1 <- 
    lapply(
      pops, 
      function(pop_name){
        #fname <- paste0(system.file("extdata", package = "diverdt"), 
        #                '/', pop_name, '_Illumina')
        fname <- paste0('/home/jonatas/EUR_AFR_EAS/',
                         pop_name, '_Axion_Human_Origins')
        load_bim_frq(fname)
      }
    )

#pop_list <- lapply(pop_list1, maf_filter)
pop_list <- pop_list1
#pop_list <- pop_list[c(1,2)]

fst_eur_afr <- 
  wc_fst(pop_list[c(1,2)])

fst_eur_afr

pbs_data <-  
  pbs_fst(pop_list[[1]], pop_list[[2]], pop_list[[3]])

pbs_data[, P_RANK := 1 - frankv(PBS, ties = 'random') / (.N + 1)]
pbs_data[, LOG_PVAL := -log(P_RANK)]

#Print the best snp for each chr
pbs_data[pbs_data[, .I[ PBS == max(PBS)], by = CHR]$V1]

pbs_data[CHR == 2 & POS > 136545415 - 1E5 & POS < 136594750 + 1E5, .N]

manhattan_plot(pbs_data, col_name = 'LOG_PVAL', fig_name = 'all_illumina.png')

pbs_mean <- 
  rolling_mean_n(pbs_data, col_name = 'PBS', window_s = 20, step_s = 5)

pbs_mean[, P_RANK := 1 - frankv(PBS_MEAN) / (.N + 1)]
pbs_mean[, LOG_PVAL := - log(P_RANK)]
pbs_mean[, 1 / .N]

peaks_mean <- 
  get_peaks(pbs_data, pbs_mean, stat_col = 'PBS', 
            score_col = 'P_RANK', score_th = 0.01, 
            window_s = 20, greater = FALSE)

peaks_mean[CHR == 2]
pbs_mean[peaks_mean[CHR == 2], on = .(CHR, CM, W_ID)]

manhattan_plot(pbs_mean, col_name = 'LOG_PVAL', fig_name = 'mean_illumina.png')

pbs_median <- 
  rolling_median_n(pbs_data, col_name = 'PBS', window_s = 20, step_s = 5)

pbs_median[, P_RANK := 1 - frankv(PBS_MEDIAN, ties = 'random') / (.N + 1)]
pbs_median[, LOG_PVAL := -log(P_RANK)]

peaks_median <- 
  get_peaks(pbs_data, pbs_median, stat_col = 'PBS', 
            score_col = 'P_RANK', score_th = 0.01, 
            window_s = 20, greater = FALSE)

pbs_median[peaks_median[CHR == 2], on = .(CHR, CM, W_ID)]

pbs_median[, summary(LOG_PVAL)]

manhattan_plot(pbs_median, col_name = 'LOG_PVAL', fig_name = 'median_illumina.png')
