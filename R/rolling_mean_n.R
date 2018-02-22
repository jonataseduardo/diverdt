
#' Evaluete the sliding window mean 
#'
#'\code{sw_mean()} estimates the branch size for every SNP 
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

sw_mean <-
  function(data_w, window_s = 3, step_s = 2){}
 data_w <- z

window_s = 4L
step_s = 2L

DT  <- data_w[c(1:11, (.N - 9):.N)]

DT <- data.table( CHR = c(rep(1,9),rep(2,10)), POS = c(sample(1:20, 9),sample(21:40,10)))
DT[, PBS := 1, by = CHR]
setkeyv(DT, c('CHR', 'POS'))
DT[, idx := .I]
DT[, idx := .I]
DT[, idd := .I]
DT[, next_idd := idd + window_s]

res <- DT[DT, .(idx = i.idx, seq = i.idx:idx),
          by = .EACHI,
          roll = Inf, 
          on = c(CHR = 'CHR', idd = 'next_idd')]

res

DT[, md := DT[res$seq, sum(PBS), by = res$idx]$V1]
DT

DT[, sid := seq(1:.N), by = CHR]
DT[, ss := sid %% step_s ]


res <- DT[DT, .(sid = i.sid, idx = i.idx, seq = i.idx:idx),
          by = .EACHI,
          roll = Inf, 
          on = c(CHR = 'CHR', idd = 'next_idd')]


res <- res[sid %% step_s == 1]

DT[sid %% step_s == 1, md := DT[, .(DT[res$seq, sum(PBS), by = res$idx]$V1)]$V1]

DT
