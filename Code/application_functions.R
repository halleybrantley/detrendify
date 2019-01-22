get_spod_signal <- function(tau, spod_trends, spodPIDs, crit=4){
  spodPeaks <- select(spodPIDs, -time) - select(spod_trends, ends_with(paste(tau)))
  spodPeaks[is.na(spodPeaks)] <- 0
  thresholds <- apply(spodPeaks, 2, 
                      function(x) median(x, na.rm=T) + 
                        crit*mean(abs(x-median(x, na.rm=T)), na.rm=T))
  spodSignal <- spodPeaks
  for (i in 1:length(thresholds)){
    spodSignal[,i] <- as.numeric(spodPeaks[,i]>thresholds[i])
  }
  return(spodSignal)
}

get_confusion <- function(spod_signal){
  spod_signal <- na.omit(spod_signal)
  for (i in 1:ncol(spod_signal)){
    spod_signal[,i] <- factor(spod_signal[,i])
  }
  mat_h1 <- confusionMatrix(spod_signal$f[spod_signal$h==1],
                            spod_signal$g[spod_signal$h==1])
  
  mat_h0 <- confusionMatrix(spod_signal$f[spod_signal$h==0],
                            spod_signal$g[spod_signal$h==0])
  
  conf_out <- cbind(mat_h0$table, mat_h1$table)
  return(conf_out)
}