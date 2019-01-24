require(mcclust)
require(aricode)

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
  
  conf_out <- data.frame(t(c(as.numeric(mat_h0$table), as.numeric(mat_h1$table))))
  names(conf_out) <- c("h0f0g0", "h0f1g0", 
                       "h0f0g1", "h0f1g1", 
                       "h1f0g0", "h1f1g0", 
                       "h1f0g1", "h1f1g1")

  return(conf_out)
}

get_NMI <- function(signal){
  data.frame(fg = round(NMI(signal$f, signal$g, variant="sqrt"),2),
             fh = round(NMI(signal$f, signal$h, variant="sqrt"),2), 
             gh = round(NMI(signal$f, signal$g, variant="sqrt"),2))
}

get_VI <- function(signal){
  data.frame(fg = round(vi.dist(signal$f, signal$g),2),
             fh = round(vi.dist(signal$f, signal$h),2), 
             gh = round(vi.dist(signal$f, signal$g),2))
}
