require(mcclust)
require(aricode)

get_spod_signal <- function(tau, trends, spodPIDs, crit=0.9){
  spodPeaks <- select(spodPIDs, -time) - select(trends, ends_with(paste(tau)))
  thresholds <- apply(spodPeaks, 2, quantile, crit, na.rm=T)
  spodSignal <- spodPeaks
  for (i in 1:length(thresholds)){
    spodSignal[,i] <- as.numeric(spodPeaks[,i]>thresholds[i])
  }
  return(spodSignal)
}

get_confusion <- function(spod_signal, nodes){
  spod_signal <- na.omit(spod_signal)
  for (i in 1:ncol(spod_signal)){
    spod_signal[,i] <- factor(spod_signal[,i])
  }
  if (length(nodes) == 3){
  mat_0 <- confusionMatrix(spod_signal[spod_signal[, nodes[1]]==0, nodes[2]],
                            spod_signal[spod_signal[, nodes[1]]==0, nodes[3]])
  
  mat_1 <- confusionMatrix(spod_signal[spod_signal[, nodes[1]]==1, nodes[2]],
                            spod_signal[spod_signal[, nodes[1]]==1, nodes[3]])
  
  conf_out <- data.frame(t(c(as.numeric(mat_0$table), as.numeric(mat_1$table))))
  names(conf_out) <- c("nodes000", "nodes010", 
                       "nodes001", "nodes011", 
                       "nodes100", "nodes110", 
                       "nodes101", "nodes111")
  } else {
   mat <- confusionMatrix(spod_signal[, nodes[1]],
                    spod_signal[, nodes[2]])
   conf_out <- data.frame(mat$table)
  }

  return(conf_out)
}

get_NMI <- function(signal, nodes){
  data.frame(nodes12 = round(NMI(signal[, nodes[1]], signal[,nodes[2]], variant="sqrt"),2),
             nodes13 = round(NMI(signal[, nodes[1]], signal[, nodes[3]], variant="sqrt"),2), 
             nodes23 = round(NMI(signal[, nodes[2]], signal[, nodes[3]], variant="sqrt"),2))
}

get_VI <- function(signal, nodes){
  if (length(nodes)==3){
  data.frame(nodes12 = round(vi.dist(signal[, nodes[1]], signal[,nodes[2]]),2),
             nodes13 = round(vi.dist(signal[, nodes[1]], signal[,nodes[3]]),2), 
             nodes23 = round(vi.dist(signal[, nodes[2]], signal[,nodes[3]]),2))
  } else {
   round(vi.dist(signal[, nodes[1]], signal[,nodes[2]]),2)
  }
}
  
getThresh <- function(x){mean(x, na.rm=T)+3*sd(x, na.rm=T)}

getCutoff <- function(peaks, node){
  thresh0 <- getThresh(peaks[,node])
  print(thresh0)
  peaks[which(peaks[,node] > thresh0), node] <- NA
  thresh <- getThresh(peaks[,node])
  while(thresh0 - thresh > 0.001){
    thresh0 <- thresh
    print(thresh0)
    peaks[which(peaks[,node] > thresh), node] <- NA
    thresh <- getThresh(peaks[,node])
  }
  thresh
}