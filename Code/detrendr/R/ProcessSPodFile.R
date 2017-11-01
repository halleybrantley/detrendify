###############################################################################
# Process SPod File
# Halley Brantley
# 07-13-2017

# 07-21-2017 JRC added sonic.w and GPS
# 8-30-2017 HLB added standard devations to output, adjusted to handle 5 nodes
# 8-30-2017 HLB changed output format for Tai
# 10-25-2017 HLB updated baseline functions, changed 10sec data to long format
###############################################################################
require(plyr)
require(tools)
require(zoo)
require(devtools)
require(Matrix)
require(Cairo)
require(ggplot2)
# Load functions from file
source("./SPodFunctions/screeningFunctions.R")
# Load package - for "getBaseline" and related function
find_rtools()
load_all("C:/Users/HBrantle/Desktop/Repositories/detrendify/Code/detrendr")
###############################################################################
# Need to adjust!
nodes <- c("f", "g", "h")
latitude <- c(30, 31, 31) # Enter GPS Latitude for each node
longitude <- c(-70, -72.2, -70) # Enter GPS Longitude for each node
sites <- c(1, 1, 1)
dataFolder <- "Raw_Data/Data_ShedRow/fgh"
# dateI <- "2017-08-23"
# filename <- paste0(dataFolder,"/", "PreDeployment-all_", dateI, ".csv")
files <- dir(dataFolder, full.names=TRUE)
pb <- txtProgressBar(min = 0, max=length(files))
for (filename in files){
  print(paste0("Processing file: ", filename))
  dateI <- substr(basename(filename), 15, 24)
  ###############################################################################
  # Load SPod file
  spod <- read.csv(filename, 
                   header=TRUE, na.strings = "N/A")
  spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                   format= "%m/%d/%Y %H:%M:%S")) 
  ###############################################################################
  # QA screens for off-scale and relative humidity swings
  
  # Create new data.frame to store processed data
  spodQC <- spod
  
  for (node in nodes){
    pid <- paste(node, "SPOD.PID..V.", sep=".")
    if (class(spod[,pid]) != "numeric"){
      spod[,pid] <- as.numeric(as.character(spod[,pid]))
    }
    rh <- paste(node,  "SPOD.Pressure..V.", sep=".")
    if (class(spod[,rh]) != "numeric"){
      spod[,rh] <- as.numeric(as.character(spod[,rh]))
    }
    if(length(which(!is.na(spodQC[,pid])))>0){
      spodQC[,pid] <- screenOffScale(spod[,pid], spod$time)
      if(length(which(!is.na(spodQC[,pid])))>0){
        spodQC[,pid] <- screenRH(spodQC[,pid], spod$time, spod[,rh])
      }
    }
    
    ########################################################
    uCol <- paste(node, "SPOD.Sonic.U", sep=".")
    vCol <- paste(node, "SPOD.Sonic.V", sep=".")
    wCol <- paste(node, "SPOD.Sonic.W", sep=".")
    spodQC[which(spodQC[, uCol]>50), uCol] <- NA
    spodQC[which(spodQC[, vCol]>50), vCol] <- NA
    spodQC[which(spodQC[, wCol]>50), wCol] <- NA
  }
  
  ###############################################################################
  # Plot screened time series
  time <- spod$time
  
  Cairo(file= paste0("./Daily_Figures/", dateI,
                    "Screened_timeSeries.png"),
        type="png",
        units = "in",
        width=4, height=6,
        pointsize = 12,
        dpi = 400)
  par(mfrow=c(length(nodes),1), mar=c(2,4,2,2))
  for (node in nodes){
    pid <- paste(node, "SPOD.PID..V.", sep=".")
    if(length(which(!is.na(spod[,pid])))>0){
      plot(spod[,pid]~time, type="l", ylab=paste(node, "PID",  "cts"))
      lines(spodQC[,pid]~time, col="red")
    }
  }
  dev.off()
  ###############################################################################
  # Aggregate and remove baseline
  timeBase <- "10 sec"
  timeBreaks <- seq(round(min(time), "min"), 
                    round(max(time), "min"), timeBase)
  spodQC$timeCut <- cut(spodQC$time, timeBreaks)
  par(mfrow=c(1,1))
  
  spodAgg <- {}
  i <- 1
  for (node in nodes){
    pidCol <- paste(node, "SPOD.PID..V.", sep=".")
    uCol <- paste(node, "SPOD.Sonic.U", sep=".")
    vCol <- paste(node, "SPOD.Sonic.V", sep=".")
    wCol <- paste(node, "SPOD.Sonic.W", sep=".")
    sTempCol <- paste(node, "SPOD.Sonic.Temp..oC.", sep=".")
    spodNode <- spodQC[, c("timeCut", pidCol, uCol, vCol, wCol, sTempCol)]
    names(spodNode)[2:6] <- c("pid", "u", "v", "w", "sonictemp")
    spodAggNode <- ddply(subset(spodNode, !is.na(timeCut)), 
                         .(timeCut), summarize, 
                     pid = mean(pid, na.rm=TRUE),
                     sonic.u = mean(u, na.rm=TRUE),
                     sonic.v = mean(v, na.rm=TRUE),
                     sonic.w = mean(w, na.rm=TRUE),
                     sonic.sonictemp = mean(sonictemp, na.rm=TRUE), 
                     meanSquares = mean(pid^2, na.rm=TRUE))

    spodAggNode$latitude <- latitude[i]
    spodAggNode$longitude <- longitude[i]
    spodAggNode$site <- sites[i]
    spodAggNode$DeviceNumber <- node

    if(length(which(!is.na(spodAggNode[,"pid"])))>0){
      print(paste0("Fitting baseline node ", node))
      spodAggNode[, "baseline"] <- 
        getBaseline(spodAggNode[, "pid"], 
                  lambda0 = 10, 
                  maxiter = 10000, tau=.1)
      spodAggNode[, "pid.corrected"] <- 
        spodAggNode[, "pid"] -
        spodAggNode[,  "baseline"]
    } else {
      spodAggNode$baseline <- NA
      spodAggNode$pid.corrected <- NA
    }
    
    spodAgg <- rbind(spodAgg, spodAggNode)
    i <- i+1
  }

  spodAgg$time <- as.POSIXct(as.character(spodAgg$timeCut))
  spodAgg <- spodAgg[order(spodAgg$time),]

  write.csv(spodAgg[,c(14, 2:13)], paste0("./Processed_Data/Baseline_removed_", 
                            timeBase, "_avg_", 
                            dateI, ".csv"), row.names=FALSE)
  #############################################################################
  # Aggregate and detect signal
  timeBase <- "5 min"
  timeBreaks <- seq(round(min(time), "min"), 
                    round(max(time), "min"), timeBase)
  spodAgg$timeCut <- cut(spodAgg$time, timeBreaks)
  
  # calculate median of corrected.pid signal by node
  # for each node
  nodeMedian <- ddply(spodAgg, .(DeviceNumber), 
                   summarize,
                   pid.median = quantile(pid.corrected, .5, na.rm=T))
  spodAgg2 <- merge(spodAgg, nodeMedian, by="DeviceNumber")
  
  # Plot raw pid with baseline
  ggplot(spodAgg, aes(y=pid, x=time, group=DeviceNumber)) +
    geom_line() +
    geom_line(aes(y=baseline), col="red") +
    facet_grid(DeviceNumber~.)#, scales="free_y")
  ggsave(file = paste0("./Daily_Figures/", dateI,
                       "baseline_fit.png"), width=10, height=10)
  
  # Plot correcte pid with threshold set at 4 times the median
  ggplot(spodAgg2, aes(y=pid.corrected, x=time, group=DeviceNumber)) +
    geom_line() +
    geom_line(aes(y=4*pid.median), col="red")+
    facet_grid(DeviceNumber~.)
  ################################################################################
  detectSignal <- function(pid, threshold){
    length(which(pid > threshold))
  }
  sdFunc <- function(y, q){
    sd(y[y<q], na.rm=TRUE)
  }
  ###############################################################################
  # Aggregate and count detects
  timeBase <- "5 min"
  timeBreaks <- seq(round(min(time), "min"), 
                    round(max(time), "min"), timeBase)
  spodAgg$timeCut <- cut(spodAgg$time, timeBreaks)
  
  # calculate standard deviation of corrected pid below 4 times median
  # for each node
  noiseSD <- ddply(spodAgg2, .(DeviceNumber), 
                   summarize,
                   pid.sd = sdFunc(pid.corrected, 4*pid.median))
  
  # Aggregate data, set detects to 5 times noise standard deviation
  spodAgg3 <- merge(spodAgg2, noiseSD)
  spodSummary <- ddply(spodAgg3, .(timeCut, DeviceNumber), summarize, 
                           pid.mean = mean(pid.corrected, na.rm=TRUE),
                           possibleDetects = length(which(!is.na(pid.corrected))),
                           detects = detectSignal(pid.corrected, 6*pid.sd),
                           sonic.u = mean(sonic.u, na.rm=TRUE),
                           sonic.v = mean(sonic.v, na.rm=TRUE),
                           sonic.w = mean(sonic.w, na.rm=TRUE),
                           sonic.sonictemp = mean(sonic.sonictemp, na.rm=TRUE),
                           pid.sd = sd(pid.corrected, na.rm=TRUE),
                           pid.sdRaw = sqrt(mean(meanSquares, na.rm=TRUE) - 
                                              mean(pid)^2), 
                           latitude = mean(latitude), 
                           longitude = mean(longitude),
                           site = sites[1])
  
  spodSummary$wd <- atan2(-spodSummary$sonic.u, -spodSummary$sonic.v)*180/pi + 180
  spodSummary$ws <- sqrt(spodSummary$sonic.u^2 + spodSummary$sonic.v^2)
  spodSummary[which(spodSummary$detects>0), c("pid.sd", "pid.sdRaw")] <- NA
  
  spodSummary <- spodSummary[which(!is.na(spodSummary$timeCut)),
                             c("site", "timeCut", "DeviceNumber", "pid.mean", 
                              "possibleDetects", "detects", "sonic.u",
                              "sonic.v",	"sonic.w",	"sonic.sonictemp",
                              "pid.sd",	"pid.sdRaw",	"wd",	"ws",	"latitude",	
                              "longitude")]
  
  spodSummary$status <- 1
  spodSummary$status[which(is.na(spodSummary$pid.mean))] <- 0
  spodSummary$QAQC <- NA
  
  ddply(spodSummary, .(DeviceNumber), summarize,
        mean.sd = mean(pid.sd, na.rm=TRUE), 
        mean.RawSD = mean(pid.sdRaw, na.rm=TRUE))
  
  spodSummary$timeCut <- as.POSIXct(as.character(spodSummary$timeCut))
  
  # Plot detects by node
  ggplot(spodSummary, aes(x=timeCut, y=detects,
                          group=DeviceNumber, col=DeviceNumber)) +
    facet_grid(DeviceNumber~.)+
    geom_line()
  ggsave(file = paste0("./Daily_Figures/", dateI,
                       "_detects.png"), width=10, height=10)
  
  spodSummary <- spodSummary[order(spodSummary$timeCut),]
  
  write.csv(spodSummary, paste0("./Processed_Data/SPod_summary_", 
                            timeBase, "_avg_", 
                            dateI, ".csv"), row.names=FALSE)
}
