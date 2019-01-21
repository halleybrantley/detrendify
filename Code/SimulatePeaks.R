library(splines)
rm(list=ls())
source("sim_generating_functions.R")

set.seed(39207498)
for(n in c(500,1000,2000,4000)){
  for (i in 1:100){
    df <- generate_peaks_design(n)
    save(df, file = sprintf("../SimData/peaks_n_%i_sim%03.0f.RData", n, i))
  }
}
