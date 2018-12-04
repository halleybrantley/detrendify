library(splines)
rm(list=ls())
source("data_generating_functions.R")

set.seed(39207491)
for(n in c(300,500,1000)){
  for (i in 1:100){
    df <- generate_mixednorm(n)
    save(df, file = sprintf("../SimData/mixednorm_n_%i_sim%03.0f.RData", n, i))
    df <- generate_gaus(n)
    save(df, file = sprintf("../SimData/gaus_n_%i_sim%03.0f.RData", n, i))
    df <- generate_shapebeta(n)
    save(df, file = sprintf("../SimData/shapebeta_n_%i_sim%03.0f.RData", n, i))
    df <- generate_peaks_design(n)
    save(df, file = sprintf("../SimData/peaks_n_%i_sim%03.0f.RData", n, i))
  }
}
