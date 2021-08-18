########################################################################################################
# Prerequisites - packages and functions
########################################################################################################

library(MASS)
library(survival)
library(bindata)

source("fn_simulcohort.R")
source("fn_samplingncc.R")
source("fn_simstudy.R")
source("fn_misc.R")


########################################################################################################
# Simulation Study - four censoring scenarios
########################################################################################################

for(lk in 1:4){
  ncc.design(nsimul = 2000, output1 = paste("res_C",lk*20,".rds", sep=""), 
             seed = 1, nfull = 5000, maxm = 5, 
             Smax.t = .95, Smax.censor = 1-(lk*.2), 
             relz1 = exp(.5), relz2 = exp(.9))
}


########################################################################################################
# Result Summary
########################################################################################################

summary_SD = summary_BIAS = summary_DV = summary_MEAN = summary_VAR = summary_RE = 
  array(NA, c(4,5,3,5), dimnames = list(cprop = c("C20%", "C40%", "C60%", "C80%"),
                                        mseq = c("m=1", "m=2", "m=3", "m=4", "m=5"),
                                        est = c("beta[1]", "beta[2]", "Lambda[0]"),
                                        mth = c("full", "standard/partial", "modified/partial", "standard/pseudo", "modified/pseudo")))


summary_VRE = summary_BRE = 
  array(NA, c(4,5,3,2), dimnames = list(cprop = c("C20%", "C40%", "C60%", "C80%"),
                                        mseq = c("m=1", "m=2", "m=3", "m=4", "m=5"),
                                        est = c("beta[1]", "beta[2]", "Lambda[0]"),
                                        mth = c("partial", "pseudo")))


for(lk in 1:4){
  res = readRDS(paste("res_C",lk*20,".rds",sep=""))
  
  for(mk in 1:5){
    tmp = summary.res(res[mk,,,])
    
    summary_BIAS[lk, mk, ,] = tmp$BIAS
    summary_MEAN[lk, mk, ,] = tmp$MEAN
    summary_VAR[lk, mk, ,] = tmp$VAR
    summary_DV[lk, mk, ,] = tmp$DV
    summary_SD[lk, mk, ,] = sqrt(tmp$VAR)
    summary_RE[lk, mk, ,] = tmp$RE
    
    ## Design Bias
    tmpB = abs(summary_MEAN[lk,mk,,-1]-summary_MEAN[lk,mk,,1]); 
    summary_BRE[lk,mk,,] = tmpB[,c(1,3)]/tmpB[,c(2,4)]
    
    ## Design Variance
    tmpV = summary_VAR[lk,mk,,-1]-summary_VAR[lk,mk,,1]; 
    summary_VRE[lk,mk,,] = tmpV[,c(1,3)]/tmpV[,c(2,4)]
  }
}
