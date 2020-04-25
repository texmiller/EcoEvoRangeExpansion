#-----------------------------------------------------------------------------#
# Mini meta on spread rates in exp-evo studies: functions
# Date created: Sept 26 2018
# When using lnCVR and sampling lnCVR codes, please cite: 
# Nakagawa et al. (2015) Methods in Ecology and Evolution, 6, 143-152.
#-----------------------------------------------------------------------------#


# 1. lnRR and sampling variance
#-----------------------------------------------------------------------------#
calclnRR <- function(mtreat, mcont){
  
  lnRR <- log(mtreat/mcont)
  
  return(lnRR)
  
}

calcsamplnRR <- function(SDtreat, SDcont, mtreat, mcont, ntreat, ncont){
  
  samplnRR <- ((SDcont^2)/(ncont * mcont^2)) + ((SDtreat^2)/(ntreat * mtreat^2))
  
  return(samplnRR)
  
}


# 2. lnCVR and sampling variance (cite: equations 11 and 12 from Nakagawa et al. 
# (2015) Methods in Ecology and Evolution, 6, 143-152)
#-----------------------------------------------------------------------------#
calclnCVR <- function(SDtreat, SDcont, mtreat, mcont, ntreat, ncont){
  
  lnCVR <- log((SDtreat/mtreat) / (SDcont/mcont)) + (1/(2 * (ntreat - 1))) - (1/(2 * (ncont - 1)))
  
  return(lnCVR)
  
}

calcsamplnCVR <- function(SDtreat, SDcont, mtreat, mcont, ntreat, ncont){
  
  samplnCVR <- (SDcont^2/(ncont * mcont^2)) + (1/(2 * (ncont - 1))) - (2 * (cor(log(mcont), log(SDcont))) * sqrt(((SDcont^2) / (ncont * mcont^2)) * (1 / (2 * (ncont - 1))))) + (SDtreat^2 / (ntreat * mtreat^2)) + (1 / (2 * (ntreat - 1))) - (2 * (cor(log(mtreat), log(SDtreat))) * sqrt(((SDtreat^2) / (ntreat * mtreat^2)) * (1 / (2 * (ntreat - 1)))))
  
  return(samplnCVR)
  
}


# 3. Calculate speed
#-----------------------------------------------------------------------------#

calcspeed <- function(dist, time){
  
  speed <- dist/time
  
  return(speed)
  
}

