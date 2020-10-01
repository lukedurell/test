#--------------------------------------------------------------
#- Filename: fvpcx_EML.R
#- Author: luke
#- Date: 9/28/2020
#- Description: fvpc with exogenous variables
#--------------------------------------------------------------

#- [0] dependencies and functions --------------------------------------------------
source("/Users/lukedurell/Box Sync/education/research/EML/scripts/EML_functions.R")
#- new functions for exogenous FPEX
source("/Users/lukedurell/Box Sync/education/research/EML/fvpc_exog/FPEX.R")
library(vars)
library(tidyverse)


#- [1] data --------------------------------------------------
EML_path <- "/Users/lukedurell/Box Sync/education/research/EML/data/clean_12_10_2019/EML_through_12_10_2019.csv"

#- read in and format data
EML <- read.csv(EML_path, stringsAsFactors = F)
EML$Date.Time <- as.POSIXct(EML$Date.Time, format = "%m/%d/%y %H:%M")

#- remove flats
flatdat <- flat_remove2(EML, "DOsat", 0)
flatdat <- flat_remove2(flatdat, "Temp", 0)


#- [2] initialize FTS --------------------------------------------------
FTS_en <- penalized_smoother(flatdat, "DOsat", lam =  .01, nbasis= 21)
FTS_ex <- penalized_smoother(flatdat, "Temp", lam =  .015, nbasis= 21)


#- [3] dhat and phat --------------------------------------------------
dps <- method.FPEX(FTS_en = FTS_en, FTS_ex = FTS_ex, D = 21, Pmax = 3)
dps #- sadly for some reason this always returns dhat = D

#- even checking the little plot doesn't answer this question
fun <- function(d) {
  n = 1000
  p = 3
  R = 1
  (n + R + (p * d)) / (n - R - (p * d))
}
#- n ~ 100 for the one week, ~ 2000 for the whole data
ds <- 1:21
vals <- sapply(ds, fun)
plot(ds, vals)

#- [4] fix dhat = 5 phat = 2 and forecast one step ahead --------------------------------------------------
data <- flatdat
hat.d <- 5
hat.p <- 2
duration <- 14
  
#- vars::predict is just not working inside the function for some reason
#- uncomment forecast call when you figure out vars::predict
fvpcx_forecast <- function(data, hat.d, hat.p, duration) {
  #- fixed dhat = 5, phat = 2
  #- one step ahead from 2 weeks
  #- avoid: forecast day a flatline (empty) obs or forecast intounobserved month
  #-
  #- need to:
  #- use two weeks to center the endogenous variable, two weeks + 2 hour to center ex
  #- pca on the two weeks endog, pca on the two weeks + 2 hour exog
  #- take two weeks scores from both endog and exog, fit VAR
  start_dat <- dplyr::filter(data, data$Date.Time < "2019-10-14 00:00:00")
  trunc_dat <- dplyr::filter(data, data$Date.Time < "2019-11-01 00:00:00")

  start <- sample(date(start_dat$Date.Time), 1)
  end <- start + duration + 1
  
  #- #- get the two-week data set and a two week + 2 hr set
  data_n <- dplyr::filter(trunc_dat, trunc_dat$Date.Time >= start & trunc_dat$Date.Time < end)
  data_n1 <- dplyr::filter(trunc_dat, trunc_dat$Date.Time >= start & trunc_dat$Date.Time <= end)
  
  if (dim(data_n)[1] == dim(data_n1)[1]) {
    message("end was flat, restart")
    return(NULL)
  }
  
  #- get data as centered FTS
  FTS_en <- penalized_smoother(data_n, "DOsat", lam =  .01, nbasis= 21) #- lambda chosen from report 3
  FTS_ex <- penalized_smoother(data_n1, "Temp", lam =  0.015, nbasis= 21) #- lambda chosen from report 3
  
  #- before removing row means get the endogenous one step ahead pred!
  #- TO BE RETURNED
  fcst_endog <- fpcv_forecast(FTS = FTS_en, hat.d = hat.d, hat.p = hat.p)
  
  means_en <- rowMeans(FTS_en[[1]])
  means_ex <- rowMeans(FTS_ex[[1]])
  
  #- Subtract mean from all observations
  FTS_en[[1]] <- (FTS_en[[1]] - means_en)  	
  FTS_ex[[1]] <- (FTS_ex[[1]] - means_ex)	
  
  #- the centered true one-step-ahead data. Notice, this always will forecast hour = 0
  true_step   <- dplyr::filter(trunc_dat, trunc_dat$Date.Time == end)
  true_FTS_en <- penalized_smoother(true_step, "DOsat", lam = 0.01, nbasis = 21)
  
  #- TO BE RETURNED
  true_FTS_en[[1]]  <- true_FTS_en[[1]] - means_en
  
  #- now, predict n+1 score
  pca_FTS_en <- pca.fd(FTS_en, nharm=21) #- nharm = D normally
  pca_FTS_ex <- pca.fd(FTS_ex, nharm=21)
  
  #- decide d' based on proportion of variance 
  vars_ex <- pca_FTS_ex$varprop
  var_amount <- 0
  dprime <- 0
  for (i in 1:length(vars_ex)) {
    if (var_amount < .9) {
      dprime <- dprime + 1
      var_amount <- var_amount + vars_ex[i]
    }
  }

  #- get hat.d scores for the n+1 prediction!
  scores_en <- as.matrix(pca_FTS_en$scores[,1:hat.d])
  scores_ex <- as.matrix(pca_FTS_ex$scores[,1:dprime])
  last <- dim(scores_ex)[1]
  scores_ex_n <- as.matrix(scores_ex[-last,]) 
  scores_ex_1 <- t(as.matrix(scores_ex[last,]))

  #- we need to give colnames to the scores 
  #- to avoid warnings from vars predict function below
  colnames(scores_en)   <- as.character(paste0("endog_", seq(1:hat.d)))
  colnames(scores_ex_n) <- as.character(paste0("exog_", seq(1:dprime)))
  colnames(scores_ex_1) <- as.character(paste0("exog_", seq(1:dprime)))
  
  #- forecast n+1 score with VAR model
  #- uses the true one step ahead exog!
  #VAR.pre <- predict(VAR(scores_en, exogen = scores_ex_n, p = hat.p), n.ahead=1, type="const", dumvar = scores_ex_1)$fcst
  fit <- VAR(scores_en, p = hat.p, type="const", exogen = scores_ex_n)
  VAR.pre <- predict(fit, n.ahead=1, dumvar = scores_ex_1)$fcst
  
  yhat <- c()
  for(i in 1:hat.d)
  {
    #- extract forecasted values
    yhat <- c(yhat, VAR.pre[[i]][1])
  }

  #- initialize and sum up KL representation
  #- harmonics are vectors and yhat are scores
  FAR.pre <- pca_FTS_en$harmonics[1]*0
  for(i in 1:hat.d)
  {
    FAR.pre <- FAR.pre+yhat[i]*pca_FTS_en$harmonics[i]
  }
  fcst_exog <- FAR.pre

  #- uncomment below when you can use vars::predict inside function
  list(fcst_exog, fcst_endog, true_FTS_en, start, end)
}
data <- flatdat
hat.d <- 5
hat.p <- 2
duration <- 14

fvpcx_forecast <- function(data, hat.d, hat.p, duration) 
fvpcx_forecast(flatdat, 5, 2, 14)
#- [5] plot all that! --------------------------------------------------
real <- true_FTS_en

mse_ex <- inprod(real - fcst_exog, real - fcst_exog)
mse_en <- inprod(real - fcst_endog, real - fcst_endog)
  
plot(real, ylim = c(-80, 80))

title(end)

lines(fcst_exog, col = "red")
lines(fcst_endog, col = "blue")

legend("topleft", cex = .25, legend = c("FTS", paste0("FPEX: ", round(mse_ex,2)), paste0("FPE: ", round(mse_en,2))),
       col = c("black", "red", "blue"), lty = 1)


  
  
#   
# onestep <- fvpcx_forecast(data = flatdat, hat.d = 5, hat.p = 2, duration = 14)
# fcst_ex <- onestep[[1]] #- fts obj
# fcst_en <- onestep[[2]] #- fts obj
# real <- onestep[[3]] #- df obj
# end <- onestep[[5]]
# mse_ex <- inprod(real-fcst_ex, real-fcst_ex)
# mse_en <- inprod(real-fcst_en, real-fcst_en)
# 
# plot(real, ylim = c(-80, 80))
# 
# title(paste0(end, "mse exog = ", round(mse_ex, 2), " / mse endog = ", round(mse_en, 2)))
# 
# lines(fcst_ex, col = "red")
# lines(fcst_en, col = "blue")
# 
# 
# 
