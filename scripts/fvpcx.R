#--------------------------------------------------------------
#- Filename: fvpcx.R
#- Author: luke
#- Date: 9/24/20
#- Description: Applies FVPC forecast with exog vars to EML data
#--------------------------------------------------------------

#- [0] dependencies --------------------------------------------------
library(vars)
# library(fda)
# library(refund)
library(tidyverse)
# library(lubridate)
# library(viridis)
fpe_path <- "/Users/lukedurell/Box Sync/education/research/EML/aue_replication/FPE2.R"
source(fpe_path)
function_path <- "/Users/lukedurell/Box Sync/education/research/EML/scripts/EML_functions.R"
source(function_path)

EML_path <- "/Users/lukedurell/Box Sync/education/research/EML/data/clean_12_10_2019/EML_through_12_10_2019.csv"

#- [1] data --------------------------------------------------

#- read in and format data
EML <- read.csv(EML_path, stringsAsFactors = F)
EML$Date.Time <- as.POSIXct(EML$Date.Time, format = "%m/%d/%y %H:%M")

#- remove flats
flatdat <- flat_remove2(EML, "DOsat", 0)
flatdat <- flat_remove2(flatdat, "Temp", 0)

#- smooth data
# dosat_ps_fts <- penalized_smoother(flatdat, "DOsat", lam =  .01, nbasis= 20) #- lambda chosen from report 3
# temp_ps_fts  <- penalized_smoother(flatdat, "Temp", lam =  0.015, nbasis= 20) #- lambda chosen from report 3




data = flatdat
hat.d = 5
hat.p = 2
duration = 14
#- [2] forecaster fixed d and p --------------------------------------------------
fpcvx_forecaster <- function(data = flatdat, hat.d = 5, hat.p = 2, duration = 14) {
  #- This function takes in full data and fixed p and d
  #- it selects two weeks and then predicts one step ahead
  #- it is all wrapped into one function in order to be able to 
  #- use the true one-step-ahead value of the exogenous Temp variable.
  #- nbasis for penalized smoother = 20..
  
  #- select two weeks worth of data.
  
  #- #- first make sure to cut off end, ensure that start dat cant pick 
  #- close to end
  
  start_dat <- dplyr::filter(data, data$Date.Time < "2019-10-14 00:00:00")
  trunc_dat <- dplyr::filter(data, data$Date.Time < "2019-11-01 00:00:00")
  
  #- #- first get a random start date
  start <- sample(date(start_dat$Date.Time), 1)
  
  #- #- end date
  end <- start + duration + 1
  
  #- #- get the two-week data set and a two week + 2 hr set
  data_durat <- dplyr::filter(trunc_dat, trunc_dat$Date.Time >= start & trunc_dat$Date.Time < end)
  data_plus1 <- dplyr::filter(trunc_dat, trunc_dat$Date.Time >= start & trunc_dat$Date.Time <= end)
  
  if (dim(data_durat)[1] == dim(data_plus1)[1]) {
    message("end was flat, restart")
    return(NULL)
    }
  
  #- get data as centered FTS
  FTS_en <- penalized_smoother(data_durat, "DOsat", lam =  .01, nbasis= 20) #- lambda chosen from report 3
  FTS_ex <- penalized_smoother(data_plus1, "Temp", lam =  0.015, nbasis= 20) #- lambda chosen from report 3
  
  # #- Determines length of FTS
  # n_en <- dim(FTS_en[[1]])[2]
  # n_ex <- dim(FTS_ex[[1]])[2]
  
  #- before removing row means get the endogenous one step ahead pred!
  fcst_endog <- fpcv_forecast(FTS = FTS_en, hat.d = hat.d, hat.p = hat.p)
  
  #- Center data
  means_en <- rowMeans(FTS_en[[1]])
  means_ex <- rowMeans(FTS_ex[[1]])
  
  #- make sure this works:
  #- mean(FTS[[1]][1,])
  #- mean(FTS[[1]][2,])
  
  #- Subtract mean from all observations
  FTS_en[[1]] <- (FTS_en[[1]] - means_en)  	
  FTS_ex[[1]] <- (FTS_ex[[1]] - means_ex)  	
  
  #- #- the centered true one-step-ahead data. Notice, this always will forecast hour = 0
  true_step <- dplyr::filter(trunc_dat, trunc_dat$Date.Time == end)
  true_FTS_en  <- penalized_smoother(true_step, "DOsat", lam = 0.01, nbasis = 20)
  
  true_FTS_en[[1]]  <- true_FTS_en[[1]] - means_en
  
  #- FTS is the smoothed grid*n functional data set 
  #- hat.d and hat.p are results from method.FPE2
  
  #- now, predict n+1 score
  pca_FTS_en <- pca.fd(FTS_en, nharm=15) #- nharm = D normally
  pca_FTS_ex <- pca.fd(FTS_ex, nharm=15)
  
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
  
  # if(hat.d == 1) { #- NOT SET UP FOR FVPCX: waiting to fix the fFPE
  #   VAR.pre <- predict(arima(pca.FTS$scores[,1],order=c(hat.p,0,0)),n.ahead=1)$pred[1]
  #   yhat <- VAR.pre
  # } else { 
    
    #- get hat.d scores for the n+1 prediction!
    scores_en <- as.matrix(pca_FTS_en$scores[,1:hat.d])
    scores_ex <- as.matrix(pca_FTS_ex$scores[,1:dprime])
    last <- dim(scores_ex)[1]
    scores_ex_n <- as.matrix(scores_ex[-last,]) 
    scores_ex_1 <- t(as.matrix(scores_ex[last,]))

    # if(hat.p == 0) { #- NOT SET UP FOR FVPCX: waiting to fix the fFPE
    #   yhat <- colMeans(vec.ts)
    # } else { #- SET UP FOR D.HAT = 5 and P.HAT = 2
    
      #- we need to give colnames to the scores 
      #- to avoid warnings from vars predict function below
      colnames(scores_en) <- as.character(paste0("endog_", seq(1:hat.d)))
      colnames(scores_ex_n) <- as.character(paste0("exog_", seq(1:dprime)))
      colnames(scores_ex_1) <- as.character(paste0("exog_", seq(1:dprime)))
      
      #- forecast n+1 score with VAR model
      fit <- VAR(scores_en, exogen = scores_ex_n, p = hat.p)
      
      #- uses the true one step ahead exog!
      VAR.pre <- predict(fit, n.ahead=1, type="const", dumvar = scores_ex_1)$fcst
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

  list(fcst_exog, fcst_endog, true_FTS_en, start, end)
}


#- returns NULL if flat dat is selected
onestep <- fpcvx_forecaster(data = flatdat, hat.d = 5, hat.p = 2, duration = 14)
fcst_ex <- onestep[[1]] #- fts obj
fcst_en <- onestep[[2]] #- fts obj
real <- onestep[[3]] #- df obj
end <- onestep[[5]]
mse_ex <- inprod(real-fcst_ex, real-fcst_ex)
mse_en <- inprod(real-fcst_en, real-fcst_en)

plot(real, ylim = c(-80, 80))

title(paste0(end, "mse exog = ", round(mse_ex, 2), " / mse endog = ", round(mse_en, 2)))

lines(fcst_ex, col = "red")
lines(fcst_en, col = "blue")
