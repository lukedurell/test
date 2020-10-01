#--------------------------------------------------------------
#- Filename: FPEX.R
#- Author: luke
#- Date: 9/28
#- Description: functions for fFPEX 
#--------------------------------------------------------------

#- [0] dependencies --------------------------------------------------
library(vars)
library(fda)

#- [1] trace function --------------------------------------------------

#- with exogenous variables!
FPEX.trace <- function(res, R = 2, p = 2){
  res <- t(t(res))
  d <- length(res[1,])
  n <- length(res[,1])
  if(d == 1)
  {
    out <- (p * d + n + R)/(n - R - p * d) * var(res)
  }
  else
  {
    out <- (p * d + n)/(n - p * d) * sum(diag(cov(res)))
  }
  out
}

#- [2] method for dhat and phat --------------------------------------------------
method.FPEX <- function(FTS_en, FTS_ex, D = 21, Pmax = 3){
  #- FTS is the smoothed grid*n functional data set 
  #- D is the max dimensions to use for the score prediction
  #- Pmax is the max dimensions to use for the VAR(p) model
  
  #- [2a] prep data, centering --------------------------------------------------

  #- Computes the mean from the first m observations
  means_en <- rowMeans(FTS_en[[1]])
  means_ex <- rowMeans(FTS_ex[[1]])
  
  #- Subtract mean from all observations
  FTS_en[[1]] <- (FTS_en[[1]] - means_en)
  FTS_ex[[1]] <- (FTS_ex[[1]] - means_ex)
  
  #- Performs FPCA on centered FTS data
  pca_FTS_en <- pca.fd(FTS_en, nharm=D)
  pca_FTS_ex <- pca.fd(FTS_ex, nharm=D)
  
  # compute the total variance
  vartot <- sum(pca_FTS_en$values)
  
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
  
  #- vector if dprime = 1 or matrix if dprime > 1
  scores_ex <- as.matrix(pca_FTS_ex$scores[,1:dprime])
  
  #- [2b] D*Pmax matrix of FPE values --------------------------------------------------
  #- the matrix below will contain the different (FPE + vartot - var.explain) values
  #- for p in 0:Pmax and d in 1:D
  
  #- D * Pmax matrix of FPE values, want min!
  values <- matrix(0, D, (Pmax+1))
  
  #- iterate over every dimension of max components!
  for(d in 1:D)
  {
    #- check out dim(pca_FTS_en$scores)
    scores_en <- as.matrix(pca_FTS_en$scores[,1:d])
    
    var.explain <- sum(pca_FTS_en$values[1:d])
    #- check out vartot - var.explain
    
    for(p in 0:Pmax)
    {
      if(d == 1)
      {
        #- this is just fitting an ar(p) and getting the residuals
        #- method = "ML" to avoid convergence errors...
        res <- arima(scores_en, xreg = scores_ex, order = c(p, 0, 0), method = "ML")$residuals
      }
      else
      {
        if(p == 0)
        {
          #- ar(0) model, just scores minus their mean, worthless b/c already centered
          #- noex_p0_mean  <- scores_en - t(matrix(rep(colMeans(scores_en), dim(scores_en)[1]), d))
          #- noex_p0_arima <- arima(scores_en, order = c(0,0,0), method = "ML")$residuals
          #- above two are same..
          
          #- below two are same (can check at d == 1)
          res <- lm(scores_en ~ scores_ex)$res
          # res <- arima(scores_en, xreg = scores_ex, order = c(0,0,0), method = "ML")$residuals
          
        }
        else
        {
          #- d > 1 and p > 0, we have a VAR model
          #- we have to give colnames to the scores, 
          #- otherwise we get warnings from vars resid function
          colnames(scores_en) <- as.character(seq(1:d))
          colnames(scores_ex) <- as.character(seq(1:dprime))
          
          res <- resid(VAR(scores_en, p=p, type="const", exogen = scores_ex))
        }
      }
      #- Fill in values matrix with FPE.trace result
      values[d, p+1] <- FPEX.trace(res=res, R = dprime, p=p) + vartot - var.explain
    }
  }
  
  # compute the estimates hat.p and hat.d for optimal order p and dimension d
  #- which.min goes column by column and finds the vector index of min val
  
  #- how many times does D go into which.min(values) - 1? gives column for p
  hat.p <- (which.min(values) - 1) %/% D
  
  #- remainder of D / which.min(values). If d = 21, you get 0 hence next if statement
  hat.d <- which.min(values) %% D
  if(hat.d == 0) {hat.d <- hat.d + D}
  c(hat.p, hat.d)
}

plot(1:21, values[,2], type = "l")
lines(1:21, values[,3], col = "red")
lines(1:21, values[,4], col = "blue")

df <- data.frame(p1 <- values[,2],
                 p2 <- values[,3],)







