# 2024-06-09: This is Hanlin's code which produces interval forecasts for compositional data analysis and uses the alpha transformation
# 2024-06-13: Updated for LC CODA as the forecast method

#################
# load R package
#################

require(psych)
require(ftsa)
require(tseries)
require(sandwich)
require(Compositional)

###################################
# Function for computing the score
###################################

# l: lower bound
# u: upper bound
# x: actual holdout data
# alpha: level of significance alpha = 0.2

cpd <- function(holdout, lb, ub, alpha)
{
  lb_ind = ifelse(holdout < lb, 1, 0) #indicator for lower bound
  ub_ind = ifelse(holdout > ub, 1, 0) #indicator for upper bound
  cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout) #percentage coverage in the holdout data
  cpd = abs(cover - (1 - alpha)) 
  rm(lb_ind); rm(ub_ind); rm(cover)
  return(cpd)
}

interval_score <- function(holdout, lb, ub, alpha)
{
  lb_ind = ifelse(holdout < lb, 1, 0)
  ub_ind = ifelse(holdout > ub, 1, 0)
  score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  rm(lb_ind); rm(ub_ind)
  return(mean(score))
}

#####################
# Interval forecasts
#####################

# object: data matrix
# alpha_val: alpha tuning parameter
# ncomp_tuning: tuning parameter used for selecting the number of components
# fh: forecast horizon
# fore_method: forecast method
# B, K: number of bootstrap samples
# alpha: level of significance

# this interval forecast is only for the alpha transformation: need to redo this for CLR and ILR transformations too

ILR_fun_int <- function(object, ncomp_method, ncomp_tuning = 0.001, fh,
                          fore_method = c("ets", "arima", "rwf"), B = 399, K = 1, sig)
  
{
  object_ILR = ilr(object) #alpha used in alpha transformation, creates the alpha-transformed matrix of the data in the real space
  rownames(object_ILR) <- as.character(years.fit.1)
  n_year = nrow(object_ILR) #number of years = number of rows (ok)
  n_age = ncol(object_ILR) #number of ages - number of columns (not ok but only the label, this should be age bands * causes)

  # ================================================================
  SVD_decomp <- svd(object_ILR, nu=3, nv=3) #adding in the other factors for LC forecasting
  U <- SVD_decomp$u 
  V <- SVD_decomp$v
  S <- diag(SVD_decomp$d)
    
  bx<- V[,1]
  kt<- S[1,1]*U[,1] 
  bx2 <- V[,2]
  kt2 <- S[2,2]*U[,2]
  bx3 <- V[,3]
  kt3 <- S[3,3]*U[,3]
  
    # to determine the number of components used in principal component analysis using the ETS method (Shang and Haberman, 2019), the cumulative percentage of variance is used
    # selecting the number of components. Per Shang and Haberman (2019), when the number of components = 1 we get back to Lee Carter
    # for this application, use number of components = 1
    
    if(ncomp_method == "eigenvalue") 
  {
    ncomp = select_K(tau = ncomp_tuning, SVD_decomp$d^2) #don't do this for the alpha-transformation CODA paper (2024-06-13)
  }
  else if(ncomp_method == "fixed")
  {
    ncomp = 3 # here the original fixed components was 6, where Shang and Haberman (2019) tested using this for sensitivity analysis
  }
  else
  {
    warning("The number of components may be determined by eigenvalue ratio or fixed.")
  }
 
  basis = as.matrix(SVD_decomp$v[,1:ncomp]) #matrix of bx (which is the age / cause coefficient)
  fit = matrix(kt,length(years.fit.1),1) %*% t(bx) + matrix(kt2,length(years.fit.1),1) %*% t(bx2) + matrix(kt3,length(years.fit.1),1) %*% t(bx3)
  score_3 = cbind(kt, kt2, kt3)
  score = t(score_3)
  
  #check fit = score * basis
  fit_1 = t(score) %*% t(basis)
  # fit_1 should be the same as fit
  resi = as.matrix(t(object_ILR - fit))
  
  fh <- ih
  # forecast of the component scores from above
  score_fit_forecast = cbind(c(kt, kt.for.one), c(kt2, kt.for.two), c(kt3, kt.for.3))
  score_forecast = cbind(c(kt.for.one), c(kt.for.two), c(kt.for.3))
  
  # ================================================================
  # michelle's edits here (i.e. to replace this with the original forecast in the alpha-transformation script)
  
  # determine in-sample forecast error for principal component scores
  # olivia is the forecast matrix of means
  # replace with LC-CODA (or add the LC-CODA option in here)
  # the forecast under "arima" only takes the first component of the SVD (and we need 3)
  # do not need this - replace with the original forecast in the alpha-transformation coda script
  
  olivia = t(score_forecast)
  # check that alpha.proj.for.one[34:63,] is the same as t(basis %*% t(score_forecast))
  # the point forecast is equal to basis %*% olivia
  
  rownames(olivia) = 1:ncomp #only 1 row if using LC-CODA
  colnames(olivia) = 1:fh #forecast horizon, should be equal to ih in the alpha-transformation code
  
  # forecast errors
  # fore is the mean of the forecast of scores over the forecast horizon?
  # forerr is the mean of the forecast errors over the forecast horizon?
  # forerr = matrix(NA, (n_year - ncomp - fh + 1), ncomp) #matrix of forecast errors (this one is important for later, dimensions when using LC-CODA is the number of years you are forecasting for because ncomp = 1)
  
  forerr = matrix(NA, fh, ncomp)
  fore = matrix(NA, 1, ncomp) #forecast matrix 
  fore = colMeans(score_forecast)
  
  for(i in 1:fh)
  {
    forerr[i,] = score_forecast[i,] - fore # this is the error in the estimation of the principal component(s), subtracts the forecast mean of residuals from the actual mean of residuals
  }
  rm(i)
  
  B <- 399
  K <- 1
  
  # create the bootstrapped matrix of residuals
  q = array(NA, dim = c(n_age, B, K, fh)) #age/cause by number of bootstrap
  for(j in 1:fh)
  {
    for(i in 1:n_age)
    {
      for(k in 1:K)
      {
        q[i,,k,j] = sample(resi[i,], size = B, replace = TRUE) #resi is the matrix of residuals 
      }
    }
  }
  rm(i); rm(j); rm(k)
  
  #bootstrap sample of the forecast for 3 components
  ny = array(NA, dim = c(ncomp, B, fh))
  for(j in 1:fh)
  {
    for(i in 1:ncomp)
    {
      ny[i,,j] = sample(forerr[,i], size = B, replace = TRUE) #creates the matrix of PC score errors by sampling from forerr array
    }
  }
  rm(i); rm(j)
  
  # adding the PC score error to the predicted score
  # fh = forecast horizon
  
  # the combined matrix is called fo, which is the sum of oli (forecast matrix of means) and ny (PC score errors)
  oli = array(NA, dim = c(ncomp, B, fh))
  for(j in 1:fh)
  {
    for(i in 1:ncomp)
    {
      oli[i,,j] = olivia[i,j] #creates the matrix of PC score errors by sampling from forerr array
    }
  }
  rm(i); rm(j)
  
  #oli = array(rep(olivia, B * fh), dim = c(ncomp, B, fh))
  fo = array(NA, dim = c(ncomp, B, fh)) #3 by 399 by 30
  for(j in 1:fh)
  {
    for(i in 1:B)
    {
      fo[,i,j] = oli[,i,j] + ny[,i,j]
    }
  }
  rm(i); rm(j)
  
  # construct bootstrapped predictions
  # q is the residuals
  pred = array(NA, dim = c(n_age, B, K, fh))
  for(j in 1:fh)
  {
    for(i in 1:B)
    {
      for(k in 1:K)
      {
        pred[,i,k,j] = basis %*% fo[,i,j] + q[,i,k,j]
      }
    }
  }
  rm(i); rm(j); rm(k)
  
  pred_resize = array(NA, dim = c(n_age, B * K, fh))
  for(j in 1:fh)
  {
    for(i in 1:B)
    {
      pred_resize[,i,j] = pred[,i,,j]
    }
  }
  rm(i); rm(j)
  
  # transform back
  d_x_t_star_fore = array(NA, dim = c(fh, B * K, (n_age+1)))
  for(iw in 1:fh)
  {
    for(ij in 1:(B * K))
    {
      d_x_t_star_fore[iw,ij,] = ilrInv(t(pred_resize[,ij,iw])) #ilrInv because we applied the CLR transformation earlier
    }
  }
  rm(iw); rm(ij)
  
  return(apply(d_x_t_star_fore, c(1, 3), quantile, c((100 - sig)/200, 1 - (100 - sig)/200), na.rm = TRUE))
  
}