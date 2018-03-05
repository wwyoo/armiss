ar.miss <-
function(data, epsilon = 0.001, order = NULL,
             max.iter = 100, sym = NA, control.optim = list(maxit = 200))
#Imputation via conditional Gaussian

#argument: data
#input: vector or time series object with numeric entries

#argument: epsilon
#input: a numeric value that controls convergence criterion
#       defaults to 0.001

#argument: order
#input: positive integer for AR order(if known)
#       defaults to NULL and use AIC to select

#argument: max.iter
#input: positive integer to control the max iteration
#       of the proposed algorithm

#argument: sym
#input: symbol used to indicate missing values,
#       defaults to NA, but can be numeric or character

#argument: contol.optim
#input: list of control options to do MLE using optim
#       set max iter for BFGS to be 200

#output: estimates of the process mean, AR coefficients, and
#        innovation variance based on MLE
{
  require(mvtnorm) 
  #library to construct multivariate normal likelihoods
  #and to sample from multivariate normal
  #case where no missing values
  if((any(is.na(data)) || any(data == sym)) == FALSE){
    cat("no missing values detected...\n")
    cat("estimate AR parameters using complete data...\n")

    autocom <- ar.mle(data, aic = TRUE, demean = TRUE)
    thetahat <- c(mu = autocom$x.mean, phi = autocom$ar,
                    sigma2 = autocom$var.pred)
    return(thetahat)
  }
  else{
    N <- length(data)  #data length 
    P <- elem(data, sym = sym)  #permutation matrix
    if(is.na(sym)){
       XO <- data[!is.na(data)]  #observed data
    }
    else{
       XO <- data[data!=sym]
    }
    tobs <- length(XO)
    PO <- P[1:tobs, ]

    #initial estimates obtained by fitting AR(1)
    #using Yule-Walker
    if(is.null(order) == TRUE){ #if no order specified
      init <- ar.yw(XO, aic = FALSE, order.max = 1)
    }
    else{
      init <- ar.yw(XO, aic = FALSE, order.max = order)
    }
    mu <- init$x.mean  #process mean
    phihat <- init$ar  #AR coefficients
    sigma2hat <- init$var.pred  #innovation variance

    #imputation step
    #algorithm stop if :
    #1. satisfy convergence criterion
    #2. exceed upper iteration limit
    #whichever of the two occur first
    for(i in 1:max.iter){
      Pmu <- mu * rep(1, N)
      muO <- Pmu[1:tobs]
      muM <- Pmu[-(1:tobs)]
      SigmaPP <- P %*% covmat(phi = phihat,
                   sigma2 = sigma2hat, N = N) %*% t(P)
      SigmaOO <- SigmaPP[1:tobs, 1:tobs]
      SigmaMM <- SigmaPP[(tobs + 1):N, (tobs + 1):N]
      SigmaOM <- SigmaPP[1:tobs, (tobs + 1):N]
      SigmaMO <- SigmaPP[(tobs + 1):N, 1:tobs]
      
      #likelihood at previous iterate
      L1 <- dmvnorm(XO, mean = muO, sigma = SigmaOO, log = TRUE)

      muMgen <- muM + SigmaMO %*% solve(SigmaOO, (XO - muO))
      SigmaMgen <- SigmaMM - SigmaMO %*% solve(SigmaOO, SigmaOM)
      
      #generate missing values
      YOplus <- rmvnorm(1, muMgen, SigmaMgen)
      Yplus <- c(XO, YOplus)  #combine with observed values
      Y <- solve(P, Yplus)  #reconstruct series
 
      #estimate AR coefficients and innovation variance using MLE
      #choose by AIC if no order specified
      if(is.null(order) == TRUE){
        auto <- ar.mle(Y, aic = TRUE, demean = TRUE,
                  control = control.optim)
      }
      else{
        auto <- ar.mle(Y, aic = FALSE, demean = TRUE,
                  order.max = order, control = control.optim)
      }

      #reestimate parameters based on reconstructed series
      mu <- auto$x.mean
      phihat <- auto$ar
      sigma2hat <- auto$var.pred
      muO <- mu * rep(1, tobs)
      SigmaOO <- PO %*% covmat(phi = phihat,
                   sigma2 = sigma2hat, N = N) %*% t(PO)
      #likelihood at current iterate
      L2 <- dmvnorm(XO, mean = muO, sigma = SigmaOO, log = TRUE)
      
      #stopping criterion with tolerance level
      normratio <- abs(L2 - L1)
      if(normratio < epsilon * abs(L1)){
        break  #exit loop if satisfy criterion
      }
    }
    thetahat <- c(mu = mu, phi = phihat, sigma2 = sigma2hat)
    return(thetahat)  #converged final estimates
  }
}
