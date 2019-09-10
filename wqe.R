# Copyright Information: This code was developed by Minh H. Nguyen at the Department of Computer 
# Science University College London in 2019. All queries are to be directed to minh.nguyen.18@ucl.co.uk
# or alternatively minguyen@tcd.ie.

library('TLMoments')
library('MASS')
library('dplyr')
library('pracma')
#library('PEIP')
library('dr')
library('ggplot2')
library('zoo')



################################ Maths functions ########################################

l2 <- function(x){return(sqrt(sum(x^2)))}

dot <- function(x,y){return( x%*%y )}

getFn <- function(data, fit = 'A', mquants = trunc(length(data)/4)) 
  {
    # Fit A is via a histogram of m boxes, fit B is distributed equidistance 
    # (A is overweight tail vs B is overweight body, see the paper for more details)
    
    ecdfdata <- ecdf(data)
    
    if (fit == 'A')
    { 
      X <- seq( min(data),max(data),length.out = mquants)
      
      Fn <- ecdfdata(X)
    }
    else if(fit == 'B')
    {
      Fn <- as.numeric(seq(from = 1/mquants , to = 1, length.out = mquants))
      X <- as.numeric(quantile(data, probs = Fn))
    }
    return(list(X = X, Fn = Fn, ecdf = ecdfdata))
  }

# This calcululates the bias of fitted quantile q compared to the sample quantile q_N 
# (in pct, i.e. 1 is no bias)
sambias <- function(data,iF,theta,q){ return( iF(q,theta)/quantile(data,q,names = FALSE) ) }

mse <- function(x,y){return(mean((x-y)^2))}

processdata <- function(data,fittype,mquants)
{
  new <- getFn(data,fit = fittype ,mquants = mquants)
  X <- new$X
  Fn <- new$Fn
  ecdf <- new$ecdf
  return(list(X = X,Fn = Fn,ecdf = ecdf))
}

pop_bias <- function(  c, results , q , Q_fitted ,  Q_true , F_true_param )   
{
  return ( fquantparams( Q_fitted, 1 - (1 -q)/(1-c), as.list(results$par) )/fquantparams(Q_true, 1 - (1 -q)/(1-c), as.list(F_true_param)  ))
}

trapzgof <- function( data, par )
{
  data <- sort(data)
  n <- length(data)
  dx <- 1/n
  
  Fn <- seq(dx,1,by = dx)
  
  Fx <- pgpd( data, loc = par[1], scale = par[2], shape = par[3] )
  
  return( trapz( Fn, abs(Fn - Fx) ) )
}


##################################### Monte Carlo VaR Functions ####################################

fnvarparams <- function(f,losses,params){return(do.call(f,c(losses,params)))}   

fncdfparams <- function(f,x,params){return(do.call(f,c(x,params)))}

fquantparams <- function( f, q, params ){ return(do.call(f,c(q,params))) }

fnvar <- function(data,cdf,fvar,period,threshold,trials,param,p)
  {
    if (trials %% 10 != 0) { trials <- round(trials , -1) } 
    
    blocktrials <- trials/10
    blocks <- rep(0,10)
    
    freq <- length(data)/period
    
    cdf_cutoff <- fncdfparams(cdf,threshold,as.list(param))
    
    if( cdf_cutoff >= p ) { return(NaN)}
    
    freqnormed <- freq/(1- cdf_cutoff)
    
    # Returns a vector of cumulative losses corresponding to the number of occurences per the elements of 
    # the vector occs. This eliminates an inner for loop.
    sumlosses <- function(occs) { return(sum(fnvarparams(fvar,occs,as.list(param)))) }
    
    for ( b in 1:10)
    {
      rand_occs <- rpois(blocktrials,freqnormed)
      trialslosses <- apply(matrix(rand_occs),MARGIN = 1, FUN = sumlosses)
      blocks[b] <- quantile(trialslosses,p)
    }
    
    return(round(mean(blocks)/1e6,10))
    
  }


############################ Numerical Analysis Functions #############################

givens <- function(a,b)
  {
    # Givens rotation in 2-D. Rotates the [a,b] vector by theta to zero out b.
    
    if ( b==0 ){costheta <- 1; sintheta <- 0}
    
    else
      {
        if (abs(b) >= abs(a))
          {
            cotantheta <- -a/b
            sintheta <- 1/sqrt(1 + cotantheta^2)
            costheta <- sintheta*cotantheta
          }
        else
          {
            tantheta <- -b/a
            costheta <- 1/sqrt(1 + tantheta^2)
            sintheta <- costheta*tantheta
          } 
      }
    return(list(cost = costheta,sint = sintheta))
  }

qrinsert <- function(R,lambda,i)
  {
    # To save computation time + space, this only takes in R_old and returns R_new, the necessary 
    # component for the direction. This returns R_1 from  A_1 = Q_1R_1 where A_1 contains A with 
    # row i of sqrt(lambda)*I_p appended on as last row of A. For more info. see van Loan 'Matrix Computations' 
    # 2014 pp. 335-338.
    
    m <- length(R[,1]) + 1
    n <- length(R[1,])
    R <- rbind(R, c(numeric(length(R[1,]))))
    R[m,i] <- lambda
    
    for (j in i:n)
      {
        rotated <- givens(R[j,j],R[m,j])
        cost <- rotated$cost
        sint <- rotated$sint
        
        for (k in j:n)
          {
            R_j <- R[j,k]
            R_m <- R[m,k]
            R[j,k] <- cost*R_j - sint*R_m
            R[m,k] <- sint*R_j + cost*R_m
          }
      }
    return(R)
  }

zoomInt <- function( f, df, alpha_low, alpha_up, c1, c2 )
  {
    # Zooming into region of step size where bisection is too close.
  
    tol <- 1e-6

    maxIter <- 100
    
    k <- 1
    
    while ( k <= maxIter )
      {
        alpha_k <- 0.5*( alpha_low + alpha_up )
        
        f_k <- f(alpha_k)
        
        if ( is.nan(f_k) ) { return(NaN) }
        
        if (abs(alpha_up - alpha_low) <= tol)
          {
            alpha <- alpha_k
            
            break
          }
        
        if ( f_k > f(0) + c1*alpha_k*df(0) | f_k >= f(alpha_low) )
          {
            alpha_up <- alpha_k
            
          }
        
        else 
          {
            df_k <- df(alpha_k)
            
            if (abs(df_k) <= -c2*df(0))
              {
                alpha <- alpha_k
                
                break
              }
            
            else if ( df_k*(alpha_up - alpha_low) >= 0 )
              {
                alpha_up <- alpha_low
              }
            
            alpha_low <- alpha_k
            
          }
          
          k <- k + 1
        
        }
    
    return(unname(alpha))
    
  }

linesearch <- function( f, df, x, dir)
{
  FACT <- 10
  c1 <- 1e-4
  c2 = 0.1
  maxIter <- 10
  
  alpha_max <- 1
  
  dir <- as.vector(t(dir))
  
  phi <- function( alpha ){ f(x + alpha*dir) }
  dphi <- function( alpha ){ dot( as.vector(t(df(x + alpha*dir))), dir) }
  
  info <-  matrix(0,nrow = maxIter + 2,ncol = 2)
  colnames(info) <- c('alphas','phis')
  
  alphas <- 0
  alphas <- c(alphas, 0.9*alpha_max)
  
  phis <- phi(0)
  
  if ( any(is.na(dir)) ) { return(NaN)}
  
  dphi0 <- dphi(0)
  
  alpha_k <- 0
  
  k <- 2
  
  while ( k <= maxIter )
  {
    phis <- c( phis , alphas[k] )
    
    dphi_k <- dphi(alpha_k)
    
    if ( phis[k] > phis[1] + c1*alphas[k]*dphi0 | 
         ( phis[k] >= phis[k-1] & k > 2 ) 
        )
      {
        alpha_k <- zoomInt( phi, dphi, alphas[k-1], alphas[k], c1, c2 )
        
        break
      } 
    
    else if ( abs( dphi_k ) <= -c2*dphi0 )
      {
        alpha_k <- alphas[k]
        
        break
      }
    
    else if ( dphi_k >= 0 )
      {
        alpha_k <- zoomInt( phi, dphi, alphas[k], alphas[k-1], c1, c2 )
        
        break
      }
    
    alphas <- c( alphas, min( FACT*alphas[k],alpha_max ) )
    
    k <- k + 1
    
  }
  
  return( alpha_k )
  
}


################################# Stage 1 ##############################

stage1 <- function(data,theta0,algo = 'LM',period = 6.5,cut_pct = 0,fittype = 'A',mquants = NULL,
                   rad0 = 500, stoptype = 'step',tol = 1e-10, maxIter = 1000,eta = 0.2, pvar = 0.999, 
                   printIter = FALSE, run_VaR = TRUE)
{
  # data:     Input a VECTOR/LIST of data.
  # theta0:   Initial guess for theta = [mu,sigma,xi] in that order.
  # algo:     LM for Levenberg-Marquadt or GN for Gauss-Newton.
  # period:   Length of the time span of the dataset.
  # cut_pct:  If the dataset is a truncated version of another, this is the percentage that has been 
  #           cut from the original.
  # fittype:  A returns equally spaced quantiles in the X space. B returns quantiles
  #           corresponding to equally spaced ECDF points.
  # mquants:  Number of quantiles to take (more quantiles yield more accurate results but takes longer,
  #           my experience indicates a quarter to a third of the length of the data set is sufficient.
  # rad0:     Initial maximum search radius. Just input a reasonable number (say 5 < rad0 < 5000) because 
  #           it is going to get altered by the maximum search algorithm.
  # stoptype: Convergence criterion. If one just wants a passing fit then 'gof' should do the trick. 
  #           Otherwise standard practice is by relative normalised step size or infinity norm of step size
  #           (see paper for more detail).
  # tol:      Convergence tolerance. Should be quite small if the residuals are small (i.e. the tail is GPD)
  # maxIter:  Maximum number of iterations. Note that one should be flexible with this if they wanted a 
  #           goof fit for a difficult problem.
  # eta:      Threshold for accuracy of the quadratic model. This is the ratio between actual decrease and 
  #           the predicted decrease through the JACobian approximation (see paper for more detail).
  # pvar:     Threshold percentile to calculate the MC VaR
  # prinIter: Set to TRUE for details of every 50th iteration.
  # run_VaR:  If TRUE then algorithm automatically runs an MC VaR using the fitted parameter values
  #           with 100,000 trials.
  
  frequency <- length(data)/period
  
  # If the data is truncated
  threshold <- unname(quantile(data,cut_pct))
  freauency_normed <- frequency/(1 - cut_pct)
  full_ecdf <- ecdf(data)
  data <- data[ data >= threshold]
  pvar <- 1 - 0.001/(1 - cut_pct)
  
  if ( is.null(mquants) ) { mquants <- trunc(length(data)/4) }

  else if (mquants > length(data)) { stop('Number of quantiles cannot be greater than length of the data set.') }
  
  dataset <- processdata(data,fittype ,mquants)
  X <- dataset$X
  Fn <- dataset$Fn
  ecdf <- dataset$ecdf
  frequency <- length(data)/period
  
  Fx <- function(theta) { 1 - (1 + (theta[3]/theta[2])*(X - theta[1]))^(-1/theta[3])  }
  
  # icdf or quantile function
  iFx <- function(q,theta){ return((theta[2]/theta[3])*((1 - q )^(-theta[3]) - 1) + theta[1] )}
  
  # Shorter version of the GoF test code. This is used to avoid for loops to boost convergence speed.
  trapzgof <- function(theta){return(trapz(Fn,res <- abs(Fn - Fx(theta))))}
  
  JAC <- function(theta) 
    { cbind((1/theta[2])*( 1 + (theta[3]/theta[2])*(X - theta[1]))^(-1/theta[3]-1)  ,
            ((X - theta[1])/theta[2]^2)*( 1 + (theta[3]/theta[2])*(X - theta[1]))^(-1/theta[3]-1) ,
            ((1/theta[3]^2)*log( 1 + (theta[3]/theta[2])*(X - theta[1])) - (X - theta[1])/(theta[2]*theta[3] +
            (theta[3]^2)*(X - theta[1])) )*( 1 + (theta[3]/theta[2])*(X - theta[1]))^(-1/theta[3]))
  }
  
  rho <- function(theta){return(Fn - Fx(theta))}
  
  phi <- function(theta){return(0.5*l2(rho(theta))^2)}
  
  dphi <- function(theta) { t(JAC(theta))%*%rho(theta) }
  
  lsqrdir <- function(theta_k,rad0,rho_k,J_k,maxIter = 100)
  {
    lambda <- 0.1
    QR0 <- qr(J_k)
    Q0 <-qr.Q(QR0)
    R0 <-qr.R(QR0)
    dtheta = numeric(length(theta_k))
    negdF <- -t(J_k)%*%rho_k
    dimQR <- dim(J_k); n <- dimQR[2];
    k <- 0
    while ( k <= maxIter & abs(l2(dtheta) - rad0) > 1e-12 & lambda > 0 )
    {
      Q_k <- Q0
      R_k <- R0
      for (i in 1:n){ R_k <- qrinsert(R_k,sqrt(lambda),i) }
      dtheta <- backsolve(R_k, backsolve(R_k,negdF,transpose = TRUE))
      eigvec <- backsolve(R_k,dtheta,transpose = TRUE)
      lambda <- max(0, lambda + ((l2(dtheta)/l2(eigvec))^2)*(l2(dtheta) - rad0)/rad0)
      k <- k + 1
    }
    if (lambda == 0)
    {
      dtheta <- backsolve(R_k, backsolve(R_k,negdF,transpose = TRUE))
    }
    return(dtheta)
  }
  
  # This grabs the approximate maximum radius depending on size of the problem (see paper, end of optimisation chapter)
  getradM <- function(theta,rad0,rho0,JAC0,maxIter = 100)
    {
      
      if ( any(is.na(JAC0))) 
      {
        # Fits spline inter/extrapolation for the Jacobian columns if there are NaNs. This is necessary for more difficult data sets.
        JAC0 <- apply( JAC0, MARGIN = 2, FUN = na.spline  )
      }
    
      dtheta <- lsqrdir(theta,rad0,rho0,JAC0,maxIter)
      
      dtheta_k <- dtheta
      
      radM <- rad0
      
      if (l2(dtheta_k) < 0.1*rad0){ radM <- 5*l2(dtheta_k)}
      while (abs(l2(dtheta_k) -radM) < 1e-5)
      {
        radM <- 5*radM
        dtheta_k <- lsqrdir(theta,rho0,radM,JAC0,maxIter)
      }
      return(radM)
    }
  
  fitLM <- function(theta0 = theta0,rad0 = rad0, stoptype = stoptype, tol = tol, maxIter = maxIter,eta = eta, printIter = printIter)
  {
    
    theta_k <- theta0
    radM <- getradM(theta0,rad0,rho(theta0),JAC(theta0))
    
    rad_k <- radM
    info <- matrix(0,nrow = 5,ncol = maxIter)
    rownames(info) <- c('gofIter','sambiasIter','tauIter','stepIter','phiIter')
    gofstat_k <- trapzgof(theta_k)
    k <- 1
    while (k <= maxIter)
    {
      
      J_k <- JAC(theta_k)
    
      if ( any(is.na(J_k))) 
      {
        J_k <- apply( J_k, MARGIN = 2, FUN = na.spline  )
      }
      
      JJ_k <- t(J_k)%*%J_k
      
      rho_k <- rho(theta_k)
      phi_k <- phi(theta_k)
      dphi_k <- t(J_k)%*%rho_k
      
      q_k <- function(delta){return(phi_k + dot(delta, dphi_k) + 0.5*dot(delta,(JJ_k%*%delta)))}
      
      dtheta_k <- as.vector(lsqrdir(theta_k,rad_k,rho_k,J_k))
      
      if (theta_k[2] + dtheta_k[2] < 0){dtheta_k[2] <- 0}
      
      if (theta_k[3] + dtheta_k[3] < 0){dtheta_k[3] <- 0}
      
      if ( X[1] - theta_k[1] - dtheta_k[1] <0 ) { dtheta_k[1] <- 0  }
      
      
      tau_k <- (phi_k - phi(theta_k + dtheta_k))/(q_k(numeric(length(dtheta_k))) - q_k(dtheta_k))
      
      if (is.nan(tau_k)) 
        { 
        rad_k <- 0.5*rad_k
        next 
        }
      
      if (tau_k < 0.25){rad_k <- 0.25*rad_k}
      
      else if((tau_k > 0.75) & (abs(l2(dtheta_k)^2 - rad_k^2) < 1e-12)){rad_k <- min(2*rad_k,radM)}
      
      theta_k_1 <- theta_k
      
      if (tau_k > eta)
        {
          theta_k <- theta_k + dtheta_k
          
          gofstat_k_1 <-gofstat_k
          
          gofstat_k <- trapzgof(theta_k)
          
          info['gofIter',k] <- gofstat_k
          info['sambiasIter',k] <- sambias(data,iFx,theta_k,pvar)
          info['tauIter',k] <- tau_k
          info['stepIter',k] <- l2(theta_k - theta_k_1)/l2(theta_k_1)
          info['phiIter',k] <- phi(theta_k)
          
          if (printIter == TRUE)
            {
              if (k%%50 == 0 )
                {
                cat(paste('Iteration: ',k,
                          '\nPhi_k ',round(info['phiIter',k],5),
                          '\nStep Size: ',round(info['stepIter',k],10),
                          '\nSample 99.9% Bias: ',round(info['sambiasIter',k],5),
                          '\nGoF Stat.: ',round(gofstat_k,5),'\n','\n'))
                }
            }
          
          if (stoptype == 'step')
            {
              if (info['stepIter',k] <= tol)
                {
                  break
                }
            }
          
          else if (stoptype == 'grad')
            {
              if (norm(dphi_k,'i') <= tol*(1 + abs(phi_k)))
                {
                  break
                }
            }
          
          else if (stoptype =='gof')
            {
              if (gofstat_k <= 0.068) 
                {
                  break
                }
            }
          
          if ( k > 100 & gofstat_k > gofstat_k_1 &  gofstat_k <= 0.068 ) 
            {
              break
            }
          
          if ( abs(info['sambiasIter',k] - 1) <= 0.05 &  gofstat_k <= 0.068 )
            {
              break
            }
          
          if (k > 10)
            {
              if (round(info['gofIter',k-10],8) - round(gofstat_k,8) == 0 & gofstat_k <= 0.068 ) 
                {
                  cat('No fit progress made in the last 10 steps. \nTerminating algorithm.',sep = '\n')
                } 
            
                break 
            }
        }
      
      else if (tau_k <= eta)
        {
          if (printIter == TRUE){
              {if (k%%50 == 0 ){print(cat('Iteration ',k,'\nWarning: No sufficiently descent direction found.\n\n',sep = '\n'))}}}
          
        if (k >1 ) {info[,k] <- info[,k-1]}
        
        if (stoptype == 'step')
          {
            if (info['stepIter',k] <= tol){break}
          }
        else if (stoptype == 'grad')
          {
            if (norm(dphi_k,'i') <= tol*(1 + abs(phi_k))){break}
          }
        else if (stoptype =='gof')
          {
            if (gofstat_k <= 0.068) {break}
          }
        
        if ( k > 100 & gofstat_k > gofstat_k_1 &  gofstat_k <= 0.068 ) {break}
        
        if ( abs(info['sambiasIter',k] - 1) <= 0.05 &  gofstat_k <= 0.068 )
          {
            break
          }
        
        if (k > 100)
          {
            if (round(info['gofIter',k-50],8) - round(gofstat_k,8) == 0 & gofstat_k <= 0.068 ) 
              { 
                cat('No fit progress made in the last 50 steps. \nTerminating algorithm.',sep = '\n')
              } 
              
              break 
          }
        
      }
      
      k <- k + 1
    }
    
    return(list(par = theta_k,
                rhohat = rho(theta_k),
                niter = k-1,
                sambias = info['sambiasIter',k-1],
                gof = info['gofIter',k-1],
                info = data.frame(phis = as.vector(info['phiIter',1:k-1]),
                                   taus = as.vector(info['tauIter',1:k-1]),
                                   stepsizes = as.vector(info['stepIter',1:k-1]), 
                                   gofs = as.vector(info['gofIter',1:k-1]),
                                   sambiases = as.vector(info['sambiasIter',1:k-1])
                                   ) 
                )
          )
  }
  
  fitGN <- function(theta0 = theta0, stoptype = stoptype, tol = tol, maxIter = maxIter,eta = eta, printIter = printIter)
    {
      
      theta_k <- theta0
      
      info <- matrix(0,nrow = 4,ncol = maxIter)
      rownames(info) <- c('gofIter','sambiasIter','stepIter','phiIter')
      
      gofstat_k <- trapzgof(theta_k)
      
      k <- 1
      
      
      while ( k <= maxIter )
        {
          
          J_k <- JAC(theta_k)
          
          rho_k <- rho(theta_k)
          
          if ( any(is.na(J_k))) 
          {
            J_k <- apply( J_k, MARGIN = 2, FUN = na.spline  )
          }
          
          JJ_k <- t(J_k)%*%J_k
          
          dtheta_k <- tryCatch(solve(JJ_k, - t(J_k)%*%rho(theta_k) ), error = function(err) NA)
          
          alpha_k <- linesearch( phi, dphi, theta_k, dtheta_k )
          
          if ( any(is.na(dtheta_k))  ) { break}
          
          if( is.nan(alpha_k) ) { alpha_k <- 1  }
          
          if (theta_k[2] + dtheta_k[2] < 0 ) { dtheta_k[2] <- 0}
          
          if (theta_k[3] + dtheta_k[3] < 0 ) { dtheta_k[3] <- 0}
          
          if ( X[1] - theta_k[1] - dtheta_k[1] <0 ) { dtheta_k[1] <- 0  }
          
          theta_k <- theta_k + alpha_k*dtheta_k
          
          theta_k_1 <- theta_k
          
          gofstat_k_1 <- gofstat_k
          
          gofstat_k <- trapzgof(theta_k)
          
          info['gofIter',k] <- gofstat_k
          info['sambiasIter',k] <- sambias(data,iFx,theta_k,pvar)
          info['stepIter',k] <- l2(theta_k - theta_k_1)/l2(theta_k_1)
          info['phiIter',k] <- phi(theta_k)
          
          
          if (printIter == TRUE)
            {
              if (k%%50 == 0 )
                {
                  cat(paste('Iteration: ',k,
                            '\nPhi_k ',round(info['phiIter',k],5),
                            '\nStep Size: ',round(info['stepIter',k],10),
                            '\nSample 99.9% Bias: ',round(info['sambiasIter',k],5),
                            '\nGoF Stat.: ',round(gofstat_k,5),'\n','\n')
                      )
                }
            }
          
          if (stoptype == 'step')
            {
              if (info['stepIter',k] <= tol){break}
            }
          
          else if (stoptype == 'grad')
            {
              if (norm(dphi(theta_k),'i') <= tol*(1 + abs(phi(theta_k)))){break}
            }
          
          else if (stoptype =='gof')
            {
              if (gofstat_k <= 0.068) {break}
            }
          
          if ( k > 100 & gofstat_k > gofstat_k_1 &  gofstat_k <= 0.068 ) {break}
          
          if ( abs(info['sambiasIter',k] - 1) <= 0.05 &  gofstat_k <= 0.068 )
            {
              break
            }
          
          if (k > 10)
            {
              if (round(info['gofIter',k-10],5) - round(gofstat_k,5) == 0  ) 
                {
                  cat('No fit progress made in the last 10 steps. \nTerminating algorithm.\n',sep = '\n')
                } 
            
                break 
            }
          
          k <- k +1
          
        }
      
      return(list(par = theta_k,
                  rhohat = rho(theta_k),
                  niter = k-1,
                  sambias = info['sambiasIter',k-1],
                  gof = gofstat_k,
                  info = data.frame(
                                    phis = as.vector(info['phiIter',1:k-1]),
                                    stepsizes = as.vector(info['stepIter',1:k-1]), 
                                    gofs = as.vector(info['gofIter',1:k-1]),
                                    sambiases = as.vector(info['sambiasIter',1:k-1])
                                    ) 
                  )
            )
    }
  
  if (algo == 'LM')
    {
      results <- fitLM(theta0 = theta0,rad0 = rad0, stoptype = stoptype, tol = tol, maxIter = maxIter,eta = eta, printIter = printIter)
    }
  
  else if (algo == 'GN')
    {
      results <- fitGN(theta0 = theta0, stoptype = stoptype, tol = tol, maxIter = maxIter,eta = eta, printIter = printIter)
    }
  
  if (results$gof <= 0.068) {  cat('\nOptimisation successful.\n',sep = '\n') }
  else { cat('\nCould not find optimal parameter values. Returning best solution found.\n',sep = '\n') }
  
  if( printIter == TRUE ) {  cat('\nNow running Monte Carlo VaR',sep = '\n') }
  
  if ( run_VaR == TRUE)
    {
     mcvar <- fnvar(data,pgpd,rgpd,period,threshold,100000,as.vector(results$par),pvar) 
    }
  
  else
    {
      mcvar <- NULL
    }
  
  xi_est <- xi_range(data,0.1)
  xi_h <- xi_est$xi_hill
  xi_h_conf <- xi_est$xi_hill_conf
  t_min_h <- match(min(xi_h),xi_h)
  t_max_h <- match(max(xi_h),xi_h)
  xi_h_low_min <- xi_h[t_min_h] - xi_h_conf[t_min_h]
  xi_h_low_max <- xi_h[t_min_h] - xi_h_conf[t_min_h]
  xi_h_high_min <- xi_h[t_max_h] - xi_h_conf[t_max_h]
  xi_h_high_max <- xi_h[t_max_h] - xi_h_conf[t_max_h]
  
  xi_p <- xi_est$xi_pickands
  xi_p_conf <- xi_est$xi_pickands_conf
  t_min_p <- match(min(xi_p),xi_p)
  t_max_p <- match(max(xi_p),xi_p)
  xi_p_low_min <- xi_p[t_min_p] - xi_p_conf[t_min_p]
  xi_p_low_max <- xi_p[t_min_p] + xi_p_conf[t_min_p]
  xi_p_high_min <- xi_p[t_max_p] - xi_p_conf[t_max_p]
  xi_p_high_max <- xi_p[t_max_p] + xi_p_conf[t_max_p]
  
  if (printIter == TRUE)
  {
    cat('\nStage 1 ', algo,
        '\nNo. of Iterations: ',results$niter,
        '\n\nGoF. Stat. : ', results$gof, 
        '\nTail (95%) Bias: ', iFx(1 - 0.05/(1 - cut_pct),results$par)/quantile(data,1 - 0.05/(1 - cut_pct)),
        '\nTail (99.9%) Bias: ',results$sambias,
        '\nMu: ',results$par[1],
        '\nSigma: ',results$par[2],
        '\nXi: ',results$par[3],
        '\nImplied Tail Threshold: ', full_ecdf(results$par[1]),
        '\nMedian Xi Hill: ', median(xi_h),
        '\nMinimum Xi Hill Function 90% Range with asymptotic normal 95% CI: ',xi_h_low_min,xi_h_low_max,
        '\nMaximum Xi Hill Function 90% Range with asymptotic normal 95% CI: ',xi_h_high_min, xi_h_high_max,
        '\nMedian Xi Pickands: ', median(xi_p),
        '\nMinimum Xi Pickands Function 90% Range with asymptotic normal 95% CI: ',xi_p_low_min , xi_p_low_max,
        '\nMaximum Xi Pickands Function 90% Range with asymptotic normal 95% CI: ',xi_p_high_min ,xi_p_high_max,
        '\n100,000 trials of MC VaR in mlns GBP at 99.9%: ',mcvar,sep = '\n')
  }
  
  F_data_frame <- data.frame( X = X,F_n = Fn, F_x = Fx(results$par) )
  
  return(list(par = results$par,
              trapzgof = trapzgof,
              F_data_frame = F_data_frame,
              fullresults = results,
              mcvar = mcvar,
              xi_h = xi_h, 
              xi_h_conf = xi_h_conf,
              xi_p = xi_p,
              xi_p_conf = xi_p_conf)
         )

}








################################ Stage 2 ODR (using Levenberg-Marquadt) ################################

stage2 <- function(data,theta0, v = 1, gamma = NULL,rhohat = NULL,epsilon0 = NULL,
                   period = 6.5,cut_pct = 0, w2 = NULL, fittype = 'A',
                   mquants = trunc(length(data)/4),rad0 = 500, stoptype = 'grad',
                   tol = 1e-10, maxIter = 1000,eta = 0.2,pvar = 0.999,printIter = FALSE, 
                   run_VaR = TRUE)
{
  # data:         Input a VECTOR/LIST of data.
  # theta0:       Initial guess for theta = [mu,sigma,xi] in that order.
  # v,gamma:      If 'log' then returns a weight vector whose last element is v times as large as
  #               the first element with logarithmically increasing values. If 'invres' returns
  #               a weight vector whos elemtns are the reciprocals of the predicted residuals' from
  #               stage 1 rhohat scaled by gamma.
  # rhohat:       Vector of stage 1 residuals. Needed for invres weights.
  # epsilon0:     Initial guess for epsilon, the vector of estimated perturbations to the 
  #               observation coordinates.
  # period:       Length of the time span of the dataset.
  # cut_pct:      If the dataset is a truncated version of another, this is the percentage that 
  #               has been cut from the original.
  # w2:           Weight vector for the perturbations to the observation coordinates. If there are 
  #               no suspected bad data points present the recommended value is w2 = ones(length(data)).
  # fittype:      Type A fits m quantiles with m + 1 equidistant intervals in [min(data),max(data)] 
  #               (this yields more tail quantiles). Type B fits m quantiles using their correspondong 
  #               m+1 equidistant intervals in F(X_i) in [1/m,1] for i = 1,...,m (this yields more 
  #               body quantiles but with less interpolation).
  # mquants:       Number of quantiles to take (more quantiles yield more accurate results but takes 
  #               longer, my experience indicates a quarter to a third of the length of the data set 
  #               is sufficient.
  # rad0:         Initial maximum search radius. Just input a reasonable number (say 5 < rad0 < 5000) 
  #               because it is going to get altered by the maximum search algorithm.
  # stoptype:     Convergence criterion. If one just wants a passing fit then 'gof' should do the 
  #               trick. Otherwise standard practice is by relative normalised step size or infinity 
  #               norm of step size (see paper for more detail).
  # tol:          Convergence tolerance. Should be quite small if the residuals are small 
  #               (i.e. the tail is GPD)
  # maxIter:      Maximum number of iterations. Note that one should be flexible with this if they 
  #               wanted a good fit for a difficult problem.
  # eta:          Threshold for accuracy of the quadratic model. This is the ratio between actual 
  #               decrease and the predicted decrease through the JACobian approximation (see 
  #               paper for more detail).
  # pvar:         Threshold percentile to calculate the MC VaR
  # prinIter:     Set to TRUE for details of every 50th iteration.
  # run_VaR:      If TRUE then algorithm automatically runs an MC VaR using the fitted parameter values
  #               with 100,000 trials.
  
  frequency <- length(data)/period
  
  # If the data is truncated
  threshold <- unname(quantile(data,cut_pct))
  freauency_normed <- frequency/(1 - cut_pct)
  pvar <- 1 - 0.001/(1 - cut_pct)
  full_ecdf <- ecdf(data)
  data <- data[ data >= threshold ]
  
  # Checking inputs
  if ( is.null(mquants) ) { mquants <- trunc(length(data)/4) }

  else if (mquants > length(data)) { stop('Number of quantiles cannot be greater than length of the data set.') }
  
  if (is.null(epsilon0)) {epsilon0 <- rep(0,mquants)}
  
  # Weight functions
  ones_m <- rep(1,mquants)
  
  if ( is.null(gamma) | is.null(rhohat) ) 
    {
      w1 <- getweights(mquants, v = v , gamma = NULL, rhohat = NULL)
    }
  
  else
    {
      w1 <- getweights(mquants, v = NULL , gamma = gamma, rhohat = rhohat)
    }
  
  if ( !is.null(gamma) & is.null(rhohat))
    {
      stop('A residual vector must be provided for invres type of weights.')
    }
  
  if (is.null(w2)) {w2 <- ones_m}
  
  else if (w2 == 0) {w2 <- rep(0,mquants)}

  dataset <- processdata(data,fittype ,mquants)
  X <- dataset$X
  Fn <- dataset$Fn
  ecdf <- dataset$ecdf
  frequency <- length(data)/period
  
  # Checking dimensions
  if ( length(epsilon0) != mquants ) { 
    stop('Epsilon vector must be the same length as number of quantiles.') } 
  
  if ( length(w1) != mquants ) { 
    stop('Weight_1 vectors must be the same length as number of quantiles.') } 
  
  if ( length(w2) != mquants ) { 
    stop('Weight_2 vectors must be the same length as number of quantiles.') } 
  
  
  
  ########### Optimisation functions #############
  
  #Fx <- function(theta,epsilon) { pgpd( X + epsilon, loc = theta[1], scale = theta[2], shape = theta[3] )  }
  
  Fx <- function(theta,epsilon) {1 - (1 + (theta[3]/theta[2])*(X + epsilon - theta[1]))^(-1/theta[3]) }
  
  # There is no need to include epsilon in this icdf because evaluations of the bias are done after the
  # shifting of the X coordinates with epsilon.
  iFx <- function(q,theta){ return((theta[2]/theta[3])*((1 - q )^(-theta[3]) - 1) + theta[1] )}

  trapzgof <- function(theta,epsilon){return(trapz(Fn,abs(Fn - Fx(theta,epsilon))))}

  # JAC is a matrix.
  JAC_theta <- function(theta,epsilon)
  {   cbind(
      (1/theta[2])*w1*( 1 + (theta[3]/theta[2])*(X + epsilon - theta[1]))^(-1/theta[3]-1)  ,
      
      ((X + epsilon - theta[1])/theta[2]^2)*w1*( 1 + (theta[3]/theta[2])*(X + epsilon - theta[1]))^(-1/theta[3]-1) ,
      
      w1*((1/theta[3]^2)*log( 1 + (theta[3]/theta[2])*(X - theta[1])) - (X - theta[1])/(theta[2]*theta[3] +
      (theta[3]^2)*(X - theta[1])) )*( 1 + (theta[3]/theta[2])*(X - theta[1]))^(-1/theta[3])        )
  }

  # jac is a vector. 
  jac_epsilon <- function(theta,epsilon)
  {-(1/theta[2])*w1*( 1 + (theta[3]/theta[2])*(X + epsilon - theta[1]))^(-1/theta[3]-1)}
  

  rho1 <- function(theta,epsilon){ return( w1*(Fn - Fx(theta,epsilon)) ) }

  rho2 <- function(epsilon){ return( w2*epsilon ) }

  rho <- function(theta,epsilon){ return( c(rho1(theta,epsilon),rho2(epsilon)) ) }

  phi <- function(theta,epsilon){ return( 0.5*l2(rho(theta,epsilon))^2 ) }

  lsqrdir <- function(theta_k,epsilon_k,rad0,J_t,j_e,maxIter = 100)
  {
    # I designed the next few bizarre functions to minimise the order of complexity. For more details
    # I direct the user to the ODR Efficient implementation section in the appendix of the paper.
    lambda <- 0.1
    rho1_k <- rho1(theta_k,epsilon_k)
    rho2_k <- rho2(epsilon_k)

    einv <- function( lambda) {as.vector( 1/( j_e^2 + w2^2 + lambda) )}
    h <- function( lambda ) { as.vector(1 - einv(lambda)*( j_e^2 ) )}
    aux <- as.vector(j_e*( j_e*rho1_k + (w2^2)*rho2_k ))
    rhostar <- function( lambda ) { rho1_k - einv(lambda)*aux }
    
    J_k <- J_t*sqrt(h(0))

    dtheta <- rep(0,length(theta_k))
  
    QR0 <- qr(J_k)
    Q0 <-qr.Q(QR0)
    R0 <-qr.R(QR0)
    dimQR <- dim(J_k); n <- dimQR[2]

    k <- 0
    
    while ( k <= maxIter & abs(l2(dtheta) - rad0) > 1e-12 & lambda > 0 )
    {
      Q_k <- Q0
      R_k <- R0
    

      for (i in 1:n){ R_k <- qrinsert(R_k,sqrt(lambda),i) }

      aux2 <- -t(J_t)%*%rhostar(lambda)
      dtheta <- backsolve(R_k, backsolve(R_k,aux2,transpose = TRUE))
      eigvec <- backsolve(R_k,dtheta,transpose = TRUE)
      
      l2_eigvec_k <- l2(eigvec)
      
      if (l2_eigvec_k == 0) { l2_eigvec_k <- eps(1) }

      lambda <- max(0, lambda + ((l2(dtheta)/l2_eigvec_k)^2)*(l2(dtheta) - rad0)/rad0)
      
      k <- k + 1
    }
    
    if (lambda == 0)
    {
      dtheta <- backsolve(R_k, backsolve(R_k,aux2,transpose = TRUE))
    }
    
    depsilon <- -einv(lambda)*( j_e*rho1_k + w2*rho2_k + j_e*( J_t%*%dtheta ))
    
    return(list(dtheta = as.double(dtheta),depsilon = depsilon))
  }
  
  getradM <- function(theta0,epsilon0,rad0,maxIter = 100)
  {
    Jt <- JAC_theta(theta0,epsilon0)
    je <- jac_epsilon(theta0,epsilon0)
    
    if ( any(is.na(je))) 
      {
        je <- na.spline(je)
      }
      
      if ( any(is.na(Jt))) 
      {
        Jt <- apply( Jt, MARGIN = 2, FUN = na.spline  )
      }
    
    dtheta_0 <- lsqrdir(theta0,epsilon0,rad0,Jt,je,maxIter)$dtheta
    
    radM <- rad0
    
    if (l2(dtheta_0) <= 0.01*rad0){ radM <- 5*l2(dtheta_0)}
    
    while (abs(l2(dtheta_0)) <= 1e-5)
    {
      radM <- 5*radM
      dtheta_0 <- lsqrdir(theta0,epsilon0,radM,Jt,je,maxIter)
    }
    
    return (radM)
  }
  
  # This function evaluates the quadratic model in O(super-linear) because of the 0's in the sub-JACobians
  quadmodel <- function(dphi_k,J_t,j_e,rho1_k,rho2_k,dtheta,depsilon)
  {
    rho_k <- c(rho1_k,rho2_k)
    delta <- c(dtheta,depsilon)
    temp1 <- as.vector(J_t%*%dtheta)
    temp2 <-2*dot(temp1, j_e*depsilon) + dot(as.vector(j_e^2 + w2^2),depsilon)
    
    q <- 0.5*l2(rho_k)^2 + dot(delta, dphi_k) + l2(temp1)^2 + temp2
    
    return(q)
  }
  
  fitODR <- function(stoptype = stoptype, tol = tol, maxIter = maxIter,
                     eta = eta, printIter = printIter)
  {
    theta_k <- theta0
    epsilon_k <- epsilon0
    
    radM <- getradM(theta0,epsilon0,rad0)
    rad_k <- radM
    
    info <- matrix(0,nrow = 5,ncol = maxIter)
    rownames(info) <- c('gofIter','sambiasIter','tauIter','stepIter','phiIter')
    
    gofstat_k <- trapzgof(theta_k,epsilon_k)
    
    gofstat_k_1 <- gofstat_k
    
    k <- 1
    
    while (k <= maxIter)
    {
      rho1_k <- rho1(theta_k,epsilon_k)
      rho2_k <- rho2(epsilon_k)
      
      Jt_k <- JAC_theta(theta_k,epsilon_k)
    
      je_k <- jac_epsilon(theta_k,epsilon_k)
      
      if ( any(is.na(je_k))) 
        {
          je_k <- na.spline(je_k)
        }
      
      if ( any(is.na(Jt_k))) 
        {
          Jt_k <- apply( Jt_k, MARGIN = 2, FUN = na.spline  )
        }
      
      dphi_k <- c(  t(Jt_k)%*%rho1_k, je_k*rho1_k + w2*rho2_k )
      
      q_k <- function(dthe,deps){ quadmodel( dphi_k, Jt_k, je_k, rho1_k, rho2_k, dthe, deps ) }
      
      direction <- lsqrdir(theta_k,epsilon_k,rad_k,Jt_k,je_k)
      
      dtheta_k <- direction$dtheta
      
      depsilon_k <- direction$depsilon
      
      phi_k <- phi(theta_k,epsilon_k)
    
      
      if (theta_k[2] + dtheta_k[2] < 0){dtheta_k[2] <- 0}
      if (theta_k[2] <= 0 ){ theta_k[2] <- theta_k_1[2] }
      
      if (theta_k[3] + dtheta_k[3] < 0){dtheta_k[3] <- 0}
      
      if ( min(X + epsilon_k + depsilon_k) - (theta_k[1] + dtheta_k[1]) <0 ) { dtheta_k[1] <- 0  }
      
      tau_k <- ( phi(theta_k,epsilon_k) + phi(theta_k + dtheta_k, epsilon_k + depsilon_k) )/( 
        q_k(rep(0,length(theta_k)),rep(0,length(epsilon_k))) - q_k(dtheta_k,depsilon_k) )
      
      if (is.nan(tau_k)) { tau_k <- 0 }
      
      if (tau_k < 0.25){ rad_k <- 0.25*rad_k }
      
      else if ( (tau_k > 0.75) & (abs(l2(dtheta_k)^2 - rad_k^2) < 1e-12) ){ rad_k <- min(2*rad_k,radM) }
      
      theta_k_1 <- theta_k
      
      if (tau_k > eta)
        {
          theta_k <- theta_k + dtheta_k
          epsilon_k <- epsilon_k + depsilon_k
          
          gofstat_k_1 <- gofstat_k
          gofstat_k <- trapzgof(theta_k,epsilon_k)
          
          X <- X + epsilon_k
          
          info['gofIter',k] <- gofstat_k
          info['sambiasIter',k] <- sambias(data,iFx,theta_k,pvar)
          info['tauIter',k] <- tau_k
          info['stepIter',k] <- l2(theta_k - theta_k_1)/l2(theta_k_1)
          info['phiIter',k] <- phi(theta_k,epsilon_k)
          
          if (printIter == TRUE)
            {
              if (k%%50 == 0 )
                {
                  cat(paste('Iteration: ',k,
                            '\nPhi_k: ',round(info['phiIter',k],5),
                            '\nStep Size: ',round(info['stepIter',k],10),
                            '\nSample 99.9% Bias: ',round(info['sambiasIter',k],5),
                            '\nGoF Stat.: ',round(gofstat_k,5),'\n','\n'))
                }
            }
        
          # Stopping conditions
          if (stoptype == 'step')
            {
              if (info['stepIter',k] <= tol) {break}
            }
          
          else if (stoptype == 'grad')
            {
              if (max(abs(dphi_k)) <= tol*(1 + abs(phi_k))) {break}
            }
          
          else if (stoptype =='gof')
            {
              if (gofstat_k <= 0.068) {break}
            }
          
          if (  k > 100 & gofstat_k > gofstat_k_1 &  gofstat_k <= 0.068 )
            {
              break
            }
          
          if ( abs(info['sambiasIter',k] - 1) <= 0.05 &  gofstat_k <= 0.068 )
          {
            break
          }
          
          if (k > 50 && round(info['gofIter',k-50],10) - round(gofstat_k,10) == 0  )
            {
              if ( printIter == TRUE)
                {
                  cat('No fit progress made in the last 50 steps. \nTerminating algorithm.',sep = '\n')
                }
  
              break
            }
        
        }
      
      else if (tau_k <= eta)
      {
        gofstat_k_1 <- gofstat_k
        gofstat_k <- trapzgof(theta_k,epsilon_k)
        
        
        if (printIter == TRUE)
            {
              {if (k%%50 == 0 )
                {
                print(cat('Iteration ',k,'\nNo sufficiently descent direction found.',sep = '\n'))
                
                break
                }
              }
            }
        
        if (k >1 ) 
          {
            info[,k] <- info[,k-1]
          }
        
        if (stoptype == 'step')
          {
            if (info['stepIter',k] <= tol) {break}
          }
        
        else if (stoptype == 'grad')
          {
            if (max(abs(dphi_k)) <= tol*(1 + abs(phi_k))) {break}
          }
        
        else if (stoptype =='gof')
          {
            if (gofstat_k <= 0.068) {break}
          }
        
        if ( k > 100 && gofstat_k > gofstat_k_1 &&  gofstat_k <= 0.068 )  
          {
            break
          }
        
        if ( abs(info['sambiasIter',k] - 1) <= 0.05 &  gofstat_k <= 0.068 )
          {
            break
          }
        
        if (k > 50 && round(info['gofIter',k-50],10) == round(gofstat_k,10)   )
          {
            if ( printIter == TRUE)
              {
                cat('No fit progress made in the last 50 steps. \nTerminating algorithm.',sep = '\n')
              }
            
              break
          }
      }
      
      k <- k + 1
      
    }
    
    return(list(par = theta_k,
                epsilon = epsilon_k,
                rho1hat = rho1(theta_k,epsilon_k),
                rho2hat = rho2(epsilon_k),
                niter = k-1,
                sambias = info['sambiasIter',k],
                gof = gofstat_k,
                info = data.frame(phis = as.vector(info['phiIter',1:k]),
                                  taus = as.vector(info['tauIter',1:k]),
                                  stepsizes = as.vector(info['stepIter',1:k]), 
                                  gofs = as.vector(info['gofIter',1:k]),
                                  sambiases = as.vector(info['sambiasIter',1:k])
                                  )
                ) 
           )
    
  }
  
  odrresults <- fitODR(stoptype = stoptype, tol = tol, maxIter = maxIter,
                 eta = eta, printIter = printIter)
  
  if (odrresults$gof <= 0.068) {  cat('\nOptimisation successful.\n',sep = '\n') }
  else { cat('\nCould not find optimal parameter values. Returning best solution found.\n',sep = '\n') }
  
  if( printIter == TRUE ) {  cat('\nNow running Monte Carlo VaR',sep = '\n') }
  
  if ( run_VaR == TRUE)
    {
     mcvar <- fnvar(data,pgpd,rgpd,period,threshold,100000,as.vector(odrresults$par),pvar) 
    }
  
  else
    {
      mcvar <- NULL
    }
  
  xi_est <- xi_range(data,0.1)
  xi_h <- xi_est$xi_hill
  xi_h_conf <- xi_est$xi_hill_conf
  t_min_h <- match(min(xi_h),xi_h)
  t_max_h <- match(max(xi_h),xi_h)
  xi_h_low_min <- xi_h[t_min_h] - xi_h_conf[t_min_h]
  xi_h_low_max <- xi_h[t_min_h] - xi_h_conf[t_min_h]
  xi_h_high_min <- xi_h[t_max_h] - xi_h_conf[t_max_h]
  xi_h_high_max <- xi_h[t_max_h] - xi_h_conf[t_max_h]
  
  xi_p <- xi_est$xi_pickands
  xi_p_conf <- xi_est$xi_pickands_conf
  t_min_p <- match(min(xi_p),xi_p)
  t_max_p <- match(max(xi_p),xi_p)
  xi_p_low_min <- xi_p[t_min_p] - xi_p_conf[t_min_p]
  xi_p_low_max <- xi_p[t_min_p] + xi_p_conf[t_min_p]
  xi_p_high_min <- xi_p[t_max_p] - xi_p_conf[t_max_p]
  xi_p_high_max <- xi_p[t_max_p] + xi_p_conf[t_max_p]
  
  if (printIter == TRUE)
    {
      cat('\nODR Stage 2.\nNo. of Iterations: ',odrresults$niter,
          '\nNo. of Quantiles: ', mquants,
          '\n\nGoF. Stat. : ', odrresults$gof, 
          '\nTail (95%) Bias: ', iFx(1 - 0.05/(1 - cut_pct),odrresults$par)/quantile(data,1 - 0.05/(1 - cut_pct)),
          '\nTail (99.9%) Bias: ',odrresults$sambias,
          '\nMu: ',odrresults$par[1],
          '\nSigma: ',odrresults$par[2],
          '\nXi: ',odrresults$par[3],
          '\nImplied Tail Threshold: ', full_ecdf(odrresults$par[1]),
          '\nMedian Xi Hill: ', median(xi_h),
          '\nMinimum Xi Hill Function 90% Range with asymptotic normal 95% CI: ',xi_h_low_min,xi_h_low_max,
          '\nMaximum Xi Hill Function 90% Range with asymptotic normal 95% CI: ',xi_h_high_min, xi_h_high_max,
          '\nMedian Xi Pickands: ', median(xi_p),
          '\nMinimum Xi Pickands Function 90% Range with asymptotic normal 95% CI: ',xi_p_low_min , xi_p_low_max,
          '\nMaximum Xi Pickands Function 90% Range with asymptotic normal 95% CI: ',xi_p_high_min ,xi_p_high_max,
          '\n100,000 trials of MC VaR in mlns GBP at 99.9%: ',mcvar,sep = '\n'
          )
    }
  
  F_data_frame <- data.frame( X = X, F_n = Fn, F_x = Fx(odrresults$par,odrresults$epsilon) )
  
  return(list(par = odrresults$par,
              F_data_frame = F_data_frame,
              fullresults = odrresults,
              X = X, 
              Fn = Fn, 
              F_pred = Fx(odrresults$par,odrresults$epsilon),
              trapzgof = trapzgof,
              mcvar = mcvar,
              xi_h = xi_h, 
              xi_h_conf = xi_h_conf,
              xi_p = xi_p,
              xi_p_conf = xi_p_conf))
  
}

getweights <- function(mquants, v = NULL , gamma = NULL, rhohat = NULL)
{
  # Log weight type increases logarithmically from 0 to param (v in paper) for according to the increasing quantiles
  # Inverse residual weight type weighs the reciprocals of the first stage residuals by param (gamma in paper)
  if (is.null(v))
  {
    if (is.null(rhohat)){ print('Error: A vector of residuals must be provided as rhohat for invres.')}
    
    if ( is.null(gamma) ){ return( (500/max((rhohat))^2)*( 1/((rhohat)^2) ) )  }
    
    else {return( gamma*( 1/((rhohat)^2 ) )) }
  }
  
  else 
  {
    return( exp(log(10)*seq(log10(1),log10(v),length.out =  mquants)) ) 
  }
  
}





################################ Tail Index Estimators ########################

xi_hill <- function(data)
{
  # This coincides with MLE for 1-param Pareto with xi = 1/alpha. The original paper Hill (1975) derives this 
  # as the gradient of the tail value frequency per the power law. for every cutoff point k.
  
  n <- length(data)
  N <- seq.int(1,n,1)
  revdata <- sort(data, decreasing = TRUE)
  logrevdata <- log(revdata)
  
  xi <- (1/N)*cumsum(logrevdata) - logrevdata
  
  xi_conf <- (1.96/sqrt(N))*xi
  
  return(data.frame(k = N,xi = xi,conf = xi_conf))
}

xi_pickands <- function(data)
{
  # Takes subsamples of tails of increasing sizes using the cutoff, the 50-pct and 75-pct of the truncated
  # datasets. See Pickands (1975) for more details.
  
  n <- length(data)
  N <- seq.int(1,n,1)
  revdata <- sort(data, decreasing = TRUE)
  
  x75 <- revdata[ seq(2,trunc(length(revdata)/4),1) ]
  x50 <- revdata[ seq(3,2*trunc(length(revdata)/4),2) ]
  x0 <- revdata[ seq(5,4*trunc(length(revdata)/4),4) ]
  
  xi <- (1/log(2))*log( (x75 - x50)/(x50 - x0 ) )
  
  sigma <- (x50 )/( (2^xi - 1)/xi)
  
  xi_conf <- (1.96/N[1:length(xi)])*sqrt( ((xi^2)*(2^(xi +1)))/(  (2*(2^xi -1)*log(2) ))^2 )
  
  gof <- rep(0,length(xi))
  
  for ( i in 1:length(xi) )
    {
      cut_data <- revdata[1:length(data[data >= x0[i]])]
      
      gof[i] <- tryCatch(trapzgof(cut_data, c( x0[i], sigma[i], xi[i] ) ), error = function(err) NA )
      
    }
  
  return(data.frame(k = 4*N[1:length(xi)],xi = xi, sigma = sigma, conf = xi_conf, gof = gof))
}

xi_range <- function(data,p)
{
  # This function selects the (100 x p)% range of the xi's produced by the Hill and Pickands functions.
  # This is done to remove the first and last few observations that produce wildly large swings.
  
  p <- p/2
  q = 1 - p
  
  xi_h_results <- xi_hill(data)
  xi_h <- xi_h_results$xi
  xi_h_conf <- xi_h_results$conf
  xi_h_ord <- sort(xi_h)
  xi_h_trunc <- xi_h[(xi_h >= xi_h_ord[trunc(p*length(data))]) & (xi_h <= xi_h_ord[trunc(q*length(data))]) ]
  
  xi_p_results <- xi_pickands(data)
  xi_p <- xi_pickands(data)$xi
  xi_p_conf <- xi_p_results$conf
  xi_p_ord <- sort(xi_p)
  xi_p_trunc <- xi_p[(xi_p >= xi_p_ord[trunc(p*length(xi_p))]) & (xi_p <= xi_p_ord[trunc(q*length(xi_p))]) ]
  
  return(list( xi_hill =  xi_h_trunc ,
               xi_hill_conf =  xi_h_conf,
               xi_pickands = xi_p_trunc,
               xi_pickands_conf =  xi_p_conf
               )
         )
  
}


######################### Testing stability of the stage 2 WQE tail ##################################

theo_tail_cut_param <- function(param,t)
{
  # This returns the theoretical parameter values of the truncated tail given an initial set of tail 
  # estimators based on the max-stable characteristic of the GPD (tail index is the same if GPD applies).
  if ( length(t) > 1)
    {
      return( cbind(t, param[2] + param[3]*(t - param[1]) , param[3]) )
    }
  
  else
    {
      return( c(t, param[2] + param[3]*(t - param[1]) , param[3]) )
    }
  
}


get_nIter_vector <- function (maxweight,cutoffs)
{
  # This is only a convenient approximate no. of iterations needed for the given weights. One should 
  # always check if the iterations are enough for convergence.
  n_cuts <- length(cutoffs)
  
  if ( maxweight <= 40  ) { return( seq.int(3000, 5000, length.out = n_cuts )) }
  
  else if ( maxweight > 40 & maxweight <= 60 )
  { return( seq.int(3000, 7000, length.out = n_cuts )) }
  
  else if ( maxweight > 60 & maxweight <= 80 )
  { return( seq.int(4000, 10000, length.out = n_cuts ))  }
  
  else if ( maxweight > 80 & maxweight<= 100 )
  { return( seq.int(5000, 13000, length.out = n_cuts )) }
  
  else if ( maxweight > 100 & maxweight <= 130 )
  { return( seq.int(8000, 16000, length.out = n_cuts )) }
  
  else if ( maxweight > 130 & maxweight <= 200 )
  { return( seq.int(10000, 30000, length.out = n_cuts )) }
  
  else if ( maxweight > 200 & maxweight <= 500 )
  { return( seq.int(10000, 40000, length.out = n_cuts )) }
  
  else { return( seq.int(10000, 30000, length.out = n_cuts )) }
}


truncate_s1_tail <- function(data,theta0 = NULL,tail_cutoffs = NULL,
                             mquants = NULL,fittype = 'A', period = 6.5,rad0 = 500,
                             stoptype = 'step',tol = 1e-10, eta = 0.2,
                             nIter_vector = NULL, print_cut_Iter = TRUE)
{
  if ( is.null(mquants)) {  mquants <- trunc(length(data)/4)}
  
  init_mquants <- mquants
  
  if (is.null(tail_cutoffs)) { 
    tail_cutoffs <- seq(0,0.8,by=0.2) 
    
    nIter_vector <- seq(10000,20000,length.out = length(tail_cutoffs)) 
    }
  
  
  thresholds <- unname(quantile(data,tail_cutoffs))
  
  init_results_GN <- stage1(data,theta0 = theta0,algo = 'GN', period = period,cut_pct = tail_cutoffs[1],
                         fittype = fittype, mquants = init_mquants, stoptype = stoptype,
                         tol = tol, maxIter = nIter_vector[1], eta = eta,pvar = 1 - 0.001/( 1 - tail_cutoff[1]), printIter = TRUE)
  
  theo_params_GN <- as.matrix(t(apply( as.matrix(thresholds),MARGIN = 1,
                                    function(t) theo_tail_cut_param(init_results_GN$par,t))))
  
  results_GN <- matrix(0,nrow = length(tail_cutoffs),ncol = 13)
  
  colnames(results_GN) <-  c('cuts',
                          'mu_f','sigma_f','xi_f','mcvar_f','gof_f','sambias_f',
                          'mu_t','sigma_t','xi_t','mcvar_t','gof_t','sambias_t')
  
  results_GN[,'cuts'] <- tail_cutoffs
  results_GN[1,c('gof_t','gof_f')] <- init_results_GN$fullresults$gof
  results_GN[1,c('sambias_t','sambias_f')] <- init_results_GN$fullresults$sambias
  results_GN[,'sambias_t'] <- init_results_GN$fullresults$sambias
  results_GN[1,c('mu_f','sigma_f','xi_f')] <- theo_params_GN[1,]
  results_GN[,c('mu_t','sigma_t','xi_t')] <- theo_params_GN[,]
  results_GN[1,c('mcvar_f','mcvar_t')] <- init_results_GN$mcvar
  
  theta_c_GN <- theo_params_GN[1,]
  
  init_results_LM <- stage1(data,theta0 = theta0,algo = 'LM', period = period,cut_pct = tail_cutoffs[1], rad0 = rad0,
                            fittype = fittype, mquants = init_mquants, stoptype = stoptype,
                            tol = tol, maxIter = nIter_vector[1], eta = eta,pvar = 1 - 0.001/( 1 - tail_cutoff[1]), printIter = TRUE)
  
  theo_params_LM <- as.matrix(t(apply( as.matrix(thresholds),MARGIN = 1,
                                       function(t) theo_tail_cut_param(init_results_LM$par,t))))
  
  results_LM <- matrix(0,nrow = length(tail_cutoffs),ncol = 13)
  
  colnames(results_LM) <-  c('cuts',
                             'mu_f','sigma_f','xi_f','mcvar_f','gof_f','sambias_f',
                             'mu_t','sigma_t','xi_t','mcvar_t','gof_t','sambias_t')
  
  results_LM[,'cuts'] <- tail_cutoffs
  results_LM[1,c('gof_t','gof_f')] <- init_results_LM$fullresults$gof
  results_LM[1,c('sambias_t','sambias_f')] <- init_results_LM$fullresults$sambias
  results_LM[,'sambias_t'] <- init_results_LM$fullresults$sambias
  results_LM[1,c('mu_f','sigma_f','xi_f')] <- theo_params_LM[1,]
  results_LM[,c('mu_t','sigma_t','xi_t')] <- theo_params_LM[,]
  results_LM[1,c('mcvar_f','mcvar_t')] <- init_results_LM$mcvar
  
  theta_c_LM <- theo_params_LM[1,]
  
  for ( i in 2:length(tail_cutoffs) )
  {
    cut <- tail_cutoffs[i]
    cut_threshold <- thresholds[i]
    
    theo_param_GN <- theo_params_GN[i,]
    theo_param_LM <- theo_params_LM[i,]
    
    cut_data <- data[data > cut_threshold]
    
    cut_mquants <- trunc((length(data[data > cut_threshold])/length(data))*init_mquants)
    
    tailresults_GN <- stage1(cut_data,theta0 = theta0, algo = 'GN',period = period,
                             cut_pct = 0,fittype = fittype, mquants = cut_mquants,stoptype = stoptype,
                             tol = tol, maxIter = nIter_vector[i],pvar = 1 - 0.001/(1-cut),printIter = FALSE)
    
    results_GN[i,c('mu_f','sigma_f','xi_f')] <- tailresults_GN$par
    
    results_GN[i,'gof_f'] <- tailresults_GN$fullresults$gof
    results_GN[i,'gof_t'] <- tailresults_GN$trapzgof( theo_param_GN)
    
    results_GN[i,'mcvar_f'] <- tailresults_GN$mcvar
    results_GN[i,'mcvar_t'] <- fnvar(data,pgpd,rgpd,period,cut_threshold,100000,theo_param_GN,1 - 0.001/(1 - cut))
    
    results_GN[i,'sambias_f'] <- tailresults_GN$fullresults$sambias
    
    if ( print_cut_Iter == TRUE)
    {
      cat('\n\nGauss-Newton\nCutting the first ',cut*100,'% of the datas set off from ',cut_threshold,
          '\nFitted 99.9% Sample Bias: ',results_GN[i,'sambias_f'],
          'Fitted 99.9% MC Var in mlns GBP: ',results_GN[i,'mcvar_f'],
          '\nPredicted 99.9% Sample Bias: ',results_GN[i,'sambias_t'],
          'Predicted 99.9% MC Var in mlns GBP: ',results_GN[i,'mcvar_t'],
          '\n',sep = '\n')
    }
    
    tailresults_LM <- stage1(cut_data,theta0 = theta0, algo = 'LM',period = period,
                             cut_pct = 0,fittype = fittype, mquants = cut_mquants,rad0 = rad0,stoptype = stoptype,
                             tol = tol, maxIter = nIter_vector[i],eta = eta, pvar = 1 - 0.001/(1-cut),printIter = FALSE)
    
    results_LM[i,c('mu_f','sigma_f','xi_f')] <- tailresults_LM$par
    
    results_LM[i,'gof_f'] <- tailresults_LM$fullresults$gof
    results_LM[i,'gof_t'] <- tailresults_LM$trapzgof( theo_param_LM)
    
    results_LM[i,'mcvar_f'] <- tailresults_LM$mcvar
    results_LM[i,'mcvar_t'] <- fnvar(data,pgpd,rgpd,period,cut_threshold,100000,theo_param_LM,1 - 0.001/(1 - cut))
    
    results_LM[i,'sambias_f'] <- tailresults_LM$fullresults$sambias
    
    if ( print_cut_Iter == TRUE)
    {
      cat('\n\n Levenberg-Marquadt\nCutting the first ',cut*100,'% of the datas set off from ',cut_threshold,
          '\nFitted 99.9% Sample Bias: ',results_LM[i,'sambias_f'],
          'Fitted 99.9% MC Var in mlns GBP: ',results_LM[i,'mcvar_f'],
          '\nPredicted 99.9% Sample Bias: ',results_LM[i,'sambias_t'],
          'Predicted 99.9% MC Var in mlns GBP: ',results_LM[i,'mcvar_t'],
          '\n',sep = '\n')
    }
    
    
  }
  
  results_GN <- data.frame(results_GN)
  results_LM <- data.frame(results_LM)
  
  return(list(results_GN = results_GN ,results_LM = results_LM))
  
}


truncate_s2_tail <- function(data,theta0 = NULL,epsilon0 = NULL,tail_cutoffs = NULL, 
                             mquants = NULL,v = 1,gamma = NULL, rhohat = NULL,
                             w2 = NULL,fittype = 'A', period = 6.5,rad0 = 500,
                             stoptype = 'step',tol = 1e-10, eta = 0.2,hold_weight_ratio = TRUE,
                             nIter_vector = NULL, print_cut_Iter = TRUE)
{
  if ( is.null(mquants)) {  mquants <- trunc(length(data)/4)}
  
  init_mquants <- mquants
  
  if (is.null(tail_cutoffs)) { 
    tail_cutoffs <- seq(0,0.8,by=0.2) 
    
    nIter_vector <- seq(10000,20000,length.out = length(tail_cutoffs)) 
    }
  
  
  if ( is.null(theta0)) { theta0 <- c(0,150,1.1) }
  
  if ( !is.null(gamma) )
  {
    v <- NULL
  }

  if ( is.null(w2) ) { w2 <- rep(1,mquants) }

  init_results <- stage2(data,theta0 = theta0,v = v,gamma = gamma,rhohat = rhohat, 
                         epsilon0 = epsilon0,period = period,cut_pct = tail_cutoffs[1],w2 = w2,
                         fittype = fittype, mquants = init_mquants,rad0 = rad0, stoptype = stoptype,
                         tol = tol, maxIter = nIter_vector[1], eta = eta,pvar = 1 - 0.001/( 1 - tail_cutoff[1]), printIter = TRUE)

  # if ( init_results$fullresults$gof > 0.068) { stop('The initial fit failed the GoF test. Try increasing
  #                         the number of iterations for the initial estimates.') }

  thresholds <- unname(quantile(data,tail_cutoffs))

  theo_params <- as.matrix(t(apply( as.matrix(thresholds),MARGIN = 1,
                                    function(t) theo_tail_cut_param(init_results$par,t))))

  results <- matrix(0,nrow = length(tail_cutoffs),ncol = 13)


  colnames(results) <-  c('cuts',
                          'mu_f','sigma_f','xi_f','mcvar_f','gof_f','sambias_f',
                          'mu_t','sigma_t','xi_t','mcvar_t','gof_t','sambias_t')

  results[,'cuts'] <- tail_cutoffs
  results[1,c('gof_t','gof_f')] <- init_results$fullresults$gof
  results[1,c('sambias_t','sambias_f')] <- init_results$fullresults$sambias
  results[,'sambias_t'] <- init_results$fullresults$sambias
  results[1,c('mu_f','sigma_f','xi_f')] <- theo_params[1,]
  results[,c('mu_t','sigma_t','xi_t')] <- theo_params[,]
  results[1,c('mcvar_f','mcvar_t')] <- init_results$mcvar

  theta_c <- theo_params[1,]

  for ( i in 2:length(tail_cutoffs) )
    {
      cut <- tail_cutoffs[i]
      cut_threshold <- thresholds[i]
      theo_param <- theo_params[i,]
      cut_data <- data[data > cut_threshold]
      
      cut_mquants <- trunc((length(data[data > cut_threshold])/length(data))*init_mquants)
      
      if (!is.null(rhohat)) 
      {
        
        cut_rhohat <- tail( rhohat , cut_mquants ) 
        
        tailresults <- stage2(cut_data,mquants = cut_mquants,theta0 = theta0,epsilon0 = NULL,period = period,
                              rhohat = cut_rhohat,cut_pct = 0,v = v,gamma = gamma, w2 = w2, fittype = fittype,
                              rad0 = rad0,stoptype = stoptype,tol = tol, maxIter = nIter_vector[i],
                              eta = eta, pvar = 1 - 0.001/(1-cut),printIter = FALSE)
        }

      # pvar is re-normalised in the function, no need to re-normalise here.
      else 
      {
        tailresults <- stage2(cut_data,theta0 = theta0,epsilon0 = NULL,period = period,rhohat = rhohat,
                              cut_pct = 0,v = v,gamma = gamma, w2 = NULL, fittype = fittype, 
                              mquants = cut_mquants,
                              rad0 = rad0,stoptype = stoptype,tol = tol, maxIter = nIter_vector[i],
                              eta = eta, pvar = 1 - 0.001/(1-cut),printIter = FALSE)
      }
      
      theta_c <- tailresults$par
      
      cut_init_epsilon <- tail( init_results$fullresults$epsilon , cut_mquants)

      results[i,c('mu_f','sigma_f','xi_f')] <- tailresults$par

      results[i,'gof_f'] <- tailresults$fullresults$gof
      results[i,'gof_t'] <- tailresults$trapzgof( theo_param,cut_init_epsilon)

      results[i,'mcvar_f'] <- tailresults$mcvar
      results[i,'mcvar_t'] <- fnvar(data,pgpd,rgpd,period,cut_threshold,100000,theo_param,1 - 0.001/(1 - cut))

      results[i,'sambias_f'] <- tailresults$fullresults$sambias

      if ( print_cut_Iter == TRUE)
        {
          cat('\n\nCutting the first ',cut*100,'% of the datas set off from ',cut_threshold,
              '\nFitted 99.9% Sample Bias: ',results[i,'sambias_f'],
              'Fitted 99.9% MC Var in mlns GBP: ',results[i,'mcvar_f'],
              'Fitted GoF Stat.: ',results[i,'gof_f'],
              '\nPredicted 99.9% Sample Bias: ',results[i,'sambias_t'],
              'Predicted 99.9% MC Var in mlns GBP: ',results[i,'mcvar_t'],
              'Predicted GoF Stat.: ',results[i,'gof_t'],
              '\n',sep = '\n')
        }

    }

  results <- data.frame(results)

  return(results)
  
}



predict_s2_tail <- function(data,theta0 = NULL,epsilon0 = NULL,tail_cutoffs = NULL, 
                             mquants = NULL,v = 1,gamma = NULL, rhohat = NULL,
                             w2 = NULL,fittype = 'A', period = 6.5,rad0 = 500,
                             stoptype = 'step',tol = 1e-10, eta = 0.2,hold_weight_ratio = TRUE,
                             nIter_vector = NULL, print_cut_Iter = TRUE)
{
  # Similar to truncated_s2_tail but does not fit at every cutoff, just extrapolates using the 
  # stability condition.
  if ( is.null(mquants)) {  mquants <- trunc(length(data)/4)}
  
  init_mquants <- mquants
  
  if (is.null(tail_cutoffs)) { 
    tail_cutoffs <- seq(0,0.8,by=0.2) 
    
    nIter_vector <- seq(10000,20000,length.out = length(tail_cutoffs)) 
    }
  
  
  if ( is.null(theta0)) { theta0 <- c(0,150,1.1) }
  
  if ( !is.null(gamma) )
  {
    v <- NULL
  }
  
  if ( is.null(w2) ) { w2 <- rep(1,mquants) }
  
  init_results <- stage2(data,theta0 = theta0,v = v,gamma = gamma,rhohat = rhohat, 
                         epsilon0 = epsilon0,period = period,cut_pct = tail_cutoffs[1],w2 = w2,
                         fittype = fittype, mquants = init_mquants,rad0 = rad0, stoptype = stoptype,
                         tol = tol, maxIter = nIter_vector[1], eta = eta,pvar = 0.999, printIter = TRUE)
  
  # if ( init_results$fullresults$gof > 0.068) { stop('The initial fit failed the GoF test. Try increasing
  #                         the number of iterations for the initial estimates.') }
  
  thresholds <- unname(quantile(data,tail_cutoffs))
  
  theo_params <- as.matrix(t(apply( as.matrix(thresholds),MARGIN = 1,
                                    function(t) theo_tail_cut_param(init_results$par,t))))
  
  results <- matrix(0,nrow = length(tail_cutoffs),ncol = 7)
  
  
  colnames(results) <-  c('cuts','mu','sigma','xi','mcvar','gof','sambias')
  
  results[,'cuts'] <- tail_cutoffs
  results[1,'gof'] <- init_results$fullresults$gof
  results[,'sambias'] <- init_results$fullresults$sambias
  results[,c('mu','sigma','xi')] <- theo_params[,]
  results[1,'mcvar'] <- init_results$mcvar

  for ( i in 2:length(tail_cutoffs) )
    {
      cut <- tail_cutoffs[i]
      cut_threshold <- thresholds[i]
      cut_data <- data[data > cut_threshold]
      cut_mquants <- trunc((length(cut_data)/length(data))*init_mquants)
      
      theo_param <- theo_params[i,]
      
      epsilon <- tail(init_results$fullresults$epsilon, cut_mquants)
      
      cut_processed_data <- processdata( cut_data ,fittype, cut_mquants )
      
      #X <- cut_processed_data$X
      Fn <- seq(1/length(cut_data),1,length.out = length(cut_data))
      
      Fx <- {1 - (1 + (theo_param[3]/theo_param[2])*(cut_data - theo_param[1]))^(-1/theo_param[3]) }
      
      cut_gof <- round(trapz( x = Fn, y = abs( Fn - Fx )  ),10)
      
      cut_mcvar <- fnvar(data,pgpd,rgpd,period,cut_threshold,100000,theo_param,1 - 0.001/(1 - cut))
      
      results[i,'gof'] <- cut_gof
      results[i,'mcvar'] <- cut_mcvar
      
      if ( print_cut_Iter == TRUE)
        {
          cat(
              '\n\nCutting the first ',cut*100,'% of the datas set off from ',cut_threshold,
              '\nGoF Stat.: ',cut_gof,
              'Predicted 99.9% MC Var in mlns GBP: ',results[i,'mcvar'],
              '\n',sep = '\n'
              )
        }
      
    }
  
  results <- data.frame(results)
  
  return(results)
  

}


get_ODR_xi_range <- function(data,theta0 = NULL,epsilon0 = NULL,tail_cutoffs = NULL, shift_data = FALSE,
                             mquants = NULL,v = 1,gamma = NULL, rhohat = NULL,
                             w2 = NULL,fittype = 'A', period = 6.5,rad0 = 500,
                             stoptype = 'step',tol = 1e-10, eta = 0.2,hold_weight_ratio = TRUE,
                             nIter_vector = NULL, print_cut_Iter = TRUE)
  {
    if ( is.null(mquants)) {  mquants <- trunc(length(data)/4)}
    
    init_mquants <- mquants
    
    if (is.null(tail_cutoffs)) { 
      tail_cutoffs <- seq(0,0.8,by=0.2) 
      
      nIter_vector <- seq(10000,20000,length.out = length(tail_cutoffs)) 
      }
    
    if ( is.null(theta0)) { theta0 <- c(0,150,1.1) }
    
    if ( !is.null(gamma) )
    {
      v <- NULL
    }
    
    if ( is.null(w2) ) { w2 <- rep(1,mquants) }
    
    init_results <- stage2(data,theta0 = theta0,v = v,gamma = gamma,rhohat = rhohat, run_VaR = FALSE,
                           epsilon0 = epsilon0,period = period,cut_pct = tail_cutoffs[1],w2 = w2,
                           fittype = fittype, mquants = init_mquants,rad0 = rad0, stoptype = stoptype,
                           tol = tol, maxIter = nIter_vector[1], eta = eta,pvar = 0.999, printIter = TRUE)
    
    
    thresholds <- unname(quantile(data,tail_cutoffs))
    
    results <- matrix(0,nrow = length(tail_cutoffs),ncol = 2)
    
    colnames(results) <-  c('cuts','xi')
    
    results[,'cuts'] <- tail_cutoffs
    results[1,'xi'] <- init_results$par[3]
    
    for ( i in 2:length(tail_cutoffs) )
      {
        cut <- tail_cutoffs[i]
        cut_threshold <- thresholds[i]
        cut_data <- data[data > cut_threshold]
        
        cut_mquants <- trunc((length(data[data > cut_threshold])/length(data))*init_mquants)
        
        if (!is.null(rhohat)) 
        {
          
          cut_rhohat <- tail( rhohat , cut_mquants ) 
          
          tailresults <- stage2(cut_data,mquants = cut_mquants,theta0 = theta0,epsilon0 = NULL,period = period,
                                rhohat = cut_rhohat,cut_pct = 0,v = v,gamma = gamma, w2 = NULL, fittype = fittype,
                                rad0 = rad0,stoptype = stoptype,tol = tol, maxIter = nIter_vector[i],run_VaR = FALSE,
                                eta = eta, pvar = 1 - 0.001/(1-cut),printIter = FALSE)
          }
  
        else 
        {
          tailresults <- stage2(cut_data,theta0 = theta0,epsilon0 = NULL,period = period,rhohat = rhohat,
                                cut_pct = 0,v = v,gamma = gamma, w2 = NULL, fittype = fittype, 
                                mquants = cut_mquants,run_VaR = FALSE,
                                rad0 = rad0,stoptype = stoptype,tol = tol, maxIter = nIter_vector[i],
                                eta = eta, pvar = 1 - 0.001/(1-cut),printIter = FALSE)
        }
        
  
        results[i,'xi'] <- tailresults$par[3]
  
        if ( print_cut_Iter == TRUE)
          {
            cat('\n\nCutting the first ',cut*100,'% of the datas set off from ',cut_threshold,
                '\nXi: ',results[i,'xi'] ,
                '\n',sep = '\n')
          }
  
      }
  
    results <- data.frame(results)
  
    return(results)
    
  }

predict_pickands_tail <- function( data, printiter = TRUE )
  {
    n <- length(data)
    N <- seq.int(1,n,1)
    data <- sort(data)
    revdata <- sort(data, decreasing = TRUE)
    
    x75 <- revdata[ seq(2,trunc(length(revdata)/4),1) ]
    x50 <- revdata[ seq(3,2*trunc(length(revdata)/4),2) ]
    x0 <- revdata[ seq(5,4*trunc(length(revdata)/4),4) ]
    
    pickands_results <- xi_pickands(data)
    
    xi <- pickands_results$xi
    sigma <- pickands_results$sigma
    gof_fitted <- pickands_results$gof
    
    thresholds <- x0
    results <- matrix(0,nrow = length(thresholds), ncol = 6)
    colnames(results) <-  c('k','mu','sigma','xi','mean_mse','mean_gof')
    
    results[,c('mu','sigma','xi')] <- cbind(x0,sigma,xi)
    results[,'k'] <- N
    
    for (t in 1:length(thresholds))
    {
      cut_data <- data[data > thresholds[t]]
      cut_n <- length(cut_data) -1
      cut_ecdf <- ecdf(cut_data)
      
      cat('Cuttin first ', round(1- cut_n/n,4),
            '\n',sep = '\n')
      
      theta_t <- c( thresholds[t], sigma[t], xi[t])
      
      params <- theo_tail_cut_param( theta_t, cut_data )
      
      mse_t <- rep(0,cut_n)
      gof_t <- rep(0,cut_n)
      
      for ( i in seq(1,cut_n,by = 4))
      {
        
        if ( params[i,2] <= 0) 
          { 
            mse_t[i] <- NaN
            gof_t[i] <- NaN
            next
        }
        
        mse_t[i] <- mse( pgpd( cut_data[cut_data > cut_data[i]] , params[i, 1], params[i,2], params[i,3]  ), 
                          cut_ecdf(cut_data[cut_data > cut_data[i]] ) )
        
        gof_t[i] <- trapzgof( cut_data[ cut_data > cut_data[i] ] ,  params[i,] )
        
      }
      
      results[t,'mean_mse'] <- mean(mse_t)
      results[t,'mean_gof'] <- mean(gof_t)
      
    }
    
    return(results)
    
  }

  
######################################### Cut-off Estimators ###################################

mean_excess <- function(data,theta)
{
  # The fitted vector is given by the analytical continuous conditional expection E[ X - t | X > t ]
  # under GPD.
  data <- sort(data)
  fitted <- ( theta[2] + theta[3]*( data - theta[1] ) )/( 1 - theta[3] )
  
  # The empirical vector is given by the normalised empirical mean by the no. observations taken.
  N <- seq(1,length(data) ,1)
  get_mean_excess <- function(t) { return( sum((data[data > data[t]] - data[t]))/length(data[data>data[t]]) ) }
  empirical <- apply(as.matrix(N),MARGIN = 1, FUN = get_mean_excess)
  
  return(data.frame(k = N,fitted = fitted, empirical = empirical))
}

d2f <- function(x,y)
{
  n <- length(x)
  dx <- 1/n
  d2ydx2 <- (y[3:n] - y[2:(n-1)]- ( y[2:(n-1)] - y[1:(n-2)]))/(dx^2)
  d2ydx2[n-2] <- NA
  return( c(0,0,d2ydx2) )
}








