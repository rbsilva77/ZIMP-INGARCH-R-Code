#-----------------------------------------
# R code for fitting ZIPIG-INGARCH model #
#-----------------------------------------

fit_zipigingarch <- function(x, order, method){
  
  # Installing necessary R packages
  pacman::p_load(tseries, gamlss.dist) 
  
  n <- length(x) # sample size
  
  #--------------------------------------------------------
  # Auxiliary functions for Poisson inverse-Gaussian case
  qsi0 <- -1 / 2
  B <- function(x) {-sqrt((-2 * x))}
  d2B <- function(x) {(2 * sqrt(2) * (-x)^(3 / 2))^(-1)}
  T <- function(x) {x}
  D <- function(x) {log(x) / 2}
  dD <- function(x) {1 / (2 * x)}
  H <- function(x) {-log(2 * pi * x^3) / 2}
  G <- function(x) {-1 / (2 * x)}
  
  # Logit function
  logit <- function(z){log(z / (1 - z))}
  
  #--------------------------------------------------------
  
  # Conditional expectations
  
  # E(B|X)
  kappa_t <- function(x){
    
    y <- x[1:n]
    mu <- x[(n + 1):(2 * n)]
    phi <- x[2 * n + 1]
    omega <- x[2 * n + 2]
    
    EB_X <- ifelse(y == 0, (1 - omega) * dPIG(y, mu, 1 / phi) / (omega + (1 - omega) * dPIG(y, mu, 1 / phi)), 1)
    
    return(EB_X)
  }
  
  # E(BZ|X)
  eta_t <- function(x){
    
    y <- x[1:n]
    mu <- x[(n + 1):(2 * n)]
    phi <- x[2 * n + 1]
    omega <- x[2 * n + 2]
    
    EZ_X <- ifelse(y == 0, (1 - omega) * (y + 1) * dPIG(y + 1, mu, 1 / phi) / (mu * (omega + (1 - omega) * dPIG(y, mu, 1 / phi))),
                   (y + 1) * dPIG(y + 1, mu, 1 / phi) / (mu * dPIG(y, mu, 1 / phi))) 
    
    return(EZ_X)
  }
  
  # E(Bg(Z)|X)
  nu_t <- function(x){
    
    y <- x[1:n]
    mu <- x[(n + 1):(2 * n)]
    phi <- x[(2 * n + 1)]
    omega <- x[2 * n + 2]
    
    yaux <- ifelse(y > 0, y - 1, y)
    
    EBgZ_X <- ifelse(y == 0, (1 - omega) * dPIG(y, mu, 1 / phi) / (omega + (1 - omega) * dPIG(y, mu, 1 / phi)) *
                       ((-0.5) * (1 + sqrt(phi * (2 * mu + phi))) / phi), (-0.5) * mu * dPIG(yaux, mu, 1 / phi) / (y * dPIG(y, mu, 1 / phi)))
    
    return(EBgZ_X)
  }
  
  #--------------------------------------------------------
  
  # R code for fitting ZIPIG-INGARCH(1,0) model
  
  fit_zipig.ingarch1 <- function(x, method){
    
    #--------------------------------------------------------------------
    # Quasi likelihood estimation (QMLE) for ZIPIG-INGARCH(1,0) processes
    
    # QMLE(1,0)
    fit.qmle1 <- function(x, method){  
      
      # Quasi log-likelihood function
      loglik.qmle <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        
        # current lambda
        lambda <- NULL
        lambda[1] <- alpha0 / (1 - alpha1)
        
        for (i in 1:(n-1)){
          lambda[i+1] <- alpha0 + alpha1 * x[i]
        }
        
        # likelihood function
        sum(x * log(lambda) - lambda)
      }
      
      # Gradient vector
      grad.qmle <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        
        # current lambda
        lambda <- NULL
        lambda[1] <- alpha0 / (1 - alpha1)
        
        for (i in 1:(n-1)){
          lambda[i+1] <- alpha0 + alpha1 * x[i]
        }
        
        # Gradient functions
        grad1 <- function(par){
          
          x <- par[1:n]
          lambda <- par[(n+1):(2*n)]
          
          sum(x[1:n] / lambda[1:n] - 1)
        }
        
        grad2 <- function(par){
          
          x <- par[1:n]
          lambda <- par[(n+1):(2*n)]
          
          sum((x[2:n] / lambda[2:n] - 1) * x[1:(n-1)]) 
        }
        
        c(grad1(c(x, lambda)), grad2(c(x, lambda)))
      }
      
      # Formulation of the constraints: ui %*% theta >= ci,
      # where theta = (alpha0, alpha1)
      ui <- diag(1,2)
      ci <- c(0, 0)
      fit<- try(constrOptim(theta = c(.1, .1), f = loglik.qmle, grad = grad.qmle, 
                            ui = ui, ci = ci, method = method, control = list(fnscale = -1)))
      
      alpha0_qmle <- fit$par[1]
      alpha1_qmle <- fit$par[2]
      
      # Estimated lambda (lambda.qmle)
      lambda_qmle <- NULL
      lambda_qmle[1] <- alpha0_qmle / (1 - alpha1_qmle)
      
      for (i in 1:(n-1)){
        lambda_qmle[i+1] <- alpha0_qmle + alpha1_qmle * x[i]
      }
      
      phi_qmle <- 1 / (sum(((x - lambda_qmle)^2 - lambda_qmle) / lambda_qmle^2) / n)
      
      est_qmle <- c(alpha0_qmle, alpha1_qmle, phi_qmle)
      
      return(list("est_qmle" = est_qmle))
    }
    #--------------------------------------------------------------------
    
    # Tolerance - stopping criteria
    tol <- 0.0001
    
    zero_prop <- length(which(x == 0)) / n # proportion of zeros
    
    #----- Initial values (QMLE estimation) -----#
    
    qmle_par <- suppressWarnings(fit.qmle1(x, method = method)$est_qmle)
    
    alpha0_qmle <- qmle_par[1]
    alpha1_qmle <- qmle_par[2]
    phi_qmle <- qmle_par[3]
    
    alpha0old <- alpha0_qmle
    alpha1old <- alpha1_qmle  
    phiold <- phi_qmle
    omegaold <- zero_prop
    
    #----- EM algorithm -----#
    
    # Begin of the estimation process
    
    tolerance <- 1
    
    while(tolerance > tol){
      
      # Previous mu
      muold <- NULL
      muold[1] <- (1 - omegaold) * alpha0old / (1 - (1 - omegaold) * alpha1old)
      
      for (i in 1:(n-1)){
        muold[i+1] <- alpha0old + alpha1old * x[i]
      }
      
      # Gradient vector
      grad.EM <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        
        # current mu
        mu <- NULL
        mu[1] <- (1 - omegaold) * alpha0 / (1 - (1 - omegaold) * alpha1)
        for (i in 1:(n-1)){
          mu[i+1] <- alpha0 + alpha1 * x[i]
        }
        
        kappa <- kappa_t(c(x, muold, phiold, omegaold))
        eta <- eta_t(c(x, muold, phiold, omegaold))
        
        # grad functions
        grad1 <- function(par){
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum(kappa * x / mu - eta)
        }
        
        grad2 <- function(par){
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum((x[2:n] * kappa[2:n] / mu[2:n] - eta[2:n]) * x[1:(n-1)])
        }
        
        c(grad1(c(x, mu, muold)), grad2(c(x, mu, muold)))
      }
      
      # Expected likelihood function   
      loglik.EM <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        
        # Current mu
        mu <- NULL
        mu[1] <- (1 - omegaold) * alpha0 / (1 - (1 - omegaold) * alpha1)
        
        for (i in 1:(n-1)){
          mu[i+1] <- alpha0 + alpha1 * x[i]
        }
        
        kappa <- kappa_t(c(x, muold, phiold, omegaold)) # E(B|X)
        
        eta <- eta_t(c(x, muold, phiold, omegaold)) # E(BZ|X)
        
        nu <- nu_t(c(x, muold, phiold, omegaold)) # E(Bg(Z)|X)
        
        # Expected likelihood function
        sum(kappa * x * log(mu) - mu * eta + kappa * D(phiold) +
              phiold * (eta * qsi0 - kappa * B(qsi0) + nu) + kappa * logit(1-omegaold) + log(omegaold))
      }
      
      # Formulation of the constraints: ui %*% theta >= ci,
      # where theta = (alpha0, alpha1)
      ui <- diag(1, 2)
      ci <- c(0, 0)
      fit <- try(constrOptim(theta = c(alpha0old, alpha1old), f = loglik.EM, grad = grad.EM, 
                             ui = ui, ci = ci, method = method,
                             control = list(fnscale = -1)))
      
      # Error handling
      if(class(fit) != "try-error"){
        
        alpha0new <- fit$par[1]
        alpha1new <- fit$par[2]
        
        omeganew <- 1 - sum(kappa_t(c(x, muold, phiold, omegaold))) / (n-1) 
        
        phinew <- sum(kappa_t(c(x, muold, phiold, omegaold))) / 
          (2 * sum(kappa_t(c(x, muold, phiold, omegaold)) * B(qsi0) - 
                     eta_t(c(x, muold, phiold, omegaold)) * qsi0 - nu_t(c(x, muold, phiold, omegaold))))
        
        # tolerance criterion
        tolerance <- max(abs(c(alpha0new, alpha1new) - c(alpha0old, alpha1old)))
        
        # updating parameters values
        alpha0old <- alpha0new
        alpha1old <- alpha1new
        phiold <- phinew
        omegaold <- omeganew
      }
      
      if(class(fit) == "try-error"){
        tolerance = tol
      }
      
    } # End of the estimation process
    
    #----- Results -----#
    
    thetahat <- round(c(alpha0new, alpha1new, phinew, omeganew), 4)
    
    # Estimated mu
    muhat <- NULL
    muhat[1] <- (1 - omeganew) * alpha0new / (1 - (1 - omeganew) * alpha1new)  
    
    for (i in 2:n){
      muhat[i] <- alpha0new + alpha1new * x[i-1]
    }
    
    loglik <- sum(log(dZIPIG(x, mu = muhat, sigma = 1 / phinew, nu = omeganew)))
    
    # AIC: -2 * l(thetahat) + 2 * (p+q+1)
    AIC <- -2 * loglik + (1 + 0 + 1) * 2
    
    # BIC: -2 * l(thetahat) + (p+q+1) * log(n)
    BIC <- -2 * loglik + (1 + 0 + 1) * log(n)
    
    #----- Output -----#
    
    return(list("Parameters" = thetahat, "Likelihood" = loglik, "AIC" = AIC, "BIC" = BIC))
    
  } 
  # End of the fit_zipig_ingarch1 function
  
  #----------------------------------------------------------------------------------
  
  # R code for fitting ZIPIG-INGARCH(2,0) model
  
  fit_zipig.ingarch2 <- function(x, method){
    
    #--------------------------------------------------------------------
    # Quasi likelihood estimation (QMLE) for ZIPIG-INGARCH(2,0) processes
    
    # QMLE(2,0)
    fit_qmle2 <- function(x, method){
      
      # Quasi-log-likelihood function
      qloglik.qmle <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        alpha2 <- z[3]
        
        # current lambda
        lambda <- NULL
        lambda[1] <- alpha0 / (1 - alpha1 - alpha2)
        lambda[2] <- alpha0 + alpha1 * x[1]
        for (i in 2:(n-1)){
          lambda[i+1] <- alpha0 + alpha1 * x[i] + alpha2 * x[i-1]
        }
        
        # likelihood function
        sum(x * log(lambda) - lambda)
      }
      
      # Gradient vector
      grad.qmle <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        alpha2 <- z[3]
        
        # current lambda
        lambda <- NULL
        lambda[1] <- alpha0 / (1 - alpha1 - alpha2)
        lambda[2] <- alpha0 + alpha1 * x[1]
        for (i in 2:(n-1)){
          lambda[i+1] <- alpha0 + alpha1 * x[i] + alpha2 * x[i-1]
        }
        
        # Gradient functions
        grad1 <- function(par){
          x <- par[1:n]
          lambda <- par[(n+1):(2*n)]
          
          sum(x[1:n] / lambda[1:n] - 1)
        }
        
        grad2 <- function(par){
          x <- par[1:n]
          lambda <- par[(n+1):(2*n)]
          
          sum((x[2:n] / lambda[2:n] - 1) * x[1:(n-1)]) 
        }
        
        grad3 <- function(par){
          x <- par[1:n]
          lambda <- par[(n+1):(2*n)]
          sum((x[3:n] / lambda[3:n] - 1) * x[1:(n-2)]) 
        }
        
        c(grad1(c(x, lambda)), grad2(c(x, lambda)), grad3(c(x, lambda)))
      }
      
      # Formulation of the constraints: ui %*% theta >= ci,
      # where theta = (alpha0, alpha1)
      ui <- diag(1,3)
      ci <- c(0, 0, 0)
      fit.qmle <- try(constrOptim(theta = c(.1, .1, .1), f = qloglik.qmle,
                                  grad = grad.qmle, ui = ui, ci = ci,
                                  method = method,
                                  control = list(fnscale = -1)))
      
      alpha0.qmle <- fit.qmle$par[1]
      alpha1.qmle <- fit.qmle$par[2]
      alpha2.qmle <- fit.qmle$par[3]
      
      # Estimated lambda (lambda.qmle)
      lambda.qmle <- NULL
      lambda.qmle[1] <- alpha0.qmle / (1 - alpha1.qmle)
      lambda.qmle[2] <- alpha0.qmle + alpha1.qmle * x[1]
      
      for (i in 2:(n-1)){
        lambda.qmle[i+1] <- alpha0.qmle + alpha1.qmle * x[i] + alpha2.qmle * x[i-1]
      }
      
      phi.qmle <- 1 / (sum(((x - lambda.qmle)^2 - lambda.qmle) / lambda.qmle^2) / n)
      
      qmle_est <- c(alpha0.qmle, alpha1.qmle, alpha2.qmle, phi.qmle)
      
      return(list("qmle_est" = qmle_est))
    } 
    #--------------------------------------------------------------------
    
    # Tolerance - stopping criteria
    tol <- 0.0001
    
    zero_prop <- length(which(x == 0)) / n # proportion of zeros
    
    #----- Initial values (QMLE estimation) -----#
    
    qmle_res <- suppressWarnings(fit_qmle2(x, method = method)$qmle_est)
    
    alpha0.qmle <- qmle_res[1]
    alpha1.qmle <- qmle_res[2]
    alpha2.qmle <- qmle_res[3]
    phi.qmle <- qmle_res[4]
    
    alpha0old <- alpha0.qmle
    alpha1old <- alpha1.qmle
    alpha2old <- alpha2.qmle  
    phiold <- phi.qmle
    omegaold <- zero_prop
    
    #----- EM algorithm -----#
    
    # Begin of the estimation process
    
    tolerance <- 1
    
    while(tolerance > tol){
      
      # Previous mu
      muold <- NULL
      muold[1] <- (1 - omegaold) * alpha0old / (1 - (1 - omegaold) * (alpha1old + alpha2old))
      muold[2] <- alpha0old + alpha1old * x[1]
      
      for (i in 2:(n-1)){
        muold[i+1] <- alpha0old + alpha1old * x[i] + alpha2old * x[i-1]
      }
      
      # Gradient vector
      grad.EM <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        alpha2 <- z[3]
        
        # Current mu
        mu <- NULL
        mu[1] <- (1 - omegaold) * alpha0 / (1 - (1 - omegaold) * (alpha1 + alpha2))
        mu[2] <- alpha0 + alpha1 * x[1] 
        
        for (i in 2:(n-1)){
          mu[i+1] <- alpha0 + alpha1 * x[i] + alpha2 * x[i-1]
        }
        
        kappa <- kappa_t(c(x, muold, phiold, omegaold))
        eta <- eta_t(c(x, muold, phiold, omegaold))
        
        # grad functions
        grad1 <- function(par){
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum(kappa * x / mu - eta)
        }
        
        grad2 <- function(par){
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum((x[2:n] * kappa[2:n] / mu[2:n] - eta[2:n]) * x[1:(n-1)])
        }
        
        grad3 <- function(par){
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum((x[3:n] * kappa[3:n] / mu[3:n] - eta[3:n]) * x[1:(n-2)])
        }
        
        c(grad1(c(x, mu, muold)), grad2(c(x, mu, muold)), grad3(c(x, mu, muold)))
      }
      
      # Expected likelihood function   
      loglik.EM <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        alpha2 <- z[3]
        
        # Current mu
        mu <- NULL
        mu[1] <- (1-omegaold)*alpha0 / (1 - (1-omegaold)*(alpha1+alpha2))
        mu[2] <- alpha0 + alpha1 * x[1]
        
        for (i in 2:(n-1)){
          mu[i+1] <- alpha0 + alpha1 * x[i] + alpha2 * x[i-1]
        }
        
        kappa <- kappa_t(c(x, muold, phiold, omegaold)) # E(B|X)
        eta <- eta_t(c(x, muold, phiold, omegaold)) # E(BZ|X)
        nu <- nu_t(c(x, muold, phiold, omegaold)) # E(Bg(Z)|X)
        
        # Expected likelihood function
        sum(kappa * x * log(mu) - mu * eta + kappa * D(phiold) +
              phiold * (eta * qsi0 - kappa * B(qsi0) + nu) + kappa * logit(1 - omegaold) + log(omegaold))
      }
      
      # Formulation of the constraints: ui %*% theta >= ci,
      # where theta = (alpha0, alpha1)
      ui <- diag(1,3)
      ci <- c(0, 0, 0)
      fit <- try(constrOptim(theta = c(alpha0old, alpha1old, alpha2old), f = loglik.EM,
                             grad = grad.EM, ui = ui, ci = ci, method = method,
                             control = list(fnscale = -1)))
      
      # Error handling
      if(class(fit) != "try-error"){
        
        alpha0new <- fit$par[1]
        alpha1new <- fit$par[2]
        alpha2new <- fit$par[3]
        
        omeganew <- 1 - sum(kappa_t(c(x, muold, phiold, omegaold))) / (n - 1) 
        
        phinew <- sum(kappa_t(c(x, muold, phiold, omegaold))) / 
          (2 * sum(kappa_t(c(x, muold, phiold, omegaold)) * B(qsi0) - 
                     eta_t(c(x, muold, phiold, omegaold)) * qsi0 - nu_t(c(x, muold, phiold, omegaold))))
        
        # tolerance criterion
        tolerance <- max(abs(c(alpha0new, alpha1new, alpha2new) - c(alpha0old, alpha1old, alpha2old)))
        
        # updating parameters values
        alpha0old <- alpha0new
        alpha1old <- alpha1new
        alpha2old <- alpha2new
        phiold <- phinew
        omegaold <- omeganew
      }
      
      if(class(fit) == "try-error"){
        tolerance = tol
      } 
    } # End of the estimation process
    
    #----- Results -----#
    
    thetahat <- round(c(alpha0new, alpha1new, alpha2new, phinew, omeganew), 4)
    
    # Estimated mu
    muhat <- NULL
    muhat[1] <- (1 - omeganew) * alpha0new / (1 - (1 - omeganew) * (alpha1new + alpha2new))
    muhat[2] <- alpha0new + alpha1new * x[1]
    
    for (i in 2:(n-1)){
      muhat[i+1] <- alpha0new + alpha1new * x[i] + alpha2new * x[i-1]
    }
    
    loglik <- sum(log(dZIPIG(x, mu = muhat, sigma = 1/phinew, nu = omeganew)))
    
    # AIC: -2 * l(thetahat) + 2 * (p + q + 1)
    AIC <- -2 * loglik + (2 + 0 + 1) * 2
    
    # BIC: -2 * l(thetahat) + (p + q + 1) * log(n)
    BIC <- -2 * loglik + (2 + 0 + 1) * log(n)
    
    #----- Output -----#
    
    return(list("Parameters" = thetahat, "Likelihood" = loglik, "AIC" = AIC, "BIC" = BIC))
    
  } 
  # End of the fit_zipig_ingarch2 function
  
  #----------------------------------------------------------------------------------
  
  # R code for fitting ZIPIG-INGARCH(1,1) model
  
  fit_zipig.ingarch11 <- function(x, method){
    
    #--------------------------------------------------------------------
    # Quasi likelihood estimation (QMLE) for ZIPIG-INGARCH(2,0) processes
    
    # QMLE(1,1)
    fit.qmle11 <- function(x, method){
      
      # Quasi-log-likelihood function
      qloglik.qmle <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        beta1 <- z[3]
        
        # current lambda
        lambda <- NULL
        lambda[1] <- alpha0 / (1 - alpha1 - beta1)
        for (i in 1:(n-1)){
          lambda[i+1] <- alpha0 + alpha1 * x[i] + beta1 * lambda[i]
        } 
        # likelihood function
        sum(x * log(lambda) - lambda)
      }
      
      # Auxiliary functions:
      
      # First derivative of lambda_t with respect to alpha_0
      dlbd.a0 <- function(z){
        a1 <- z[1]
        a2 <- z[2]
        b1 <- z[3]
        dmut <- NULL
        dmut[1] <- 1 / (1 - a1 - a2 - b1)^2
        for(i in 2:n){
          dmut[i] <- 1 + b1 * dmut[i-1]
        }  
        dmut
      }
      
      # First derivative of lambda_t with respect to alpha_1
      dlbd.a1 <- function(z){
        a0 <- z[1]
        a1 <- z[2]
        a2 <- z[3]
        b1 <- z[4]
        dmut <- NULL
        dmut[1] <- a0 / (1 - a1 - a2 - b1)^2
        for(i in 2:n){
          dmut[i] <- x[i-1] + b1 * dmut[i-1]
        }  
        dmut
      }
      
      # Gradient vector
      grad.qmle <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        beta1 <- z[3]
        
        # current lambda
        lambda <- NULL
        lambda[1] <- alpha0 / (1 - alpha1 - beta1)
        for (i in 2:n){
          lambda[i] <- alpha0 + alpha1 * x[i-1] + beta1 * lambda[i-1]
        }
        
        # grad functions
        grad1 <- function(par){
          x <- par[1:n]
          lambda <- par[(n+1):(2*n)]
          sum((x[2:n] / lambda[2:n] - 1) * dlbd.a0(c(alpha1, 0, beta1))[1:(n-1)])
        }
        
        grad2 <- function(par){
          x <- par[1:n]
          lambda <- par[(n+1):(2*n)]
          sum((x[2:n] / lambda[2:n] - 1) * dlbd.a1(c(alpha0, alpha1, 0, beta1))[1:(n-1)]) 
        }
        
        # First derivative of lambda_t with respect to beta1
        dlbd.b1 <- function(z){
          a0 <- z[1]
          a1 <- z[2]
          b1 <- z[3]
          dmut <- NULL
          dmut[1] <- a0 / (1 - a1 - b1)^2
          for(i in 2:n){
            dmut[i] <- lambda[i-1] + b1 * dmut[i-1]
          }  
          dmut
        }
        
        grad3 <- function(par){
          x <- par[1:n]
          lambda <- par[(n+1):(2*n)]
          sum((x[2:n] / lambda[2:n] - 1) * dlbd.b1(c(alpha0, alpha1, beta1))[1:(n-1)]) 
        }
        
        c(grad1(c(x, lambda)), grad2(c(x, lambda)), grad3(c(x, lambda)))
      }
      
      # Formulation of the constraints: ui %*% theta >= ci,
      # where theta = (alpha0, alpha1)
      ui <- diag(1,3)
      ci <- c(0, 0, 0)
      fit.qmle <- try(constrOptim(theta = c(.1, .1, .1), f = qloglik.qmle, grad = grad.qmle, ui = ui, ci = ci,
                                  method = method, control = list(fnscale = -1)))
      
      alpha0.qmle <- fit.qmle$par[1]
      alpha1.qmle <- fit.qmle$par[2]
      beta1.qmle <- fit.qmle$par[3]
      
      # Estimated lambda (lambda.qmle)
      lambda.qmle <- NULL
      lambda.qmle[1] <- alpha0.qmle / (1 - alpha1.qmle - beta1.qmle)
      
      for(i in 2:n){
        lambda.qmle[i] <- alpha0.qmle + alpha1.qmle * x[i-1] + beta1.qmle * lambda.qmle[i-1]
      }
      
      phi.qmle <- 1 / (sum(((x - lambda.qmle)^2 - lambda.qmle) / lambda.qmle^2) / n)
      
      qmle_est <- c(alpha0.qmle, alpha1.qmle, beta1.qmle, phi.qmle)
      
      return(list("qmle_est" = qmle_est))
    }# End of QMLE(1,1)
    #-------------------------------------------------------------------- 
    
    # Tolerance - stopping criteria
    tol <- 0.0001
    
    zero_prop <- length(which(x == 0)) / n # proportion of zeros
    
    #----- Initial values (QMLE estimation) -----#
    
    qmle_fit <- suppressWarnings(fit.qmle11(x, method = method)$qmle_est)
    
    alpha0.qmle <- qmle_fit[1]
    alpha1.qmle <- qmle_fit[2]
    beta1.qmle <- qmle_fit[3]
    phi.qmle <- qmle_fit[4]
  
    alpha0old <- alpha0.qmle
    alpha1old <- alpha1.qmle
    beta1old <- beta1.qmle
    phiold <- phi.qmle
    omegaold <- zero_prop
    
    #----- EM algorithm -----#
    
    # Begin of the estimation process
    
    tolerance <- 1
    
    while(tolerance > tol){
      
      # Previous mu
      muold <- NULL
      muold[1] <- (1-omegaold) * alpha0old / (1 - (1-omegaold) * alpha1old - beta1old)  
      for (i in 2:n){
        muold[i] <- alpha0old + alpha1old * x[i-1] + beta1old * muold[i-1]
      }
      
      # first derivative of mu_t with respect to alpha_0
      dmut.a0 <- function(z){
        a1 <- z[1]
        b1 <- z[2]
        dmut <- array(0, c(n,1))
        dmut[1] <- (1-omegaold)*(1-(1-omegaold)*a1-b1)^(-1)
        dmut[2:n] <- 1 + b1 * dmut[1:(n-1)]
        dmut
      }
      # first derivative of mu_t with respect to alpha_1
      dmut.a1 <- function(z){
        a0 <- z[1]
        a1 <- z[2]
        b1 <- z[3]
        dmut <- array(0, c(n,1))
        dmut[1] <- (1-omegaold)^2 * a0 * (1-(1-omegaold)*a1-b1)^(-2)
        dmut[2:n] <- x[1:(n-1)] + b1 * dmut[1:(n-1)]
        dmut
      }
      # first derivative of mu_t with respect to beta1
      dmut.b1 <- function(z){
        a0 <- z[1]
        a1 <- z[2]
        b1 <- z[3]
        dmut <- array(0, c(n,1))
        dmut[1] <- (1-omegaold) * a0 * (1-(1-omegaold)*a1-b1)^(-2)
        dmut[2:n] <- muold[1:(n-1)] + b1 * dmut[1:(n-1)]
        dmut
      }
      
      # Gradient vector
      grad.EM <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        beta1 <- z[3]
        
        # Current mu
        mu <- NULL
        mu[1] <- (1-omegaold) * alpha0 / (1 - (1-omegaold) * alpha1 - beta1)
        for (i in 2:n){
          mu[i] <- alpha0 + alpha1 * x[i-1] + beta1*mu[i-1]
        }
        
        kappa <- kappa_t(c(x, muold, phiold, omegaold))
        eta <- eta_t(c(x, muold, phiold, omegaold))
        
        # Grad functions
        grad1 <- function(par){
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum((kappa * x / mu - eta) * dmut.a0(c(alpha1, beta1))[1:n])
        }
        
        grad2 <- function(par){
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum((x[2:n]*kappa[2:n] / mu[2:n] - eta[2:n]) * 
                dmut.a1(c(alpha0, alpha1, beta1))[1:n])
        }
        
        grad3 <- function(par){
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum((x[2:n]*kappa[2:n] / mu[2:n] - eta[2:n]) *
                dmut.b1(c(alpha0, alpha1, beta1))[1:n])
        }
        
        c(grad1(c(x, mu, muold)), grad2(c(x, mu, muold)), grad3(c(x,mu,muold)))
      }
      
      # expected likelihood function   
      loglik.EM <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        beta1 <- z[3]
        
        # Current mu
        mu <- NULL
        mu[1] <- (1 - omegaold) * alpha0 / (1 - (1 - omegaold) * alpha1 - beta1)
        for (i in 2:n){
          mu[i] <- alpha0 + alpha1 * x[i-1] + beta1 * mu[i-1]
        }
        
        kappa <- kappa_t(c(x, muold, phiold, omegaold)) # EB\X
        eta <- eta_t(c(x, muold, phiold, omegaold)) # EZ\X
        nu <- nu_t(c(x, muold, phiold, omegaold)) # EgZ\X
        
        # Expected likelihood function
        sum(kappa * x * log(mu) - mu * eta + kappa * D(phiold) +
              phiold * (eta * qsi0 - kappa * B(qsi0) + nu) + kappa * logit(1-omegaold) + log(omegaold))
      }
      
      # Formulation of the constraints: ui %*% theta >= ci,
      # where theta = (alpha0, alpha1)
      ui <- diag(1,3)
      ci <- c(0, 0, 0)
      fit <- suppressWarnings(try(constrOptim(theta = c(alpha0old, alpha1old, beta1old), 
                                              f = loglik.EM, grad = grad.EM, ui = ui, ci = ci, 
                                              method = method, control = list(fnscale = -1))))
      
      if(class(fit) != "try-error"){
        
        parhat <- fit$par
        
        alpha0new <- parhat[1]
        alpha1new <- parhat[2]
        beta1new <- parhat[3]
        
        omeganew <- 1 - sum(kappa_t(c(x, muold, phiold, omegaold)))/(n-1) 
        
        phinew <- sum(kappa_t(c(x, muold, phiold, omegaold))) / 
          (2*sum(kappa_t(c(x, muold, phiold, omegaold))*B(qsi0) - 
                   eta_t(c(x, muold, phiold, omegaold))*qsi0 - nu_t(c(x, muold, phiold, omegaold))))
        
        # Tolerance criterion
        tolerance <- max(abs(c(alpha0new, alpha1new, beta1new) - c(alpha0old, alpha1old, beta1old)))
        
        # Updating parameters values
        alpha0old <- alpha0new
        alpha1old <- alpha1new
        beta1old <- beta1new
        phiold <- phinew
        omegaold <- omeganew
      }
      
      if(class(fit) == "try-error"){
        tolerance = tol
        stop("Error in the estimation process!")
      } 
      
    } # End of the estimation process
    
    #----- Results -----#
    
    thetahat <- round(c(alpha0new, alpha1new, beta1new, phinew, omeganew), 4)
    
    # Estimated mu
    muhat <- NULL
    muhat[1] <- (1 - omeganew) * alpha0new / (1 - (1 - omeganew) * alpha1new - beta1new)
    for (i in 2:n){
      muhat[i] <- alpha0new + alpha1new * x[i-1] + beta1new * muhat[i-1]
    }
    
    loglik <- sum(log(dZIPIG(x, mu = muhat, sigma = 1 / phinew, nu = omeganew)))
    
    # AIC: -2 * l(thetahat) + 2 * (p + q + 1)
    AIC <- -2 * loglik + 2 * (1 + 1 + 1)
    
    # BIC: -2 * l(thetahat) + (p + q + 1) * log(n)
    BIC <- -2 * loglik + (1 + 1 + 1) * log(n)
    
    #----- Output -----#
    
    return(list("Parameters" = thetahat, "Likelihood" = loglik, "AIC" = AIC, "BIC" = BIC))
    
  } 
  # End of the fit_zipig_ingarch11 function
  
  #----------------------------------------------------------------------------------
  
  order <- as.character(order)
  
  switch(order,
         "1" = fit_zipig.ingarch1(x, method),
         "2" = fit_zipig.ingarch2(x, method),
         "11" = fit_zipig.ingarch11(x, method),
         stop("Order not expected!"))
  
} 
# End of the fit_zipig.ingarch function
