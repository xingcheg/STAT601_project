######### For binomial random component and Pregibon’s one-parameter link
gi <- function(b,x) b[1]+x*b[2]

mui <- function(la,b,x){
  gi <- gi(b,x)
  result <- 1 - ( 1 + la*exp(gi) )^(-1/la)
  return(result)
}

g_mui <- function(la,b,x){
  mui <- mui(la,b,x)
  result <- la/( (1-mui) * (1-(1-mui)^la) )
  return(result)
}

Vi <- function(la,b,x){
  mui <- mui(la,b,x)
  result <- mui*(1-mui)
  return(result)
}

wi <- function(la,b,x,m){
  Vi <- Vi(la,b,x)
  g_mui <- g_mui(la,b,x)
  result <- m * (Vi * (g_mui)^2)^(-1)
  return(result)
}

g_lai <- function(la,b,x){
  mui <- mui(la,b,x)
  result <- (-log(1-mui))/(1-(1-mui)^la) - 1/la
  return(result)
}


zi <- function(la,b,x,y){
  mui <- mui(la,b,x)
  g_mui <- g_mui(la,b,x)
  result <- (y-mui)*g_mui
  return(result)
}

W <- function(la,b,X,M){
  n <- length(X)
  w <- NULL
  for (i in 1:n){
    w <- c(w, wi(la,b,X[i],M[i]))
  }
  return(diag(w))
}

Z <- function(la,b,X,Y){
  n <- length(X)
  z <- NULL
  for (i in 1:n){
    z <- c(z, zi(la,b,X[i],Y[i]))
  }
  return(z)
}

XA <- function(la,b,X){
  n <- length(X)
  Xa <- NULL
  for (i in 1:n){
    xa <- c(1,X[i],-g_lai(la,b,X[i]))
    Xa <- rbind(Xa, xa)
  }
  return(Xa)
}


grad <- function(la,b,X,Y,M){
  XA <- XA(la,b,X)
  W <- W(la,b,X,M)
  Z <- Z(la,b,X,Y)
  return(t(XA)%*%W%*%Z)
}

Hei <- function(la,b,X,M){
  XA <- XA(la,b,X)
  W <- W(la,b,X,M)
  return(t(XA)%*%W%*%XA)
}



log_lik <- function(la,b,X,Y,M){
  n <- length(X)
  sum <- 0
  for (i in 1:n){
    mu <- mui(la,b,X[i])
    theta <- log(mu/(1-mu))
    b_theta <- -log(1-mu)
    li <- M[i]*(Y[i]*theta-b_theta)
    sum <- sum + li
  }
  return(sum)
}



sat_log_lik <- function(Y,M){
  n <- length(Y)
  sum <- 0
  for (i in 1:n){
    mu <- Y[i]
    theta <- log(mu/(1-mu))
    b_theta <- -log(1-mu)
    li <- M[i]*(Y[i]*theta-b_theta)
    sum <- sum + li
  }
  return(sum)
}



dev_resid <- function(la,b,X,Y,M){
  n <- length(X)
  resid <- NULL
  for (i in 1:n){
    mu <- mui(la,b,X[i])
    theta <- log(mu/(1-mu))
    b_theta <- -log(1-mu)
    l_mu <- M[i]*(Y[i]*theta-b_theta)

    theta_y <- log(Y[i]/(1-Y[i]))
    b_theta_y <- -log(1-Y[i])
    l_y <- M[i]*(Y[i]*theta_y-b_theta_y)
    
    r <- sign(Y[i]-mu)*sqrt(-2*(l_mu-l_y))
    resid <- c(resid,r)
  }
  return(resid)
}



para_glm_binom <- function(X,Y,M, startb, startla){
  p <- c(startb,startla)
  iter <- 0
  while(1){
    b <- p[1:2]
    la <- p[3]
    p1 <- p + solve(Hei(la,b,X,M)) %*% grad(la,b,X,Y,M)
    iter <- iter +1
    cat("iteration =", iter, "\n")
    cat("b0 = ", p1[1],",b1 = ", p1[2], ",lambda = ", p1[3], "\n")
    if( norm(p1-p, type = "F") < 1e-4 || norm(grad(p1[3],p1[1:2],X,Y,M), type = "F")<1e-4 ) break
    p <- p1
  }
  b1 <- p1[1:2]
  la1 <- p1[3]
  cat("\n\n")
  cat("MLE : b0 = ",b1[1],",b1 = ",b1[2],",lambda = ",la1,"\n\n")
  inv_Hei1 <- solve( Hei(la1,b1,X,M) )
  cat("Inverse Information:\n")
  print(inv_Hei1)
  loglik <- log_lik(la1,b1,X,Y,M)
  cat("\n Log-likelihood = ", loglik, "\n\n")
  sat_loglik <- sat_log_lik(Y,M)
  cat("\n Saturate Log-likelihood = ", sat_loglik, "\n\n")
  dev_res <- dev_resid(la1,b1,X,Y,M)
  Output <- list(estb = b1, estla = la1, inv_Hei = inv_Hei1, 
                 loglik = loglik, sat_loglik = sat_loglik, dev_res = dev_res)
  return(Output)
}









