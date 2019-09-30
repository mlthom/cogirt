#Newton-Raphson Algorythm for posterior estimates of theta
newtonrhaphson_post <- function(cand_theta, tau, alpha, y,  mu, sigma, ftn_post, tol=1e-9, max.iter=100, verbose=F, points=F, fischer=F){
  map_theta <- cand_theta
  fx <- ftn_post(map_theta, tau, alpha, y, mu, sigma, fischer)
  iter <- 0
  if(points==T){text(map_theta[1],map_theta[2],labels=format(iter),col="purple")}
  while(any(abs(fx$fpder) > tol) && (iter < max.iter)){
    map_theta <- map_theta - fx$fpder %*% solve(fx$spder)
    if(any(abs(map_theta)>(diag(sigma)*5))){ #bounded solution to help convergence
      map_theta[which(abs(map_theta)>(diag(sigma)*5))] <- (rnorm(length(map_theta)) + sign(map_theta)*(diag(sigma)*5))[which(abs(map_theta)>(diag(sigma)*5))]
    }
    fx <- ftn_post(map_theta, tau, alpha, y, mu, sigma, fischer)
    iter <- iter + 1
    #Sys.sleep(.1)
    if(verbose==T){cat("\r","At iteration", iter, "MAP theta is", format(round(map_theta,3),nsmall=3), " logLik is", format(round(m2plloglikSum(map_theta,tau,alpha,y),6),nsmall=6),"")}
    if(points==T){text(map_theta[1],map_theta[2],labels=format(iter),col="purple")}
  }
  if(verbose==T){
    if(any(abs(fx$fpder) > tol)){
      cat("Algorithm failed to converge\n")
    } else {
      cat("Algorithm converged\n")
    }
  }
  return(list(map_theta,ifelse(fischer,-1,1)*fx$spder)) #note that I have to multiply by -1 if I using Fisher scoring
}
