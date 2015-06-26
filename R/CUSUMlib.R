hitprob_sim <- function(c,nrep=1000,n=1e3,robs){
  mean(replicate(nrep,{
    R <- cumsum(robs(n))
    S <- R-cummin(R)
    any(S>=c)
  }))
}

calibratehitprob_sim <- function(hprob=0.01,nrep=1000,n=1e3,robs){
  x <- (replicate(nrep,{
    R <- cumsum(robs(n))
    S <- R-cummin(R)
    max(S)
  }))
  quantile(x,1-hprob)
}



ARL_sim <- function(c,nrep=1000,maxsteps=1e4,robs){
  mean(replicate(nrep,{
    R <- cumsum(robs(maxsteps))
    S <- R-cummin(R)
    w <- S>=c
    if (any(w))
      return(min(which(w)))
    else
      return(maxsteps*1.5)
  }))
}

getQ <- function(c,gridpoints,pobs){
  p <- c(0,(0.5+(0:(gridpoints-2)))*c/(gridpoints-1))
  ptarget <- p+c/(gridpoints-1)/2
  sapply(p,function(x) {res <- pobs(ptarget-x); c(res[1],diff(res))})
}

ARL_Markovapprox <- function(c,gridpoints=100,pobs){
  Q <- getQ(c,gridpoints,pobs)
  tryCatch(1+rep(1,gridpoints)%*%solve(diag(rep(1,gridpoints))-Q,c(1,rep(0,gridpoints-1))),error=function(e) Inf)
}


calibrateARL_Markovapprox<- function(ARL=1000,f=ARL_Markovapprox,pobs,...){
  cmax <- 1;
  while(f(pobs=pobs,cmax,...)<ARL) cmax <- cmax*2
  cmin <- 1
  while(f(pobs=pobs,cmin,...)>ARL&&cmin>1e-10) cmin <- cmin/2
  if (cmin<=1e-10){
    warning(paste("cmin=",cmin,"\n"))
    stop()
  }
  uniroot(function(x) f(pobs=pobs,c=x,...)-ARL,lower=cmin,upper=cmax)$root
}


matrix.power <- function(mat, n)
{
  # test if mat is a square matrix
  # treat n < 0 and n = 0 -- this is left as an exercise
  # trap non-integer n and return an error
  if (n == 1) return(mat)
  result <- diag(1, ncol(mat))
  while (n > 0) {
    if (n %% 2 != 0) {
      result <- result %*% mat
      n <- n - 1
    }
    mat <- mat %*% mat
    n <- n / 2
  }
  return(result)
}


hitprob_Markovapprox <- function(pobs,c,n,gridpoints=100){
  Q <- getQ(c,gridpoints=gridpoints,pobs=pobs)
  Q <- matrix.power(Q,n)
  1-rep(1,gridpoints)%*%Q%*%c(1,rep(0,gridpoints-1))
}


calibratehitprob_Markovapprox<- function(hprob=0.01,n,f=hitprob_Markovapprox,pobs,...){
  cmax <- 1;
  while(f(pobs=pobs,n=n,cmax,...)>hprob) cmax <- cmax*2
  cmin <- 1
  while(f(pobs=pobs,n=n,cmin,...)<hprob&&cmin>1e-10) cmin <- cmin/2
  if (cmin<=1e-10){
    warning(paste("cmin=",cmin,"\n"))
    stop()
  }
  uniroot(function(x) f(pobs=pobs,n=n,c=x,...)-hprob,lower=cmin,upper=cmax)$root
}


