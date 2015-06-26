## ------------------------------------------------------------------------
library(spcadjust)

## ------------------------------------------------------------------------
chart <- new("SPCCUSUMNormal",Delta=1);

## ------------------------------------------------------------------------
X <-  rnorm(250)
xihat <- xiofdata(chart,X)
str(xihat)

## ----fig.width=8,fig.height=4--------------------------------------------
plot(runchart(chart, newdata=rnorm(100),xi=xihat),ylab=expression(S[t]),xlab="t",type="b")

## ----fig.width=8,fig.height=4--------------------------------------------
plot(runchart(chart, newdata=rnorm(100,mean=c(rep(0,50),rep(1,50))),xi=xihat),ylab=expression(S[t]),xlab="t",type="b")

## ------------------------------------------------------------------------
SPCproperty(data=X,nrep=50,
            property=new("calARLCUSUM",chart=chart,target=100))
SPCproperty(data=X,nrep=50,
            property=new("calhitprobCUSUM",chart=chart,target=0.05,nsteps=1000))

## ------------------------------------------------------------------------
SPCproperty(dat=X,nrep=50,
            property=new("ARLCUSUM",chart=chart,threshold=3),
            covprob=c(0.8,0.9))
SPCproperty(dat=X,nrep=50,
            property=new("hitprobCUSUM",chart=chart,threshold=5,nsteps=100),
            covprob=c(0.8,0.9))

## ------------------------------------------------------------------------
cal <- SPCproperty(data=X,nrep=50,
            property=new("calARLCUSUM",chart=chart,target=100))
newX <- rnorm(100)
S <- runchart(chart, newdata=newX,xi=xihat)

## ----fig.width=8,fig.height=4--------------------------------------------
plot(S,ylab=expression(S[t]),xlab="t",type="b",ylim=range(S,cal@res+1,cal@raw))
lines(c(0,100),rep(cal@res,2),col="red")
lines(c(0,100),rep(cal@raw,2),col="blue")
legend("topleft",c("Adjusted Threshold","Unadjusted Threshold"),col=c("red","blue"),lty=1)

## ------------------------------------------------------------------------
chartnp <- new("SPCCUSUMNonparCenterScale",Delta=1)
SPCproperty(data=X,
            nrep=100,property=new("calARLCUSUM",chart=chartnp,target=100))

## ------------------------------------------------------------------------
chartShew <- new("SPCShewNormalCenterScale",twosided=TRUE)

## ------------------------------------------------------------------------
X <-  rnorm(250)
xihat <- xiofdata(chartShew,X)

## ----fig.width=8,fig.height=4--------------------------------------------
plot(runchart(chartShew, newdata=X,xi=xihat),ylab=expression(S[t]),xlab="t")

## ------------------------------------------------------------------------
SPCproperty(data=X,nrep=100,
            property=new("calARLShew",chart=chartShew,target=741))

SPCproperty(data=X,nrep=100,
            property=new("ARLShew",chart=chartShew,threshold=3))

SPCproperty(data=X,nrep=100,
            property=new("hitprobShew",chart=chartShew,nsteps=100,threshold=3))

SPCproperty(data=X,nrep=100,
            property=new("calhitprobShew",chart=chartShew,target=0.01,nsteps=100))

## ------------------------------------------------------------------------
cal <- SPCproperty(data=X,nrep=100,
            property=new("calARLShew",chart=chartShew,target=741))
S <- runchart(chartShew, newdata=newX,xi=xihat)

## ----fig.width=8,fig.height=4--------------------------------------------
plot(S,ylab=expression(S[t]),xlab="t",type="b",ylim=range(S,cal@res+2,cal@raw,-cal@res-1,-cal@raw))
lines(c(0,100),rep(cal@res,2),col="red")
lines(c(0,100),rep(cal@raw,2),col="blue")
lines(c(0,100),-rep(cal@res,2),col="red")
lines(c(0,100),-rep(cal@raw,2),col="blue")
legend("topleft",c("Adjusted Threshold","Unadjusted Threshold"),col=c("red","blue"),lty=1)

## ------------------------------------------------------------------------
n <- 1000
Xlinreg <- data.frame(x1= rbinom(n,1,0.4), x2= runif(n,0,1), x3= rnorm(n))
Xlinreg$y <- 2 + Xlinreg$x1 + Xlinreg$x2 + Xlinreg$x3 + rnorm(n)

## ------------------------------------------------------------------------
chartlinreg <- new("SPCCUSUMlm",Delta=1,formula="y~x1+x2+x3")
SPCproperty(data=Xlinreg,
            nrep=100,
            property=new("calARLCUSUM",chart=chartlinreg,target=100))

## ------------------------------------------------------------------------
chartlinreg <- new("SPCCUSUMlm",Delta=1,formula="y~x1")
SPCproperty(data=Xlinreg,
            nrep=100,
            property=new("calARLCUSUM",chart=chartlinreg,target=100))


## ------------------------------------------------------------------------
xihat <- xiofdata(chartlinreg,Xlinreg)
cal <- SPCproperty(data=Xlinreg,
             nrep=100,
             property=new("calARLCUSUM",chart=chartlinreg,target=100))
n <- 100
newXlinreg <- data.frame(x1= rbinom(n,1,0.4), x2= runif(n,0,1), x3=rnorm(n))
newXlinreg$y <- 2 + newXlinreg$x1 + newXlinreg$x2 + newXlinreg$x3 + rnorm(n)+c(rep(0,50),rep(1,50))
S <- runchart(chartlinreg, newdata=newXlinreg,xi=xihat)

## ----fig.width=8,fig.height=4--------------------------------------------
plot(S,ylab=expression(S[t]),xlab="t",type="b",ylim=range(S,cal@res+1,cal@raw))
lines(c(0,100),rep(cal@res,2),col="red")
lines(c(0,100),rep(cal@raw,2),col="blue")
legend("topleft",c("Adjusted Threshold","Unadjusted Threshold"),col=c("red","blue"),lty=1)

## ------------------------------------------------------------------------
n <- 1000
Xlogreg <- data.frame(x1=rbinom(n,1,0.4), x2=runif(n,0,1), x3=rnorm(n))
xbeta <- -1+Xlogreg$x1*100+Xlogreg$x2+Xlogreg$x3
Xlogreg$y <- rbinom(n,1,exp(xbeta)/(1+exp(xbeta)))

## ------------------------------------------------------------------------
chartlogreg <- new("SPCCUSUMlogreg",Delta= 1, formula="y~x1+x2+x3")
SPCproperty(data=Xlogreg,
            nrep=100,
            property=new("calARLCUSUM",chart=chartlogreg,target=100))

## ----results='hide'------------------------------------------------------
setClass("SPCCUSUMNormalROBUST", contains="SPCCUSUMNormal")
setMethod("Pofdata", "SPCCUSUMNormalROBUST",
          function(chart,data){
              list(mu= median(data), sd= mad(data), m=length(data))
          })

## ------------------------------------------------------------------------
X <-  rnorm(100)
chartrobust <- new("SPCCUSUMNormalROBUST",Delta=1)
SPCproperty(data=X,nrep=50,
            property=new("calARLCUSUM",chart=chartrobust,target=100))

## ----results='hide'------------------------------------------------------
setClass("SPCCUSUMExponential",contains="SPCCUSUM",representation(Delta="numeric"))
setMethod("Pofdata", "SPCCUSUMExponential",
          function(chart,data){
              list(lambda=1/mean(data), n=length(data))
          })
setMethod("xiofP", "SPCCUSUMExponential",
          function(chart,P) P$lambda)
setMethod("resample", "SPCCUSUMExponential",
          function(chart,P) rexp(P$n,rate=P$lambda))
setMethod("getcdfupdates", "SPCCUSUMExponential",
          function(chart, P, xi) {
              ; function(x){ if(chart@Delta<1)
                pmax(0,1-exp(-P$lambda*(x-log(chart@Delta))/(xi*(1-chart@Delta))))
              else
                pmin(1,exp(-P$lambda*(log(chart@Delta)-x)/(xi*(chart@Delta-1))))
            }
          })
setMethod("updates","SPCCUSUMExponential",
          function(chart,xi,data) log(chart@Delta)-xi*(chart@Delta-1)*data
)

ExpCUSUMchart=new("SPCCUSUMExponential",Delta=1.25)

## ----fig.width=8,fig.height=4--------------------------------------------
X <- rexp(1000)
plot(runchart(ExpCUSUMchart, newdata=rexp(100),xi=xiofdata(ExpCUSUMchart,X)),
     ylab=expression(S[t]),xlab="t",type="b")

## ------------------------------------------------------------------------
SPCproperty(data=X,
            nrep=100,
            property=new("hitprobCUSUM",chart=ExpCUSUMchart,
              threshold=1,nsteps=100),covprob=c(0.5,0.9))
SPCproperty(data=X,
            nrep=100,
            property=new("ARLCUSUM",chart=ExpCUSUMchart,
              threshold=3),covprob=c(0.5,0.9))
SPCproperty(data=X,
            nrep=100,
            property=new("calARLCUSUM",chart=ExpCUSUMchart,
              target=1000),covprob=c(0.5,0.9))

## ----eval=FALSE----------------------------------------------------------
#  vignette("spcadjust-intro")

