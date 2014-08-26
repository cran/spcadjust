### R code from vignette source 'spcadjust-intro.Rnw'

###################################################
### code chunk number 1: spcadjust-intro.Rnw:44-45
###################################################
library(spcadjust)


###################################################
### code chunk number 2: spcadjust-intro.Rnw:70-71
###################################################
chart <- new("SPCCUSUMNormal",Delta=1);


###################################################
### code chunk number 3: spcadjust-intro.Rnw:76-79
###################################################
X <-  rnorm(250)
xihat <- xiofdata(chart,X)
str(xihat)


###################################################
### code chunk number 4: spcadjust-intro.Rnw:83-84
###################################################
plot(runchart(chart, newdata=rnorm(100),xi=xihat),ylab=expression(S[t]),xlab="t",type="b")


###################################################
### code chunk number 5: spcadjust-intro.Rnw:88-89
###################################################
plot(runchart(chart, newdata=rnorm(100,mean=c(rep(0,50),rep(1,50))),xi=xihat),ylab=expression(S[t]),xlab="t",type="b")


###################################################
### code chunk number 6: spcadjust-intro.Rnw:105-109
###################################################
SPCproperty(data=X,nrep=50,
            property=new("calARLCUSUM",chart=chart,target=100))
SPCproperty(data=X,nrep=50,
            property=new("calhitprobCUSUM",chart=chart,target=0.05,nsteps=1000))


###################################################
### code chunk number 7: spcadjust-intro.Rnw:113-119
###################################################
SPCproperty(dat=X,nrep=50,
            property=new("ARLCUSUM",chart=chart,threshold=3),
            covprob=c(0.8,0.9))
SPCproperty(dat=X,nrep=50,
            property=new("hitprobCUSUM",chart=chart,threshold=5,nsteps=100),
            covprob=c(0.8,0.9))


###################################################
### code chunk number 8: spcadjust-intro.Rnw:125-128
###################################################
chartnp <- new("SPCCUSUMNonparCenterScale",Delta=1)
SPCproperty(data=X,
            nrep=100,property=new("calARLCUSUM",chart=chartnp,target=100))


###################################################
### code chunk number 9: spcadjust-intro.Rnw:142-143
###################################################
chartShew <- new("SPCShewNormalCenterScale")


###################################################
### code chunk number 10: spcadjust-intro.Rnw:147-148
###################################################
plot(runchart(chartShew, newdata=rnorm(100),xi=xiofdata(chart,X)),ylab=expression(S[t]),xlab="t")


###################################################
### code chunk number 11: spcadjust-intro.Rnw:152-163
###################################################
SPCproperty(data=X,nrep=100,
            property=new("calARLShew",chart=chartShew,target=100))

SPCproperty(data=X,nrep=500,
            property=new("ARLShew",chart=chartShew,threshold=4))

SPCproperty(data=X,nrep=500,
            property=new("hitprobShew",chart=chartShew,nsteps=100,threshold=4))

SPCproperty(data=X,nrep=500,
            property=new("calhitprobShew",chart=chartShew,target=0.01,nsteps=100))


###################################################
### code chunk number 12: spcadjust-intro.Rnw:188-191
###################################################
n <- 1000
Xlinreg <- data.frame(x1= rbinom(n,1,0.4), x2= runif(n,0,1), x3= rnorm(n))
Xlinreg$y <- 2 + Xlinreg$x1 + Xlinreg$x2 + Xlinreg$x3 + rnorm(n)


###################################################
### code chunk number 13: spcadjust-intro.Rnw:195-199
###################################################
chartlinreg <- new("SPCCUSUMlm",Delta=1,formula="y~x1+x2+x3")
SPCproperty(data=Xlinreg,
            nrep=100,
            property=new("calARLCUSUM",chart=chartlinreg,target=100))


###################################################
### code chunk number 14: spcadjust-intro.Rnw:202-207
###################################################
chartlinreg <- new("SPCCUSUMlm",Delta=1,formula="y~x1")
SPCproperty(data=Xlinreg,
            nrep=100,
            property=new("calARLCUSUM",chart=chartlinreg,target=100))



###################################################
### code chunk number 15: spcadjust-intro.Rnw:235-239
###################################################
n <- 1000
Xlogreg <- data.frame(x1=rbinom(n,1,0.4), x2=runif(n,0,1), x3=rnorm(n))
xbeta <- -1+Xlogreg$x1*100+Xlogreg$x2+Xlogreg$x3
Xlogreg$y <- rbinom(n,1,exp(xbeta)/(1+exp(xbeta)))


###################################################
### code chunk number 16: spcadjust-intro.Rnw:243-247
###################################################
chartlogreg <- new("SPCCUSUMlogreg",Delta= 1, formula="y~x1+x2+x3")
SPCproperty(data=Xlogreg,
            nrep=100,
            property=new("calARLCUSUM",chart=chartlogreg,target=100))


###################################################
### code chunk number 17: spcadjust-intro.Rnw:271-276
###################################################
setClass("SPCCUSUMNormalROBUST", contains="SPCCUSUMNormal")
setMethod("Pofdata", "SPCCUSUMNormalROBUST",
          function(chart,data){
              list(mu= median(data), sd= mad(data), m=length(data))
          })


###################################################
### code chunk number 18: spcadjust-intro.Rnw:279-283
###################################################
X <-  rnorm(100)
chartrobust <- new("SPCCUSUMNormalROBUST",Delta=1)
SPCproperty(data=X,nrep=50,
            property=new("calARLCUSUM",chart=chartrobust,target=100))


###################################################
### code chunk number 19: spcadjust-intro.Rnw:290-312
###################################################
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


###################################################
### code chunk number 20: spcadjust-intro.Rnw:317-320
###################################################
X <- rexp(1000)
plot(runchart(ExpCUSUMchart, newdata=rexp(100),xi=xiofdata(ExpCUSUMchart,X)),
     ylab=expression(S[t]),xlab="t",type="b")


###################################################
### code chunk number 21: spcadjust-intro.Rnw:324-336
###################################################
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


###################################################
### code chunk number 22: spcadjust-intro.Rnw:394-395 (eval = FALSE)
###################################################
## vignette("spcadjust-intro")


