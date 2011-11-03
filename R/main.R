#definition of charts
setClass("SPCchart")

setGeneric("xiofdata",def=function(chart,data){standardGeneric("xiofdata")})
setMethod("xiofdata", signature="SPCchart", function(chart,data){
    xiofP(chart,Pofdata(chart,data))
})

setGeneric("Pofdata",def=function(chart,data){standardGeneric("Pofdata")})
setGeneric("xiofP",def=function(chart,P){standardGeneric("xiofP")})
setGeneric("resample",def=function(chart,P){standardGeneric("resample")})


setGeneric("updates", def= function(chart,xi,data)  standardGeneric("updates"))
setGeneric("getcdfupdates", def= function(chart, P, xi) standardGeneric("getcdfupdates"))

setClass("SPCCUSUM",contains="SPCchart",representation("VIRTUAL"))
setGeneric("runchart",def=function(chart,newdata,xi){standardGeneric("runchart")})
setMethod("runchart", signature="SPCCUSUM", function(chart,newdata,xi){
    R <- cumsum(updates(chart,xi=xi, data=newdata))
    R - cummin(R)
})

############ #setting up CUSUM charts....
setClass("SPCCUSUMNormal", contains="SPCCUSUM",representation=list(Delta="numeric"))
setMethod("updates", signature="SPCCUSUMNormal",
          function(chart,xi,data) (data-xi$mu-chart@Delta/2)/xi$sd)
setMethod("Pofdata", signature="SPCCUSUMNormal",
          function(chart,data){
              list(mu= mean(data), sd= sd(data), m=length(data))
          })
setMethod("xiofP", signature="SPCCUSUMNormal",
          function(chart,P) P)
setMethod("resample", signature="SPCCUSUMNormal",
          function(chart,P) P$sd*rnorm(P$m)+P$mu)
setMethod("getcdfupdates", signature="SPCCUSUMNormal",
          function(chart, P, xi) {
              muupd <- (P$mu-xi$mu-chart@Delta/2)/xi$sd
              sdsqupd <- P$sd/xi$sd
              function(x) pnorm(x, mean=muupd, sd=sdsqupd)
          })

setClass("SPCCUSUMNonpar", contains="SPCCUSUM", representation("VIRTUAL"))
setMethod("Pofdata", signature="SPCCUSUMNonpar", function(chart,data) data)
setMethod("resample", signature="SPCCUSUMNonpar",
          function(chart,P){
              if (is.vector(P))
                  sample(P,replace=TRUE)
              else
                  P[sample.int(dim(P)[1],replace=TRUE),]
          })
setMethod("getcdfupdates", signature="SPCCUSUMNonpar",
          function(chart, P, xi) ecdf(updates(chart,xi=xi,data=P)))



setClass("SPCCUSUMNonparCenterScale", contains="SPCCUSUMNonpar",
         representation=list(Delta="numeric"))
setMethod("xiofP", signature="SPCCUSUMNonparCenterScale",
          function(chart,P) list(mu=mean(P),sd=sd(P)))
setMethod("updates", signature="SPCCUSUMNonparCenterScale",
          function(chart,xi,data) (data-xi$mu-chart@Delta/2)/xi$sd)

##### CUSUM chart for linear models
setClass("SPCCUSUMlm", contains="SPCCUSUMNonpar",
          representation=list(Delta="numeric",formula="character"))
setMethod("xiofP", "SPCCUSUMlm",
           function(chart,P)
           lm(chart@formula,data=P)
           )
setMethod("updates", "SPCCUSUMlm",
           function(chart,xi,data){
               response <-  model.response( model.frame( chart@formula,data=data))
               response - predict(xi,newdata=data) - chart@Delta/2
           }
           )

#### CUSUM chart for logistic regression
setClass("SPCCUSUMlogreg", contains="SPCCUSUMNonpar",
         representation(Delta="numeric",formula="character"))
setMethod("xiofP","SPCCUSUMlogreg",
          function(chart,P)
          glm(chart@formula,data=P,family=binomial("logit"))
          )
setMethod("updates", "SPCCUSUMlogreg",
          function(chart,xi,data){
              xbeta <- predict.glm(xi,newdata=data)
              response <-  model.response( model.frame( chart@formula,data=data))
              chart@Delta*response+ log(1+exp(xbeta))- log(1+exp(chart@Delta+xbeta))
            })



###############################
##properties - general definitions
setClass("SPCproperty",
         representation(chart="SPCchart",lowerconf="logical","VIRTUAL")
         )

setGeneric("SPCq", def=function(property, P, xi) standardGeneric("SPCq"))
setGeneric("qtrafo", def=function(property, x) standardGeneric("qtrafo"))
setGeneric("SPCoutput", def=function(property,  result) standardGeneric("SPCoutput"))



### general definition for ARL
setClass("SPCARL",
         representation(threshold="numeric","VIRTUAL"),
         prototype=list(lowerconf=FALSE),
         contains=c("SPCproperty")
         )
setMethod("SPCoutput", signature="SPCARL",
          function(property,result)
          paste("A threshold of  ", property@threshold, " gives an in-control ARL of at least ", format(result,digits=4), ".", sep="",collapse="")
          )
setMethod("qtrafo", signature="SPCARL", function(property,x) exp(x))

### general definition for hitting probabilities
setClass("SPChitprob",
         representation(threshold="numeric",nsteps="numeric","VIRTUAL"),
         prototype=list(lowerconf=TRUE),
         contains=c("SPCproperty"))
setMethod("SPCoutput", signature="SPChitprob",
          function(property,result)
          paste("A threshold of  ", property@threshold, " gives an in-control false alarm probability of at most ", format(result,digits=4), " within ",property@nsteps," steps.", sep="",collapse="")
          )
setMethod("qtrafo", signature="SPChitprob", function(property,x) exp(x)/(1+exp(x)))

##### general definitions for calibrating to a given ARL
setClass("SPCcalARL", representation(target="numeric","VIRTUAL"),
         prototype=list(lowerconf=TRUE),
         contains=c("SPCproperty"))
setMethod("qtrafo", signature="SPCcalARL", function(property,x) exp(x))
setMethod("SPCoutput", signature="SPCcalARL",
          function(property,result)
          paste("A threshold of ", format(result,digits=4), " gives an in-control ARL of at least ", property@target, ".", sep="",collapse="")
          )
##### general definitions for calibrating to a given hitting probability
setClass("SPCcalhitprob",
         representation(target="numeric",nsteps="numeric","VIRTUAL"),
         prototype=list(lowerconf=TRUE),
         contains=c("SPCproperty"))
setMethod("SPCoutput", signature="SPCcalhitprob",
          function(property,result)
          paste("A threshold of ", format(result,digits=4), " gives an in-control false alarm probability of at most ", property@target, " within ",property@nsteps, " steps.", sep="",collapse="")
          )
setMethod("qtrafo", signature="SPCcalhitprob", function(property,x) exp(x))


#############################
#properties for CUSUM charts
############################

#CUSUM calibrate ARL
setClass("calARLCUSUM", representation=list(gridpoints="numeric"),
         prototype=list(gridpoints=75),
         contains=c("SPCcalARL"))
setMethod("SPCq", signature="calARLCUSUM",
          function(property,P,xi)
          log(calibrateARL_Markovapprox(pobs=getcdfupdates(property@chart,xi=xi, P=P),ARL=property@target,gridpoints=property@gridpoints))
          )

#CUSUM calibrate hitting probability
setClass("calhitprobCUSUM",
         representation=list(gridpoints="numeric"),
         prototype=list(gridpoints=75),
         contains=c("SPCcalhitprob"))
setMethod("SPCq", signature="calhitprobCUSUM",
          function(property,P,xi)
          log(calibratehitprob_Markovapprox(pobs=getcdfupdates(property@chart,xi=xi, P=P),hprob=property@target,n=property@nsteps,gridpoints=property@gridpoints))
          )

#CUSUM - ARL
setClass("ARLCUSUM", representation=list(gridpoints="numeric",chart="SPCCUSUM"),
         prototype=list(gridpoints=75),
         contains=c("SPCARL"))
setMethod("SPCq", signature="ARLCUSUM",
          function(property,P,xi)
          as.double(log(ARL_Markovapprox(c=property@threshold,pobs=getcdfupdates(property@chart,xi=xi, P=P),gridpoints=property@gridpoints)))
          )

#CUSUM - hitting probability
setClass("hitprobCUSUM",
         representation=list(gridpoints="numeric"),
         prototype=list(gridpoints=75),
         contains=c("SPChitprob"))
setMethod("SPCq", signature="hitprobCUSUM",
          function(property,P,xi){
              res <- hitprob_Markovapprox(c=property@threshold,pobs=getcdfupdates(property@chart,xi=xi, P=P),n=property@nsteps,gridpoints=property@gridpoints);
              as.double(log(res/(1-res)))
          })


######
############ #definitions for Shewhart charts....

setClass("SPCShew",contains="SPCchart",representation="VIRTUAL")
setMethod("runchart", signature="SPCShew", function(chart,newdata,xi){
    updates(chart,xi=xi, data=newdata)
})

#parametric normal
setClass("SPCShewNormalCenterScale",contains="SPCShew")
setMethod("Pofdata", signature="SPCShewNormalCenterScale",
          function(chart,data){
              list(mu= mean(data), sd= sd(data), m=length(data))
          })
setMethod("xiofP", signature="SPCShewNormalCenterScale",
          function(chart,P) P)
setMethod("resample", signature="SPCShewNormalCenterScale",
          function(chart,P) P$sd*rnorm(P$m)+P$mu)
setMethod("getcdfupdates", signature="SPCShewNormalCenterScale",
          function(chart, P, xi) {
              muupd <- (P$mu-xi$mu)/xi$sd
              sdsqupd <- P$sd/xi$sd
              function(x) pnorm(x, mean=muupd, sd=sdsqupd)
          })
setMethod("updates",signature="SPCShewNormalCenterScale",
          function(chart,xi,data) (data-xi$mu)/xi$sd
)


#nonparametric center and scaling
setClass("SPCShewNonparCenterScale",contains="SPCShew")
setMethod("Pofdata", signature="SPCShewNonparCenterScale",
          function(chart,data) data
          )
setMethod("xiofP", signature="SPCShewNonparCenterScale",
          function(chart,P) list(mu=mean(P),sd=sd(P)))
setMethod("resample", signature="SPCShewNonparCenterScale",
          function(chart,P){
              if (is.vector(P))
                  sample(P,replace=TRUE)
              else
                  P[sample.int(dim(P)[1],replace=TRUE),]
          })
setMethod("getcdfupdates", signature="SPCShewNonparCenterScale",
          function(chart, P, xi) ecdf(updates(chart, xi=xi,data=P)))
setMethod("updates",signature="SPCShewNonparCenterScale",
          function(chart,xi,data) (data-xi$mu)/xi$sd)

#####################
#properties for Shewhart charts
#####################

###ARL for Shewhart charts
setClass("ARLShew", contains="SPCARL")
setMethod("SPCq", signature="ARLShew",
          function(property,P,xi)
           -log(1-getcdfupdates(property@chart, xi=xi, P=P)(property@threshold))
          )

#calibrating hitting probability for Shew charts
qShewcalibrateARL <- function(pobs,target){
    cmax <- 1;
    while((1/(1-pobs(cmax)))<target) cmax <- cmax*2
    cmin <- 1
    while((1/(1-pobs(cmin)))>target&&cmin>1e-10) cmin <- cmin/2
    if (cmin<=1e-10){
        warning(paste("cmin=",cmin,"\n"))
        stop()
    }
 uniroot(function(x) target-(1/(1-pobs(x))),lower=cmin,upper=cmax)$root
}
qShewlogcalibrateARL <- function(pobs, target) log(qShewcalibrateARL(pobs,target))
setClass("calARLShew",contains=c("SPCcalARL"))
setMethod("SPCq", signature="calARLShew",
          function(property,P,xi)
          qShewlogcalibrateARL(pobs=getcdfupdates(property@chart, xi=xi, P=P),target=property@target)
          )

###########computing hitting probabilties for Shewhart charts
setClass("hitprobShew", contains=c("SPChitprob"))
setMethod("SPCq", signature="hitprobShew",
          function(property,P,xi){
              survprob1step <- getcdfupdates(property@chart, xi=xi, P=P)(property@threshold)
              log(1-survprob1step^property@nsteps)-property@nsteps*log(survprob1step)
          }
          )

###########calibrating threshold for a desired  hitting probabilties for Shewhart charts
setClass("calhitprobShew",contains=c("SPCcalhitprob"))
qShewcalibratehitprob <- function(pobs,target,nsteps){
  cmax <- 1;
  while((1-pobs(cmax)^nsteps)>target) cmax <- cmax*2
  cmin <- 1
  while((1-pobs(cmin)^nsteps)<target&&cmin>1e-10) cmin <- cmin/2
  if (cmin<=1e-10){
    warning(paste("cmin=",cmin,"\n"))
    stop()
  }
 uniroot(function(x) target-(1-pobs(x)^nsteps),lower=cmin,upper=cmax)$root
}
qShewlogcalibratehitprob <- function(pobs, target, nsteps) log(qShewcalibratehitprob(pobs,target,nsteps))
setMethod("SPCq", signature="calhitprobShew",
          function(property,P,xi)
          qShewlogcalibratehitprob(pobs=getcdfupdates(property@chart, xi=xi, P=P),target=property@target,nsteps=property@nsteps)
          )





#################
#Computational routines
SPCproperty <- function(data, nrep=500,property, covprob=0.9){
    Phat <-  Pofdata(property@chart,data)
    qdiff <- replicate(nrep,{
        Phatstar<-Pofdata(property@chart,resample(property@chart,Phat))
        xihatstar <- xiofP(property@chart,Phatstar)
        SPCq(property,xi=xihatstar,P=Phatstar)-SPCq(property,xi=xihatstar,P=Phat)
    })
    raw <- SPCq(property,xi=xiofP(property@chart,Phat),P=Phat)
    if (property@lowerconf)
        res <- qtrafo(property,raw+quantile(-qdiff,covprob))
    else
        res <- qtrafo(property,raw-quantile(qdiff,covprob))
    new("SPCpropertyres", res=res, raw=qtrafo(property,raw),covprob=covprob,property=property,nrep=nrep)
}

setClass("SPCpropertyres", representation=list(nrep="numeric", property="SPCproperty", covprob="numeric",res="numeric",raw="numeric"))
setMethod("show", "SPCpropertyres",  function(object){
    for (i in 1:length(object@covprob)){
        cat(paste(strwrap(exdent=2,paste(object@covprob[i]*100,"% CI: ",SPCoutput(object@property,object@res[i]))),collapse="\n"),"\n")
    }
    cat("Unadjusted result: ", format(object@raw,digits=4),"\n")
    cat("Based on ",object@nrep, "bootstrap repetitions.\n")
})


SPC2sidedconfint <- function(data, nrep=500,property, covprob=0.9){
    sort(SPCproperty(data, nrep=nrep, property=property,  covprob=c((1-covprob)/2,1-(1-covprob)/2))@res)
}
