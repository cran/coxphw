\name{plotfp}
\alias{plotfp}
\title{Compute and plot nonlinear effect estimated by fractional polynomials}
\description{
             This function enables visualizing a nonlinear effect estimated by the fractional polynomial functionality of coxphw. It plots the log
             relative hazard versus values of a continuous covariable.
             }
\usage{
plotfp(obj, plot.x = NA, ref = NA, plot = TRUE, treatment=NULL, ref.type="value", variable=NULL, xlab=NULL, ylab=NULL, ...)}
\arguments{
  \item{obj}{an output object of \code{coxphw}, which includes \code{fp()} terms}
  \item{plot.x}{The data values for the continuous variable (e. g., plot.x=30:70)}
  \item{ref}{A reference value. The log relative hazard at this value will be 0.}
  \item{plot}{If TRUE, requests a plot. Otherwise, only an output object with x values and log relative hazards will be computed.}
  \item{treatment}{variable which is in interaction with fp variable (use "")}
  \item{ref.type}{"value" for simple fp term (no interaction),
  "interaction.time" for interaction for treatment with fp(time),
  "interaction.treat" for interaction of treatment with fp(variable),
  or any value to specify the level of the treatment variable for which the fp of "variable" should be plotted}
  \item{variable}{name of fp variable (use "")}
  \item{xlab}{label for x-axis of plot, uses variable specified in "variable" as default}
  \item{ylab}{label for y-axis of plot, uses "log relative hazard" as default}
  \item{...}{further parameters, to be used for plots (e.g., scaling of axes)}
    }
\details{
 After having fit an fp model, this function can be used to depict the nonlinearity. It supports simple non-linear effects
 as well as interaction effects of continuous variables with binary covariates (see examples section).
}
\value{
 \item{obj}{A vector with log relative hazard values.}
}
\author{Georg Heinze and Meinhard Ploner}
\seealso{coxphw}
\examples{
#Example 1
set.seed(30091)
n<-300
x<-1:n
true.func<-function(x) 3*(x/100)^{2}-log(x/100)-3*x/100
x<-round(rnorm(x)*10+40,0)
#x<-round(runif(x)*60+10,0)
time<-rexp(n,rate=1)/exp(true.func(x))
event<-rep(1,n)
futime<-runif(n,0,309000)
event<-(time<futime)+0
time[event==0]<-futime[event==0]
my.data<-data.frame(x,time,event)

fit<-coxphw(data=my.data, Surv(time,event)~fp(x),alpha.fp=c(1,0.05,0.05), AHR=TRUE, maxit=200)
fit1<-coxphw(data=my.data, Surv(time,event)~fp(x),alpha.fp=c(1,0.05,0.05), PH=TRUE)

# estimated function
y<-plotfp(fit, plot.x=10:70, ref=40, variable="x", ref.type="value")
y1<-plotfp(fit1, plot.x=10:70, ref=40, variable="x", ref.type="value")
# true function
x<-10:70
lines(x, true.func(x)-true.func(40))
lines(x, y1, lty=2)
legend(x=40, y=(max(true.func(x)-true.func(40))+min(true.func(x)-true.func(40)))/2, lty=c(NA,2,1),pch=c(1,NA,NA),legend=c("AHR estimates","PH estimates", "true"))

#Example 2
n<-200
x<-1:n
true.func<-function(x) 2.5*log(x)-2
#x<-round(rnorm(x)*10+40,0)
x<-round(runif(x)*60+10,0)

time<-round(100000*rexp(n,rate=1)/exp(true.func(x)),1)
event<-rep(1,n)
my.data<-data.frame(x,time,event)

# select best fp(1) model, without testing in RA2 algorithm
fit<-coxphw(data=my.data, Surv(time,event)~fp(x),alpha.fp=c(1,1,0.05), fp.max=1, AHR=TRUE)

# estimated function
y<-plotfp(fit,plot.x=10:70, ref=40, ref.type="value", variable="x", type="p")
# true function
x<-10:70
lines(x, true.func(x)-true.func(40))
legend(x=40, y=(max(true.func(x)-true.func(40))+min(true.func(x)-true.func(40)))/2, lty=c(NA,1),pch=c(1,NA),legend=c("estimates","true"))


#interaction of continuous variable with binary "treatment" variable

n<-200
trt<-rbinom(n,1,0.5)
x<-1:n
true.func<-function(x) 2.5*log(x)-2
#x<-round(rnorm(x)*10+40,0)
x<-round(runif(x)*60+10,0)

  time<-100*rexp(n,rate=1)/exp(true.func(x)/4*trt-(true.func(x)/4)^2*(trt==0))
  event<-rep(1,n)
  my.data<-data.frame(x,trt,time,event)
  fit<-coxphw(data=my.data, Surv(time,event)~fp(x)*trt,alpha.fp=c(1,1,0.05), AHR=TRUE)

# plots the interaction of trt with x (the effect of trt dependent on the values of x):
y<-plotfp(fit,variable="x",treatment="trt",ref.type="interaction.treat",plot.x=10:70)

# plots the effect of x in subjects with trt=0:
y0<-plotfp(fit,variable="x",treatment="trt",ref.type=0,plot.x=10:70, ref=40)

# plots the effect of x in subjects with trt=1:
y1<-plotfp(fit,variable="x",treatment="trt",ref.type=1,plot.x=10:70, ref=40)

#interaction of binary "treatment" variable with function of time

  time<-100*rexp(n,rate=1)/exp((true.func(x)/10)^2/2000*trt+trt)
  event<-rep(1,n)
  my.data<-data.frame(x,trt,time,event)


fit<-coxphw(data=my.data, Surv(time,event)~fp(time)*trt+x, alpha.fp=c(1,1,0.05), fp.max=2, AHR=TRUE)
y<-plotfp(fit,variable="time", treatment="trt", plot.x=seq(1,100,1), ref.type="interaction.time", type="l")

fit2<-coxphw(data=my.data, Surv(time,event)~fp(time)*trt+x, alpha.fp=c(1,1,0.05), fp.max=2, PH=TRUE)
y2<-plotfp(fit2,variable="time", treatment="trt", plot.x=seq(1,100,1), ref.type="interaction.time", type="l")

plot(x=1:100, xlab="time", y=y, ylab="log relative hazard", pch=1, ylim=c(min(y,y2),max(y,y2)))
lines(x=1:100, y=y2, lty=1)
legend(x=60, y=max(y,y2), lty=c(NA,1), pch=c(1,NA), legend=c("AHR","PH"))
}
\keyword{survival}

