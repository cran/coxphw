\name{plotshape}

\alias{plotshape}

\title{ Compute and Plot the Relative or Log Relative Hazard Versus Values of a Continuous Covariable }

\description{ This function either visualizes a nonlinear or a time-dependent effect of a \code{coxphw} 
              object (possibly estimated by the fractional polynomial functionality of \code{coxphw}) 
              or obtains the effect estimates at specified values of this variables. It plots/prints 
              the relative or log relative hazard versus values of a continuous covariable. }

\usage{ 
plotshape(obj, plot.x = NA, ref = NA, plot = TRUE, variable = NULL,
          treatment = NULL, ref.type = "value", addci = TRUE, exp = FALSE, 
          xlab = NULL, ylab = NULL, pval = FALSE, digits = 4, \dots) 
}
        
\arguments{
  \item{obj}{ an output object of \code{coxphw}, which may include \code{fp()}-terms. }
  \item{plot.x}{ the data values for the continuous variable (e.g., \code{plot.x = 30:70}). }
  \item{ref}{ a reference value. The log relative hazard at this value will be 0. (e.g., 
              \code{ref= 50}).}
  \item{plot}{ if set to TRUE (default), requests a plot, otherwise, a matrix with relative or 
               log relative hazards at values specified in \code{plot.x} are printed. }
  \item{variable}{ name of the variable with the nonlinear effect (use ""). }
  \item{treatment}{ variable which is in interaction with \code{variable} (or \code{fp(variable)}) (use ""). }
  \item{ref.type}{ "value" for a simple nonlinear effect of \code{variable} without interaction, \cr
                   "interaction.time" for interaction of \code{treatment} with \code{time} (or \code{fp(time)}), \cr
                   "interaction.treat" for interaction of \code{treatment} with \code{fp(variable)}, or \cr
                   any value to specify the level of the \code{treatment} variable for which the fractional 
                   polynomial of \code{variable} should be plotted. }
  \item{addci}{ confidence intervalls are obtained. Default is TRUE. }
  \item{exp}{ if set to TRUE (default), the log relative hazard is given, otherwise the relative hazard
              is requested. }
  \item{xlab}{ label for x-axis of plot, uses variable specified in \code{variable} as default. }
  \item{ylab}{ label for y-axis of plot, uses "relative hazard" or "log relative hazard" as default. }
  \item{pval}{ if \code{plot = FALSE}, add Wald-test p-values to effect estimates at values of 
               \code{plot.x}. Default is set to FALSE. }
  \item{digits}{ if \code{plot = FALSE}, number of printed digits. Default is 4. }
  \item{\dots}{ further parameters, to be used for plots (e.g., scaling of axes). }
}
    
\details{ This function can be used to depict the estimated nonlinear or time-dependent 
          effect of an object of class \code{coxphw} (possibly using fractional polynomials 
          \code{fp()}). It supports simple nonlinear effects as well as interaction effects of 
          continuous variables with binary covariates (see examples section). 
          
%          \code{plot.x} must include only positive numbers greater than 0. 
}

\value{ A list with the follwing components:
 \item{xbeta}{ a matrix with estimates of (log) relative hazard. }
 \item{cilower}{ lower confidence limits. }
 \item{ciupper}{ upper confidence limits. }
 \item{p}{ p-value, only if \code{plot = FALSE} and \code{pval = TRUE}. }
 \item{alpha}{ the significance level = 1 - confidence level. }
 \item{plot.x}{ the data values of the continuous variable. }
 \item{exp}{ an indicator if log relative hazard or relative hazard was obtained. }
}

\author{Georg Heinze, Meinhard Ploner, Daniela Dunkler}

\references{ Royston P and Altman D (1994). Regression Using Fractional Polynomials of Continuous 
Covariates: Parsimonious Parametric Modelling. \emph{Applied Statistics} \bold{43}, 429-467.

Royston P and Sauerbrei W (2008). \emph{Multivariable Model-building. A pragmatic approach to regression 
analysis based on fractional polynomials for modelling continuous variables.} Wiley, Chichester, UK. }

\seealso{ \code{\link{coxphw}} }

\examples{
# Example 1
set.seed(30091)
n <- 300
x <- 1:n
true.func <- function(x) 3 * (x / 100)^{2} - log(x / 100) - 3 * x / 100
x <- round(rnorm(n = x) * 10 + 40, digits = 0)
time <- rexp(n = n, rate = 1) / exp(true.func(x))
event <- rep(x = 1, times = n)
futime <- runif(n = n, min = 0, max = 309000)
event <- (time < futime) + 0
time[event == 0] <- futime[event == 0]
my.data <- data.frame(x, time, event)

fitahr <- coxphw(Surv(time, event) ~ fp(x), alpha.fp = c(1, 0.05, 0.05), data = my.data,
                 template = "AHR", maxit = 200, print = TRUE)   

# estimated function
plot.x <- quantile(x, probs = 0.05):quantile(x, probs = 0.95)
yahr <- plotshape(fitahr, plot.x = plot.x, ref = median(x), variable = "x", 
                  ref.type = "value", addci = FALSE, plot = TRUE, type = "l")

# true function
lines(x = plot.x, true.func(plot.x) - true.func(median(plot.x)), lty = 2)

legend("topright", legend=c("AHR estimates", "true"), bty = "n", lty = 1:2, inset = 0.05)



# Example 2: select best fp(1) model, without testing in RA2 algorithm
set.seed(512364)
n <- 200
x <- 1:n
true.func <- function(x) 2.5 * log(x) - 2
#x <- round(rnorm(x) * 10 + 40, digits = 0)
x <- round(runif(x) * 60 + 10, digits = 0)

time <- round(100000 * rexp(n= n, rate = 1) / exp(true.func(x)), digits = 1)
event <- rep(x = 1, times = n)
my.data <- data.frame(x,time,event)

fit <- coxphw(Surv(time, event) ~ fp(x), alpha.fp = c(1, 1, 1), data = my.data, fp.max = 1, 
              template = "AHR", print = FALSE)

# estimated function
plot.x <- quantile(x, probs = 0.05):quantile(x, probs = 0.95)
y <- plotshape(fit, plot.x = plot.x, ref = 40, ref.type = "value", variable = "x", type = "p")

# true function
lines(x = plot.x, y = true.func(plot.x) - true.func(40), col = "blue")

legend("topleft", legend=c("estimates", "true"), col = c("black", "blue"), bty = "n", 
       inset = 0.1, lty = c(NA, 1), pch=c(1, NA))



# Example 3: interaction of continuous variable with binary "treatment" variable
set.seed(75315)
n <- 200
trt <- rbinom(n = n, size = 1, prob = 0.5)
x <- 1:n
true.func <- function(x) 2.5*log(x)-2
#x <- round(rnorm(n = x) * 10 + 40, digits = 0)
x <- round(runif(n = x) * 60 + 10, digits = 0)
time <- 100 * rexp(n = n, rate = 1) / exp(true.func(x) / 
        4 * trt - (true.func(x) / 4)^2 * (trt==0))
event <- rep(x = 1, times = n)
my.data <- data.frame(x, trt, time, event)

fit <- coxphw(Surv(time, event) ~ fp(x) * trt, data = my.data, alpha.fp=c(1, 1, 0.05), 
              template = "AHR", print = FALSE)

# plots the interaction of trt with x (the effect of trt dependent on the values of x)
plot.x <- quantile(x, probs = 0.05):quantile(x, probs = 0.95)
y <- plotshape(fit, variable = "x", treatment = "trt", ref.type = "interaction.treat", 
               plot.x = plot.x, main = "interaction of trt with x", type = "l")

# plots the effect of x in subjects with trt = 0
y0 <- plotshape(fit, variable = "x", treatment = "trt", ref.type = 0, plot.x = plot.x, 
                type = "l", ref = median(x), main = "effect of x in subjects with trt = 0")

# plots the effect of x in subjects with trt = 1
y1 <- plotshape(fit, variable = "x", treatment = "trt", ref.type = 1, plot.x = plot.x, 
                type = "l", ref = median(x), main = "effect of x in subjects with trt = 1")



# Example 4: interaction of binary "treatment" variable with function of time
# set.seed(23917)   
# time <- 100 * rexp(n = n, rate = 1) / exp((true.func(x) / 10)^2 / 2000 * trt + trt)
# event <- rep(x = 1, times = n)
# my.data <- data.frame(x, trt, time, event)
# plot.x <- seq(from = 1, to = 100, by = 1)
#
# fitahr <- coxphw(Surv(time, event) ~ fp(time) * trt + x, alpha.fp = c(1, 1, 0.05), 
#                  data = my.data, fp.max = 2, template = "AHR")
# yahr <- plotshape(fitahr, variable = "time", treatment = "trt", 
#                   ref.type = "interaction.time", plot.x = plot.x, plot = FALSE)

# fitph <- coxphw(Surv(time, event) ~ fp(time) * trt + x, alpha.fp = c(1, 1, 0.05), 
#                 data = my.data, fp.max = 2, template = "PH")
# yph <- plotshape(fitph, variable = "time", treatment = "trt", 
#                  ref.type = "interaction.time", plot.x = plot.x, plot = FALSE)
#
# plot(x = plot.x, y = yahr$xbeta, xlab = "time", ylab = "log relative hazard", pch = 1, 
#      ylim = c(min(yahr$xbeta, yph$xbeta), max(yahr$xbeta, yph$xbeta)))
# lines(x = plot.x, y = yph$xbeta, lty = 1)
# legend("bottomright", legend = c("AHR", "PH"), bty = "n", lty = c(NA, 1), pch = c(1, NA), 
#        inset = 0.05)
}

\keyword{survival}
