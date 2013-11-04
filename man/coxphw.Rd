\name{coxphw}
\alias{coxphw}
\title{Weighted estimation in Cox regression}
\description{
             This package implements weighted estimation in Cox regression, as proposed by Schemper, Wakounig and Heinze (2009).
             Weighted estimation provides unbiased average hazard ratio estimates in case of non-proportional hazards.
             An experimental feature allows the estimation of nonlinear effects using fractional polynomials. This feature can also be used
             to estimate the interaction of a covariate with a nonlinear function of survival time.
}
\usage{
coxphw( formula=attr(data, "formula"), data=sys.parent(), breslow=NA,
 prentice=NA, taroneware=NA, id=NULL, robust=FALSE,
 jack=FALSE, normalize=TRUE, scale.weights=1, offset=NULL, 
 alpha=0.05, alpha.fp=c(0.2, 0.05, 0.05), fp.iter=10, fp.max=2, maxit=200,
 maxhs=5, xconv=1e-4, gconv=1e-4, maxstep=1, x=TRUE, y=TRUE,  
 censcorr=FALSE, trunc.weights=1, round.times.to=0.00001, add.constant=0,
 print=TRUE, sorted=FALSE, AHR=TRUE, ARE=FALSE, PH=FALSE, AHR.norobust=FALSE,
 pc=TRUE, pc.time=TRUE,\dots)     
}
\arguments{
  \item{formula}{a formula object, with the response on the left of the  operator, and the
    model terms on the right. The response must be a survival object as returned by the 'Surv' function.
    \code{formula} may include \code{fp()}-terms (see section on fractional polynomials).
    }
  \item{data}{a data.frame in which to interpret the variables named in the \code{formula} argument. }
  \item{prentice}{a righthand formula (e. g., ~A+B) with the terms which should be estimated using Prentice (survival function) weights,
   or use prentice=TRUE to apply these weights to all model terms}
  \item{breslow}{a righthand formula with the terms which should be estimated using Breslow (N at risk) weights, or TRUE}
  \item{taroneware}{a righthand formula with the terms which should be estimated using Tarone-Ware (square root of N at risk) weights, or TRUE. 
     Note that specifications in prentice overrule breslow, and breslow overrules taroneware.}
  \item{alpha.fp}{A vector of length 2 which specifies p-values to include an \code{fp()}-term: \code{alpha.fp[1]} defines the threshold p-value
    for keeping at least the linear term in the model, \code{alpha.fp[2]} defines the threshold p-value for the nonlinear term. Default is \code{c(0.20, 0.05, 0.05)}.}
  \item{fp.iter}{Number of iterations for big fp-loop, default=10. This is important if several \code{fp()}-terms 
    are included in the model as inclusion or exclusion of one term based on the threshold p-value 
    may increase or decrease the significance of others.}
  \item{fp.max}{Highest power of fp( ) terms (applies to all fp-terms). Select 1 or 2 (default=2).}
  \item{normalize}{if TRUE, weights are normalized such that their sum is equal to the number of events. May speed up or enable convergence if for some variables no weighting is used.}
  \item{alpha}{the significance level (1-\eqn{\alpha} = the confidence level), 0.05 as default.}
  \item{maxit}{maximum number of iterations (default value is 200)}
  \item{maxhs}{maximum number of step-halvings per iterations (default value is 5). 
     The increments of the parameter vector in one Newton-Rhaphson iteration step are halved, 
     unless the new likelihood is greater than the old one, maximally doing \code{maxhs} halvings.}
  \item{xconv}{specifies the maximum allowed change in standardized parameter estimates to declare convergence. Default value is 0.0001.}
  \item{gconv}{specifies the maximum allowed change in score function to declare convergence. Default value is 0.0001.}
  \item{maxstep}{specifies the maximum change of (standardized) parameter values allowed
    in one iteration. Default value is 1.}
  \item{x}{requests copying explanatory variables into output object}
  \item{y}{requests copying survival information into output object}
  \item{id}{a vector of patient identification numbers, must be integers starting from 1. These IDs are used for computing the 
    robust covariance matrix. If id=NA (the default) the program assumes that each line of the data set refers to a distinct individual.}
  \item{robust}{set to TRUE if the robust variance estimate should be computed. This robust estimate is computed from V=D'D where D is the matrix
   of dfbeta residuals which are computed by \eqn{D'D=A^{-1}U'UA^{-1}} with \eqn{A} denoting the weighted sum of contributions to the second 
   derivative of the log likelihood and \eqn{U} denoting the nxp matrix of score residuals.}
   \item{jack}{set to TRUE if the variance should be based on a complete jackknife. Each individual (as identified by distinct values of the
   variable specified in id) is left out in turn. The resulting matrix of dfbeta residuals D is then used to compute the variance matrix: V=D'D.}
   \item{offset}{specifies a variable which is included in the model but its parameter estimate is fixed at 1.}
   \item{scale.weights}{specifies a scaling factor for the weights}
   \item{trunc.weights}{specifies a quantile at which the weights are to be truncated. Can be used to increase the precision of the estimates,
         particularly if \code{censcorr=TRUE} is used. Default=1 (no truncation). Recommended value is 0.95 for mild truncation.}
   \item{censcorr}{If set to TRUE, the weights are multiplied by the inverse marginal Kaplan-Meier estimator with reverse status indicator (default=F).}
   \item{AHR}{If set to TRUE, applies a template for average hazard ratio estimation (Schemper, Wakounig and Heinze, 2009): it sets prentice=TRUE, censcorr=TRUE, robust=TRUE}
   \item{AHR.norobust}{If set to TRUE, applies a template for average hazard ratio estimation but with the Lin-Sasieni covariance: it sets prentice=TRUE, censcorr=TRUE, robust=FALSE}
   \item{ARE}{If set to TRUE, applies a template for average regression effect estimation (Xu and O'Quigley, 2000): it sets prentice=NA, breslow=NA, taroneware=NA, censcorr=TRUE, robust=TRUE}
   \item{PH}{If set to TRUE, applies a template for proportional hazards regression (similar to coxph): it sets prentice=NA, breslow=NA, taroneware=NA, censcorr=FALSE, robust=TRUE}
   \item{round.times.to}{rounds survival times to the nearest number that is a multiple of \code{round.times.to}. May improve numerical stability if very low survival times are used (mostly in simulations).}
   \item{add.constant}{this number will be added to stop times, and \code{add.constant}/2 will be added to start times}
   \item{pc}{if set to TRUE, will transform model matrix to speed up convergence during fp mode. default=TRUE.}
   \item{pc.time}{if set to TRUE, will transform the time dependent covariates (interactions of covariates with functions of time)
           to speed up convergence. default=TRUE.}
   \item{print}{requests echoing of intermediate results (particularly, for fp estimation)}
   \item{sorted}{if set to TRUE, the data set will not be sorted prior to passing it to FORTRAN. May speed up computations.}
   \item{\dots}{additional arguments, e.g., \code{km} can be used as synonym of \code{prentice}, or \code{N} for \code{breslow}}
}
\details{
       If Cox's proportional hazards regression model is used in the presence of non-proportional hazards, 
i.e., with underlying time-dependent hazard ratios of prognostic factors, the average relative risk 
for such a factor is under- or overestimated and testing power for the corresponding regression parameter 
is reduced. In such a situation weighted estimation of this parameter provides a parsimonious alternative 
to more elaborate modelling of time-dependent effects. Weighted estimation in Cox regression (WCR) extends 
the tests by Breslow and Prentice to a multi-covariate situation as does the Cox model to Mantel's logrank 
test. WCR can also be seen as a robust alternative to the standard Cox estimator, reducing the influence 
of outlying survival times on parameter estimates. Schemper (1992) first demonstrated the suitability of 
WCR for estimating average hazard ratios when hazards are non-proportional and Sasieni (1993) extensively 
investigated the favorable properties of WCR.

Weighted estimation assigns weights to risk sets, according to the survival function estimates (Prentice weights), the number of subjects at risk (Breslow weights), or according to their 
square roots (Tarone-Ware weights). These
weights are applied to the summands of the score function. The final estimate is the vector of parameter values which
equates the score function to 0. Since there is one score function corresponding to each parameter
of the model, weights may be applied to some but not necessarily to all parameters of a model.    

Breslow's tie-handling method is used by the program, other methods to handle ties are currently not available.

By default, the program estimates the covariance matrix using Lin (1991) and Sasieni (1993) sandwich estimate \eqn{A^{-1}BA^{-1}} with \eqn{-A} and \eqn{-B} denoting the sum of contriubtions 
to the second derivative of the log likelihood, weighted by \eqn{w(t_j)} and \eqn{w(t_j)^2}, respectively. This estimate is independent from
the scaling of the weights and reduces to the inverse of the information matrix in case of no weighting. 

As an experimental feature, the RA-2 algorithm of Royston and Altman (1994) (cf. also Royston and Sauerbrei, 2008) 
to select \emph{fractional polynomial} transformations for continuous variables
was implemented. This allows to obtain for flexible nonlinear estimation. 
The method can also be applied to survival time in a time by covariate interaction.  
For some examples see the documentation of the \code{plotfp} function.
    }
\value{
 \item{coefficients}{the parameter estimates}
 \item{alpha}{the significance level = 1 - confidence level}
 \item{var}{the estimated covariance matrix}
 \item{df}{the degrees of freedom}
 \item{loglik}{the null and maximimized (penalized) log likelihood}
 \item{method.ties}{the ties handling method}
 \item{iter}{the number of iterations needed to converge}
 \item{n}{the number of observations}
 \item{y}{the response}
 \item{formula}{the model formula}
 \item{means}{the means of the covariates}
 \item{linear.predictors}{the linear predictors}
 \item{method}{the estimation method (usually weighted estimation)}
 \item{method.ci}{the confidence interval estimation method (Profile Likelihood or Wald)}
 \item{ci.lower}{the lower confidence limits}
 \item{ci.upper}{the upper confidence limits}
 \item{prob}{the p-values}
 \item{call}{the function call}
 \item{cov.ls}{the covariance matrix computed by the Lin-Sasieni method (the default method)}
 \item{cov.lw}{the covariance matrix computed by the Lin-Wei method (robust covariance, only computed if robust==T)}
 \item{cov.j}{the covariance matrix computed by the jackknife method (only computed if jack==T)}
 \item{cov.method}{the method used to compute the (displayed) covariance matrix and the standard errors. This method
 is either "jack" if jack==T, or "Lin-Wei" if robust==T and jack==F, or "Lin-Sasieni" if robust==F and jack==F.}
 \item{w.matrix}{A matrix with 4 columns and rows according to the number of uncensored failure times. The first column contains the failure times,
  the remaining columns (labeled w.raw, w.obskm, and w) contain the raw weights, the weights according to the inverse of the Kaplan-Meier estimates 
  with reverse status indicator and the normalized product of both.}
 \item{dataline}{The first dataline of the input data set (needed for \code{plotfp} function)}
}
\references{
Lin D and Wei L (1989). The robust inference for the Cox proportional hazards model. \emph{Journal of the American Statistical Association} 84, 1074-1078.

Lin D (1991). Goodness-of-fit analysis for the Cox regression model based on a class of parameter estimators. \emph{Journal of the American Statistical Association} 86, 725-728.

Royston P and Altman D (1994). Regression Using Fractional Polynomials of Continuous Covariates: Parsimonious Parametric Modelling. \emph{Applied Statistics} 43, 429-467.

Royston P and Sauerbrei W (2008). \emph{Multivariable Model-building. A pragmatic approach to regression analysis based on fractional polynomials
  for modelling continuous variables.} Wiley, Chichester, UK.

Sasieni P (1993). Maximum Weighted Partial Likelihood Estimators for the Cox Model. \emph{Journal of the American Statistical Association} 88, 144-152.

Schemper M (1992).  Cox analysis of survival data with non-proportional hazard functions. \emph{The Statistician} 41, 455-465.

Schemper M, Wakounig S and Heinze G (2009). The estimation of average hazard ratios by weighted Cox regression. \emph{Statistics in Medicine} 28, 2473-2489.

Xu R and O'Quigley J (2000). Estimating average regression effect under non-proportional hazards. \emph{Biostatistics} 1, 423-439.
}
\author{Georg Heinze and Meinhard Ploner}
\note{
A SAS macro WCM with similar functionality (but without fp options) is offered for download at \url{http://www.meduniwien.ac.at/cemsiis/kb/}
}
\seealso{coxph}
\examples{
# gastric cancer data set
gastric <-
  structure(list(patnr = as.integer(c(46, 1, 2, 3, 4, 5, 47, 6,
                   7, 8, 9, 48, 10, 11, 49, 12, 13, 14, 50, 15, 16, 17, 18, 19,
                   20, 51, 21, 22, 52, 23, 53, 54, 55, 24, 25, 56, 57, 58, 59, 60,
                   61, 62, 63, 64, 26, 65, 27, 66, 28, 29, 67, 68, 69, 70, 30, 71,
                   31, 72, 32, 73, 33, 34, 74, 75, 76, 77, 78, 35, 79, 36, 80, 81,
                   82, 37, 38, 39, 83, 84, 40, 85, 41, 86, 87, 88, 42, 43, 44, 89,
                   90, 45)),
                 treat = as.integer(c(0, 1, 1, 1, 1, 1, 0, 1, 1, 1,
                   1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0,
                   0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0,
                   0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
                   1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1)),
                 time = as.integer(c(1,
                   17, 42, 44, 48, 60, 63, 72, 74, 95, 103, 105, 108, 122, 125,
                   144, 167, 170, 182, 183, 185, 193, 195, 197, 208, 216, 234, 235,
                   250, 254, 262, 301, 301, 307, 315, 342, 354, 356, 358, 380, 383,
                   383, 388, 394, 401, 408, 445, 460, 464, 484, 489, 499, 523, 524,
                   528, 535, 542, 562, 567, 569, 577, 580, 675, 676, 748, 778, 786,
                   795, 797, 855, 955, 968, 977, 1174, 1214, 1232, 1245, 1271, 1366,
                   1420, 1455, 1460, 1516, 1551, 1585, 1622, 1626, 1690, 1694, 1736
                   )),
                 status = as.integer(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                   0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))),
            .Names = c("patnr",
              "treat", "time", "status"), class = "data.frame",
            row.names = c("1",
              "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
              "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
              "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35",
              "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46",
              "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57",
              "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68",
              "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79",
              "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90"
              ))
  


# weighted estimation of average hazard ratio
  
fit2<-coxphw(data=gastric, Surv(time,status)~treat, AHR=TRUE)
# equivalent: fit2<-coxphw(data=gastric, Surv(time,status)~treat)

summary(fit2)
fit2$cov.lw     # robust covariance
fit2$cov.ls     # Lin-Sasieni covariance

# unweighted estimation, include interaction with years ('treat' must be included
# in formula!) 

gastric$years<-gastric$time/365.25
fit3<-coxphw(data=gastric, Surv(time,status)~treat+years:treat, PH=TRUE)
summary(fit3)

# select best fp(1) for interaction with time using RA2 algorithm (don't include
# 'treat' in formula if using fp(time)!!)
fit4<-coxphw(data=gastric, Surv(time,status)~fp(time):treat, 
             alpha.fp=c(1,0.05,0.05), fp.max=1, AHR=FALSE, robust=FALSE)
summary(fit4)


}
\keyword{survival}

