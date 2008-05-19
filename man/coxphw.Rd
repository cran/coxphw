\name{coxphw}
\alias{coxphw}
\title{Weighted estimation in Cox regression}
\description{
             This package implements weighted estimation in Cox regression, as proposed by Schemper (1992) and
             Sasieni (1993). Weighted estimation provides unbiased average hazard ratio estimates in case of
             non-proportional hazards (Schemper, Wakounig and Heinze (2008)).
}
\usage{coxphw(formula, 
 data=parent.frame(),
 breslow=NA, 
 prentice=NA, 
 taroneware=NA, 
 id=NULL,
 robust=FALSE,
 jack=FALSE,
 normalize=TRUE,
 scale.weights=1,
 offset=NULL,
 alpha=0.05,
 maxit=50,
 maxhs=5, 
 epsilon=1e-6, 
 maxstep=2.5, 
 censcorr=FALSE,
 x=TRUE,
 \dots)
}
\arguments{
  \item{formula}{a formula object, with the response on the left of the  operator, and the
    model terms on the right. The response must be a survival object as returned by the 'Surv' function.
    }
  \item{data}{a data.frame in which to interpret the variables named in the 'formula' argument. }
  \item{breslow}{a righthand formula with the terms which should be estimated using Breslow (N at risk) weights}
  \item{prentice}{a righthand formula with the terms which should be estimated using Prentice (survival function) weights}
  \item{taroneware}{a righthand formula with the terms which should be estimated using Tarone-Ware (square of N at risk) weights}
  \item{normalize}{if T, weights are normalized such that their sum is equal to the number of events. May speed up or enable convergence if for some variables no weighting is used.}
  \item{alpha}{the significance level (1-\eqn{\alpha} = the confidence level), 0.05 as default.}
  \item{maxit}{maximum number of iterations (default value is 50)}
  \item{maxhs}{maximum number of step-halvings per iterations (default value is 5). 
     The increments of the parameter vector in one Newton-Rhaphson iteration step are halved, 
     unless the new likelihood is greater than the old one, maximally doing \code{maxhs} halvings.}
  \item{epsilon}{specifies the maximum allowed change in penalized log likelihood to
    declare convergence. Default value is 0.0001.}
  \item{maxstep}{specifies the maximum change of (standardized) parameter values allowed
    in one iteration. Default value is 2.5.}
  \item{id}{a vector of patient identification numbers, must be integers starting from 1. These IDs are used for computing the 
    robust covariance matrix. If id=NA (the default) the program assumes that each line of the data set refers to a distinct individual.}
  \item{robust}{set to TRUE if the robust variance estimate should be computed. This robust estimate is computed from V=D'D where D is the matrix
   of dfbeta residuals which are computed by \eqn{D'D=A^{-1}U'UA^{-1}} with \eqn{A} denoting the weighted sum of contributions to the second 
   derivative of the log likelihood and \eqn{U} denoting the nxp matrix of score residuals.}
   \item{jack}{set to TRUE if the variance should be based on a complete jackknife. Each individual (as identified by distinct values of the
   variable specified in id) is left out in turn. The resulting matrix of dfbeta residuals D is then used to compute the variance matrix: V=D'D.}
   \item{offset}{specifies a variable which is included in the model but its parameter estimate is fixed at 1.}
   \item{scale.weights}{specifies a scaling factor for the weights}
   \item{censcorr}{If set to TRUE, the weights are multiplied by the inverse marginal Kaplan-Meier estimator with reverse status indicator (default=F).}
   \item{x}{If set to TRUE, adds the model matrix to the output object.}
   \item{\dots}{additional parameters}
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

Weighted estimation assigns weights to risk sets, according to the number of subjects at risk (Breslow weights), according to their 
square roots (Tarone-Ware weights) or according to the survival function estimates (Prentice weights). These
weights are applied to the summands of the score function. The final estimate is the vector of parameter values which
equates the score function to 0. Since there is one score function corresponding to each parameter
of the model, weights may be applied to some but not necessarily to all parameters of a model.    

Breslow's tie-handling method is used by the program, other methods to handle ties are currently not available.

By default, the program estimates the covariance matrix using Lin (1991) and Sasieni (1993) sandwich estimate \eqn{A^{-1}BA^{-1}} with \eqn{-A} and \eqn{-B} 
denoting the sum of contributions to the second derivative of the log likelihood, weighted by \eqn{w(t_j)} and \eqn{w(t_j)^2}, respectively. 
This estimate is independent from the scaling of the weights and reduces to the inverse of the information matrix in case of no weighting. 
    }
\value{
 \item{coefficients}{the parameter estimates}
 \item{alpha}{the significance level = 1 - confidence level}
 \item{var}{the estimated covariance matrix}
 \item{df}{the degrees of freedom}
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
 \item{dfbeta}{the dfbeta residuals, only computed if robust=T, or if jack=T. For jack=T, the dfbeta residuals are the change in 
 parameter estimates if each individual in turn is left out. For robust=T, the dfbeta residuals are computed
  via \eqn{UA^{-1}} with \eqn{U} and \eqn{A} denoting the matrix of score residuals and minus the weighted sum of 
  contributions to the second
  derivative of the likelihood. If jack=T   or robust=T, the covariance matrix
  of the parameter estimates is computed from dfbeta'dfbeta.}
 \item{cov.ls}{the covariance matrix computed by the Lin-Sasieni method (the default method)}
 \item{cov.lw}{the covariance matrix computed by the Lin-Wei method (robust covariance, only computed if robust==T)}
 \item{cov.j}{the covariance matrix computed by the jackknife method (only computed if jack==T)}
 \item{cov.method}{the method used to compute the (displayed) covariance matrix and the standard errors. This method
 is either "jack" if jack==T, or "Lin-Wei" if robust==T and jack==F, or "Lin-Sasieni" if robust==F and jack==F.}
 \item{w.matrix}{A matrix with 4 columns and rows according to the number of uncensored failure times. The first column contains the failure times,
  the remaining columns (labeled w.raw, w.obskm, and w) contain the raw weights, the weights according to the inverse of the Kaplan-Meier estimates 
  with reverse status indicator and the normalized product of both.}
}
\references{
Lin D and Wei L (1989). The robust inference for the Cox proportional hazards model. \emph{Journal of the American Statistical Association} 84, 1074-1078.

Lin D (1991). Goodness-of-fit analysis for the Cox regression model based on a class of parameter estimators. \emph{Journal of the American Statistical Association} 86, 725-728.

Schemper M (1992).  Cox analysis of survival data with non-proportional hazard functions. \emph{The Statistician} 41, 455-465.

Sasieni P (1993). Maximum Weighted Partial Likelihood Estimators for the Cox Model. \emph{Journal of the American Statistical Association} 88, 144-152.

Schemper M, Wakounig S and Heinze G (2008). Weighted estimation for Cox regression revisited. \emph{submitted}.
}
\author{Georg Heinze and Meinhard Ploner}
\note{
A SAS macro WCM with similar functionality is offered for download at \url{http://www.meduniwien.ac.at/msi/biometrie}
}
\seealso{coxph}
\examples{
# gastric cancer data set
`gastric` <-
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
  


# weighted estimation, prentice (survival curve, Kaplan-Meier) weights
  
fit2<-coxphw(data=gastric, Surv(time,status)~treat, prentice= ~treat, jack=TRUE, censcorr=TRUE)
# equivalent: fit2<-coxphw(data=gastric, Surv(time,status)~treat, km= ~treat)

summary(fit2)
fit2$cov.lw
fit2$cov.ls
fit2$cov.j
plotw(fit2)

# weighted estimation, breslow (N at risk) weights                  

fit3<-coxphw(data=gastric, Surv(time,status)~treat, prentice= ~treat,  robust=TRUE, censcorr=TRUE)
# equivalent: fit2<-coxphw(data=gastric, Surv(time,status)~treat, km= ~treat)

summary(fit3)
fit3$cov.lw
fit3$cov.ls
fit3$cov.j


}
\keyword{survival}

