".First.lib" <- function(lib, pkg) 
  library.dynam("coxphw", pkg, lib)


coxphw <- function
(
 formula=attr(data, "formula"), # formula, may contain time-interactions
 data=sys.parent(),                     #
 
 breslow=NA, # righthand formula, if breslow weighted terms, e.g. ~ A + A:C, or TRUE to weight all model effects
 prentice=NA, # righthand formula, if prentice weighted terms, e.g. ~ A + A:C, or TRUE to weight all model effects
 taroneware=NA, # righthand formula, if tarone-ware weighted terms, e.g. ~ A + A:C, or TRUE to weight all model effects
 id=NULL,                               # identifier: numeric or character or factor
 robust=FALSE,
 jack=FALSE,
 normalize=TRUE,
 scale.weights=1,
 offset=NULL,
 
 alpha=0.05,                            # confidence limit
  alpha.fp=c(0.20, 0.05, 0.05),          # alpha levels: [1] for fp(fp.max) vs. null, [2] for fp(fp.max) vs. linear, [3] for fp(2) vs. fp(1)
 fp.iter=10,                            # maximum number of iterations of large <fp> loop
 fp.max=2,                              # highest fp power (choose between 1 or 2)
 maxit=200,                              # max. iterations
 maxhs=5,                               # half steps
 xconv=1e-4,                            # convergence criterion for standardized parameter estimates
 gconv=1e-4,                            # convergence criterion for first derivatives of log likelihood
 maxstep=1,                             #
 x=TRUE,
 y=TRUE, 
 censcorr=FALSE,
 trunc.weights=1,
 round.times.to=0.00001,                 # whether times should be rounded (can be safer with FORTRAN interface), set to 0 for no rounding
 add.constant=0,                  # add 1*add.constant to start times, 2*add.constant to stop times
 
 print=TRUE,                            # print fitting information on screen
 sorted=FALSE,                          # if data is sorted by stoptime and -cens
 AHR=TRUE,                              # template to set prentice=TRUE, robust=TRUE, censcorr=TRUE
 ARE=FALSE,                             # template to set prentice=NA, breslow=NA, taroneware=NA, robust=TRUE, censcorr=TRUE
 PH=FALSE,                              # template to set prentice=NA, breslow=NA, taroneware=NA, robust=TRUE, censcorr=FALSE
 AHR.norobust=FALSE,                    #  template to set   prentice=TRUE, robust=FALSE, censcorr=TRUE
 pc=TRUE,                               # transforms data by principal components to speed up convergence when evaluating fp models
 pc.time=TRUE,                          # transforms time variables by principal components to speed up convergence in models
 ...                                    # N -> breslow; km -> prentice
 )
### by MP und GH, 2009
### Surv-Objekt entweder usual  (z.B. Surv(time,event)~A+B+C+D) oder
### Counting-Process-Style!! (z.B. Surv(start,stop,event)~A+B+C+D)
### event: 1=dead/event, 0=censored
{
        ## evt. use of abbrev. parameter names
        l <- list(...)
        if("N" %in% names(l)) breslow <- l$N
        if("km" %in% names(l)) prentice <- l$km
        
        
        # apply templates
        if (AHR.norobust==TRUE) {
          prentice<-TRUE
          censcorr<-TRUE
          robust<-FALSE
        } else if (PH==TRUE) {
          prentice<-NA
          breslow<-NA
          taroneware<-NA
          robust=TRUE
          censcorr<-FALSE
        } else if (ARE==TRUE) {
          prentice<-NA
          breslow<-NA
          taroneware<-NA
          robust=TRUE
          censcorr<-TRUE
        }   else  if (AHR==TRUE) {
          prentice<-TRUE
          censcorr<-TRUE
          robust<-TRUE
          }


        
        # preserve some hierarchy of weighting rules
        if (!is.na(prentice)){
         if (prentice==TRUE) {
            breslow<-NA
            taroneware<-NA
            }
         }
        if (!is.na(breslow)){
          if (breslow==TRUE) {
            taroneware<-NA
            }
          }
        
        n <- nrow(data)
        pl <- FALSE
        robcov <- if(robust) 1 else 0
        robcov <- if(jack) robcov+2 else robcov
        
        ## generate or reorder id's such that values are within 1:n
        if(is.null(id))
          id <- 1:n
        else
          id <- as.numeric(as.factor(id))
        maxid <- max(id)
        
        ## here only ONCE the full model matrix is spanned with all possible fp-effects
	      obj.full <- decomposeSurv(formula, data, sort=FALSE, offset)
        PTcoefs <- obj.full$PTcoefs
        if(round.times.to > 0) {
           obj.full$resp[, 1] <- round(obj.full$resp[, 1],-log10(round.times.to))
           obj.full$resp[, 2] <- round(obj.full$resp[, 2],-log10(round.times.to))
        }
        obj.full$resp[, 1] <- obj.full$resp[, 1] + add.constant         ### specify in option! round.times=automatic,or specify value
        obj.full$resp[, 2] <- obj.full$resp[, 2] + 2*add.constant
        
        ## calculate weights ...
        W <- coxphw.wei(formula, data, obj.full,
                        censcorr, normalize, breslow, prentice, taroneware,
                        scale.weights, trunc.weights, id, obj.full$covnames, obj.full$NFP)          ### include trunc.weights, id, NFP
        
        
        fit <- function(obj, pc=TRUE) {
                ind.offset <- sum(length(offset) != 0)
                kk <- ncol(obj.full$mm1) # should be k2 - NTDE
                obj$mm1 <- obj.full$mm1[, obj$ind[1:kk], drop=FALSE]
                if (pc==TRUE & sum(obj$ind[1:kk])>1) {
                   pc<-prcomp(obj$mm1)
                   obj$mm1<-predict(pc,obj$mm1)
                }
                obj$covnames <- obj.full$covnames[obj$ind]
                
                obj$timeind  <- obj.full$timeind[obj$ind[-(1:kk)]]
                obj$timedata <- obj.full$timedata[, obj$ind[-(1:kk)], drop=FALSE]
                ## re-index $timeind
                obj$timeind <- match(obj$timeind, (1:kk)[obj$ind[1:kk]])
                obj$NTDE <- length(obj$timeind)
                
                k <- ncol(obj$mm1) # number covariates w/o time-dep effects
                k2 <- k + obj$NTDE          #
                if(k2 == 0) return(0) # the null model !!??!!??! else crash ...
                
                WW <- W$weights[, obj$ind, drop=FALSE]
                NGV <- sum(WW[1, ] != 1 | WW[n, ] != 1)
                
                PARMS <- c(n, k, robcov, maxit, maxhs, maxstep, xconv, gconv, 0, 0, 0, 0, NGV, obj$NTDE, maxid, ind.offset)
                
                ## **************** fit model and return Wald *********************
                if(print) cat(paste(obj.full$covnames[obj$ind], collapse=", "), "\t\t")
                interim <- coxphw.fit(obj, id, WW, PARMS, sorted=sorted, pc=FALSE)
                wald<-interim$outpar[9] # FIT
                iterations<-interim$outpar[10]
                if (iterations>=maxit) warn<-" No convergence!"
                else warn<- paste(" OK (", iterations, " It.)")
                if(print) cat("Wald chi-square=", round(wald, 4), warn, " sum(abs(score))= ", interim$abs.score ,"\n")
                
                if(!is.finite(wald) | iterations>=maxit)
                  return(list(plr.statistic=0))
                # wald
                list(plr.statistic=wald)
        }
        
        collect<-NULL
        obj <- obj.full
        if(fp.max==2) test.stats<-matrix(0,0,6)
        else if(fp.max==1) test.stats<-matrix(0,0,5)
        else test.stats<-numeric(0)
        if(obj.full$NFP == 1) fp.iter<-1      # only 1 cycle needed for 1 fp term
        if(obj.full$NFP != 0) {
                if(print) cat("Entering fractional polynomial mode.\n")
                if(pl & print) cat("Profile penalized likelihood confidence intervals turned off.\n")
                pl<-FALSE

                obj$ind[] <- TRUE # !!??!!?? is it allowed? usually FALSE were terms like G when using G:log(time)

                ## start model: each fp() in model, but only linear
                obj$ind[apply(obj.full$fpind > 1, 2, any)] <- FALSE
                obj$ind[apply(obj.full$fpind == 1, 2, any)] <- TRUE

                ## ******* large loop until best model is found ('cycle') ********
                for(iter in 1:fp.iter) {
                        if(print) cat("\n** Iteration number", iter, "**\n")

                        ## remember inds for stop condition
                        ind.old <- obj$ind

                        ## ********* loop over the NFP fp()-variables **********
                        for(IFP in 1:obj$NFP) {
                                if(print) cat("\n** test fp() part number", IFP, "**\n")

                                ## ************* perform 'ra2'-algorithm **************
                                FPIND <- obj$fpind[IFP, ]
                                inds <- which(FPIND != 0) # where are the 8+8 MM columns
                                n.inds <- max(FPIND) # largest value (16)
                                if(n.inds != 16)
                                  stop("well, usually there should be 8 terms + 8 rep.powers ...")
                                nfp <- n.inds / 2 # nr terms w/o repeated powers

                                ## search best 2nd-degree model
                                if(fp.max==2){
                                 best2nd <- list(wald=0, i=NA, j=NA)
                                 for(iFP in 1:nfp)
                                   for(jFP in iFP:nfp) {
                                          ## for repeated power change 2nd index
                                          jstar <- if(iFP==jFP) jFP+nfp else jFP

                                          obj$ind[inds] <- FALSE # reset all indicators
                                          SEL <- which(FPIND %in% c(iFP, jstar))
                                          obj$ind[SEL] <- TRUE
                                          this.fit<-fit(obj,pc)
                                          wald <- this.fit$plr.statistic
#                                          collect<-rbind(collect,this.fit$line.collect)
                                          if(wald > best2nd$wald)
                                            best2nd <- list(wald=wald, i=iFP, j=jstar)
                                  }
                                 if(print) cat("Chi-Square of best fp(2)-model:", best2nd$wald, "\n")

                                ## null model (variable not in model)
                                 obj$ind[inds] <- FALSE # reset all indicators
                                 if(sum(obj$ind)!=0){
                                   this.fit <- fit(obj,pc)
                                   wald0 <- this.fit$plr.statistic
#                                   collect<-rbind(collect,this.fit$line.collect)
                                  }
                                  else wald0<-0

                                 ccond<-(1 - pchisq(best2nd$wald - wald0, df=4) > alpha.fp[1])
                                 if(ccond) test.stats<-rbind(test.stats, c(fp.iter, IFP, wald0, NA, NA, best2nd$wald))
                                 if(ccond) next # obj$ind is just right

                                 ## linear term
                                 SEL <- which(FPIND == 1)
                                 obj$ind[SEL] <- TRUE # set linear (1st) term
                                 this.fit<- fit(obj,pc)
                                 waldL <- this.fit$plr.statistic
#                                 collect<-rbind(collect,this.fit$line.collect)
                                 best1st <- list(wald=waldL, i=1)
                                 ccond<-((1 - pchisq(best2nd$wald - waldL, df=3) > alpha.fp[2])|(waldL>best2nd$wald))
                                 if(ccond) test.stats<-rbind(test.stats, c(fp.iter, IFP, wald0, waldL, NA, best2nd$wald))
                                 if(ccond) next # obj$ind is just right
 #                                  }
                                }
                                ## search best 1st-degree model
                                if(fp.max==1) best1st <- list(wald=0, i=NA)
                                for(iFP in fp.max:nfp) {
                                        obj$ind[inds] <- FALSE # reset all indicators
                                        SEL <- which(FPIND == iFP)
                                        obj$ind[SEL] <- TRUE # set the chosen
                                 this.fit<-fit(obj,pc)
                                 wald<-this.fit$plr.statistic
#                                 collect<-rbind(collect,this.fit$line.collect)

                                 if(wald > best1st$wald)
                                          best1st <- list(wald=wald, i=iFP)
                                }
                                if(print) cat("Chi-Square of best fp(1)-model:", best1st$wald, "\n")

                                if(fp.max==2){
                                 if((1 - pchisq(best2nd$wald - best1st$wald, df=2) > alpha.fp[3])|(best1st$wald>best2nd$wald))  {
                                        ## choose the 1st degree model as the best
                                        obj$ind[inds] <- FALSE # reset all indicators
                                        # changed next two lines 100121
#                                        obj$ind[inds[best1st$i]] <- TRUE # set the chosen
                                        obj$ind[FPIND %in% best1st$i] <- TRUE # set the chosen
                                 } else {
                                        ## choose the 2nd degree model as the best
                                        obj$ind[inds] <- FALSE # reset all indicators
                                        # changed next two lines 100121
#                                        obj$ind[inds[c(best2nd$i, best2nd$j)]] <- TRUE # set the 2 chosen
                                        obj$ind[FPIND %in% c(best2nd$i, best2nd$j)] <- TRUE # set the 2 chosen
                                 }
                                } else {
                                 obj$ind[inds] <- FALSE # reset all indicators
                                 if(sum(obj$ind)!=0){
                                  this.fit<-fit(obj,pc)
                                  wald0 <- this.fit$plr.statistic
#                                  collect<-rbind(collect,this.fit$line.collect)
                                 }
                                 else wald0<-0

                                 ccond<-(1-pchisq(best1st$wald - wald0, df=2) > alpha.fp[1])
                                 if(ccond) test.stats<-rbind(test.stats, c(fp.iter, IFP, wald0, waldL, best1st$wald, NA))
                                 if(ccond) next       # go for null model
                                 SEL <- which(FPIND == 1)
                                 obj$ind[SEL] <- TRUE # set linear (1st) term
                                 this.fit<-fit(obj,pc)
                                 waldL <- this.fit$plr.statistic
#                                 collect<-rbind(collect,this.fit$line.collect)
                                 ccond<-(1 - pchisq(best1st$wald - waldL, df=1) > alpha.fp[2])
                                 if(ccond) test.stats[nrow(test.stats),5]<-best1st$wald
                                 if(ccond) next# go for linear term
                                 #otherwise, go for fp(1) model
                                 obj$ind[inds] <- FALSE   # reset all indicators
                                 ## new 100122:
                                 SEL <- which(FPIND %in% best1st$i)
                                 obj$ind[SEL] <- TRUE
                                 # turned off 100122:
#                                 obj$ind[inds[best1st$i]] <- TRUE  # set the indicator of the best fp(1) power
                                }
                          if(fp.max==1) best2nd<-list(NA,NA)      
                          test.stats<-rbind(test.stats, c(fp.iter, IFP, wald0, waldL, best1st$wald, best2nd$wald))
                        }

                        ## check stop condition (nothing is changed)
                        if(all(ind.old == obj$ind))
                          break
                }
        if(print){
         cat("\nPretransformations (z+shift)/scale:\n")
         print(obj.full$PTcoefs)
         cat("\n")
         cat("Algorithm 'RA2' has selected:", obj.full$covnames[obj$ind], "\n\n")
         }
        if(!any(obj$ind)){
          stop("Null model is returned!")
        }                                  # change 091002: moved this line from line 190 to line 197
        }

        
        kk <- ncol(obj.full$mm1) # should be k2 - NTDE
        obj$mm1 <- obj.full$mm1[, obj$ind[1:kk], drop=FALSE]
        obj$covnames <- obj.full$covnames[obj$ind]
        
        obj$timeind  <- obj.full$timeind[obj$ind[-(1:kk)]]
        obj$timedata <- obj.full$timedata[, obj$ind[-(1:kk)], drop=FALSE]
        ## re-index $timeind
        obj$timeind <- match(obj$timeind, (1:kk)[obj$ind[1:kk]])
        NTDE <- obj$NTDE <- length(obj$timeind)
        
        W <- coxphw.wei(formula, data, obj,
                        censcorr, normalize, breslow, prentice, taroneware,
                        scale.weights, trunc.weights, id, obj$covnames, obj.full$NFP)
        
        ind.offset <- sum(length(offset) != 0)      
        
        k <- ncol(obj$mm1)    # number covariates w/o time-dep effects
        k2 <- k + NTDE                  #
        n <- nrow(data)
        
        ## generate or reorder id's such that values are within 1:n
        if(is.null(id))
          id <- 1:n
        else
          id <- as.numeric(as.factor(id))
        maxid <- max(id)
        
        PARMS <- c(n, k, robcov, maxit, maxhs, maxstep, xconv, gconv, 0, 0, 0, 0, W$NGV, NTDE, maxid, ind.offset)
        ##   if (offset) {
        ##    IOARRAY[1,1]<-0    # first variable is offset
        ##    IOARRAY[2,1]<-Z.sd[1]    # first variable is offset
        ##   }
        ## **************** fit model *********************
        value0 <- coxphw.fit(obj, id, W$weights, PARMS, sorted=sorted, pc=pc, pc.time=pc.time)
        cov.ls <- value0$cov.ls
        if(!robust)  { 
          cov.lw<-NULL
          dfbeta.resid<-NULL
        }
        if(!jack) {
          cov.j<-NULL
          if(!robust) dfbeta.resid<-NULL
        }
#        cov.lw <- cov.j <- dfbeta.resid <- NULL
        cov.method <- "Lin-Sasieni"
        covs <- cov.ls
        
        if(robust) {
                cov.lw <- value0$cov.lw
                if(!jack) dfbeta.resid <- value0$dfbeta.resid 
                covs <- cov.lw
                cov.method <- "Lin-Wei"
        } 
        if(jack) { 
                cov.j <- value0$cov.j
               dfbeta.resid <- value0$dfbeta.resid 
                covs <- cov.j
                cov.method <- "Jackknife"
        } 
        dimnames(covs) <- list(obj$covnames, obj$covnames)
        vars <- diag(covs)

        
        ## DECIDE THE MODEL: USE PVALS ##############
        probs <- 1 - pchisq((value0$coefs^2/vars), 1)
        
        
        
        ## ########## NOW ONLY FINAL MODEL IS CONSIDERED ############
        names(value0$coefs) <- obj$covnames
        if(value0$outpar[10]>=maxit)
          cat("No convergence attained in ", value0$outpar[10], " iterations.\n", sep="")
        
        Means <- colMeans(value0$mmm)
        
        ## return object
        fit <- list(coefficients = value0$coefs, # coefficients of the fit
                    cards    = value0$cards, #
                    parms    = value0$outpar,
                    ioarray  = value0$outtab,
                    dfbeta.resid = dfbeta.resid,
                    alpha    = alpha,   # significance level
                    var      = covs,    # covariance matrix
                    df       = k2,     # degrees of freedom (k + NTDE)
                    iter     = value0$outpar[10],
                    method.ties = "breslow", #
                    n = n,              # number observations
                    ##terms = terms(formula),
                    y = obj$resp,       # responses
                    formula = formula,  # original formula
                    exit.code=value0$outpar[8],
                    
                    fpind=obj$fpind,
                    PTcoefs=obj.full$PTcoefs,
                    ind=obj$ind,
                    
                    call    = match.call(),
                    offset=offset,
                    cov.j   = cov.j,    # jackknife covariance matrix
                    cov.lw  = cov.lw,   #
                    cov.ls  = cov.ls,   #
                    cov.method=cov.method,
                    w.matrix= W$w.matrix, # weight matrix
                    Wald    = value0$outpar[9],
                    means   = Means,    # means of <X> (model matrix)
                    linear.predictors= as.vector(scale(value0$mmm, Means, scale=FALSE) %*% value0$coefs),
                    method  = "Weighted Estimation",
                    method.ci= "Wald",  #
                    ci.lower= exp(value0$coefs + qnorm(alpha/2) * vars^0.5), #
                    ci.upper= exp(value0$coefs + qnorm(1 - alpha/2) * vars^0.5), #
                    prob    = 1 - pchisq((value0$coefs^2/vars), 1), # p-values
                    offset.values= obj$offset.values, # 
                    dataline = data[1,],
                    x       = if(x) obj$mm1 else NA # return original model matrix if requested
                    )
        ## if all weights are the same ...
        #if(W$const) {
        #        fit$loglik <- value0$outpar[12:11]  # no longer supported by fortran program!
        #        fit$score <- value0$outpar[7]       # no longer supported by fortran program!
        #}
        
        ##        if(offset) {
        ##          obj$covnames[1]<-as.character(paste("offset(",obj$covnames[1],")"))
        ##          fit$Wald <- t(coefs[2:k2]) %*% solve(covs[2:(k2),2:(k2)]) %*% coefs[2:(k2)]
        ##          if(x) {
        ##           fit$x <- mm1.orig[,2:k2]
        ##           fit$offset <- mm1.orig[,1]
        ##          }
        ##        }  
        names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- obj$covnames
        attr(fit, "class") <- c("coxphw.fp","coxphw", "coxph")
        fit
}

plotw <- function
(
 x,    # object of class coxphw
 rank=FALSE,
 log=FALSE,
  ...            # dummy
 )
 {
    w.matrix<-x$w.matrix[order(x$w.matrix[,1]),]
   if(rank) {
    time<-rank(w.matrix[,1])
    label<-"Ranked time"
    }
   else {
    time<-w.matrix[,1]
    label<-"Time"
    }
   if(log) {
    weights<-log(w.matrix[,2:4])
    wlabel<-"Log of weight"
   }
   else {
    weights<-w.matrix[,2:4]
    wlabel<-"Weight"
    }
   plot(time, weights[,1],type="l",lty=1,ylim=c(min(weights),max(weights)), xlab=label, ylab=wlabel)
   lines(time, weights[,2],lty=2)
   lines(time, weights[,3],lty=3)
   legend(min(time),0.95*max(weights),c("Raw weight","Censoring weight","Normalized total weight"),lty=1:3, lwd=1)
 }
 

print.coxphw <- function
(
  x,     # object of class coxphw
  ...            # dummy
 )
### MP and GH
### overall test is score test
### 2007-07
{
        print(x$call)
        cat("Model fitted by", x$method, "\n\n")
        se<- diag(x$var)^0.5
        out <- cbind(x$coefficients, se, exp(x$coefficients),
                     x$ci.lower, x$ci.upper, x$coefficients/se, x$prob)
        dimnames(out) <- list(names(x$coefficients),
                              c("coef", "se(coef)", "exp(coef)",
                                paste(c("lower", "upper"), 1 - x$alpha), "z", "p"))
        
        if (x$method.ci != "Wald")
          dimnames(out)[[2]][6] <- "Chisq"
        print(out)
        
#        if("loglik" %in% names(x)) cat("\nScore test=", x$score, " on ", x$df,
#            " df, p=", 1 - pchisq(x$score, x$df), ", n=", x$n, "", sep = "")
        cat("\nWald Chi-square=", x$Wald, " on ", x$df, "df, p=", 1 - pchisq(x$Wald, x$df), ", n=", x$n, "\n\n", sep = "")
        
        invisible(x)
}


summary.coxphw <- function
(
 object,              # object of class coxphf
 ...                  # dummy
 )
### MP and GH
### 2007-07
{
        print(object$call)
        cat("\nModel fitted by", object$method, "\n\n")
        ##cat("Confidence intervals and p-values by", object$method.ci, "\n\n")
        se<-diag(object$var)^0.5
        out <- cbind(object$coefficients, se,
                     exp(object$coefficients),
                     object$ci.lower, object$ci.upper,
                     object$coefficients/se, object$prob)
        dimnames(out) <- list(names(object$coefficients),
                              c("coef", "se(coef)", "exp(coef)",
                                paste(c("lower", "upper"), 1 - object$alpha), "z", "p"))
        if (object$method.ci != "Wald")
          dimnames(out)[[2]][6] <- "Chisq"
        
        print(out)
        
        if("loglik" %in% names(object)) {
                LL <- 2 * diff(object$loglik)
                cat("Likelihood ratio test=", LL, " on ", object$df,
                    " df, p=", 1 - pchisq(LL, object$df), ", n=", object$n, "\n", sep = "")
                cat("\nScore test=", object$score, " on ", object$df,
            " df, p=", 1 - pchisq(object$score, object$df), ", n=", object$n, "\n", sep = "")
        }
#        wald.z <- t(coef(object)) %*% solve(object$var) %*% coef(object)
        cat("Wald Chi-square =", object$Wald, "on", object$df, " df, p =",
            1 - pchisq(object$Wald, object$df))
        cat("\n\nCovariance-Matrix:\n")
        print(object$var)
        
        invisible(object)
}
