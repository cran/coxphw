".First.lib" <-
function(lib, pkg) 
  library.dynam("coxphw", pkg, lib)

"decomposeSurv" <-
function(
 formula, 
 data,
 sort=FALSE,
 offset=NULL
)
### decomposes complex survival formula
### with time-factor interactions
### 2006-10
{
	orig.formula <- formula

	## expand formula if needed:
	
	repeat {
		terms <- terms(formula, data=data)
		fac <- attr(terms, "factors")

		needed <- rownames(fac)[!rownames(fac) %in% colnames(fac)][-1]
                
		if(length(needed) == 0) break
                
		formula <- as.formula(paste(as.character(formula)[2], "~", 
                                            as.character(formula)[3], "+",
                                            paste(needed, sep="+")))
	}
        
	## construct 3-col response:
	resp <- model.extract(model.frame(formula, data = data), response)
	if(ncol(resp) == 2)
          resp <- cbind(start=rep(0, nrow(resp)), resp)
        
	## sortieren nach STOPzeit und -Cens
	if(sort) {
                sort <- order(resp[, 2],  -resp[, 3])
                data <- data[sort, , drop=FALSE]
		resp <- resp[sort, ]
	}
        
	mm <- model.matrix(formula, data = data) ## Model-Matrix

	mm1 <- mm[, -1, drop=FALSE]	# w/o intercept
	if (length(offset) != 0) {
    offset.values<-offset
  }
  else offset.values<-NA
	terms <- terms(formula, data=data)
	fac <- attr(terms, "factors")
	labels <- attr(terms, "term.labels")
	
	## splittes by special chars
	f <- function(str)
          for(chars in c("(", ")", ":", " ", ",", "*", "^"))
            str <- unlist(strsplit(str, split=chars, fixed=TRUE))
        
	rowSplit <- sapply(rownames(fac), f)	# splitted effects 
	stopName <- tail(rowSplit[[1]], 2)[1]	# name of stoptime
	rowInter <- unlist(lapply(rowSplit[-1], function(z) any(z == stopName)))
#  rowOffset <- unlist(lapply(rowSplit[-1], function(z) any(z == "offset")))
#  rowOffset

        
	fac <- fac[-1, , drop=FALSE]	# omit Surv
        
	colSplit <- sapply(colnames(fac), f)
	colInter <- unlist(lapply(colSplit, function(z) any(z == stopName)))
        
	nTimes <- colSums(fac[rowInter, , drop=FALSE])
	nFac   <- colSums(fac[!rowInter, , drop=FALSE])
        
	inters <- (nFac>0) & (nTimes>0)
	NTDE <- sum(inters)
        
	
	timedata <- matrix(0, nrow(data), 0)
	timeind <- c()
        
	## loop for (time x effect)
	for(i in which(inters)) {
		## search pure time:
		ind <- (colSums(fac[rowInter, i] != fac[rowInter, , drop=FALSE]) == 0) & (nFac==0)
		timedata <- cbind(timedata, mm1[, ind, drop=FALSE])
                
		## search pure effect:
		ind <- (colSums(fac[!rowInter, i] != fac[!rowInter, , drop=FALSE]) == 0) & (nTimes == 0)
		timeind <- c(timeind, which(ind[!colInter]))
	}
	mm1 <- mm1[, !colInter, drop=FALSE]
	
	covnames <- c(colnames(mm1),
                      paste(colnames(timedata), colnames(mm1)[timeind], sep=":")
                      )
        alter <- c(colnames(mm1),
                   paste(colnames(mm1)[timeind], colnames(timedata), sep=":")
                   )
        
	## indicator to identify the original formula:
	ind <- covnames %in% colnames(attr(terms(orig.formula, data=data), "factors")) |
               alter %in% colnames(attr(terms(orig.formula, data=data), "factors"))
	
	list(NTDE=NTDE, 			# number time dep. effects
             fac=fac,	# factor matrix ..
             resp=resp, 			# N x 3 - response matrix
             mm1=mm1,	# model matrix without time effects
             timedata=timedata, 	# matrix with time functions as columns
             timeind=timeind, 		# indicator of time-dependend effect
             covnames=covnames,	# names of covariates
             ind=ind,			# indicator if some terms of not in formula
             offset.values=offset.values # offset values
             )
}


`coxphw` <- function
(
 formula, # formula, may include time-by-covariate interactions
 data=parent.frame(),                     #

 breslow=NA, # righthand formula, if breslow weighted terms, e.g. ~ A + A:C 
 prentice=NA, # righthand formula, if prentice weighted terms, e.g. ~ A + A:C
 taroneware=NA, # righthand formula, if tarone-ware weighted terms, e.g. ~ A + A:C
 id=NULL,
 robust=FALSE,
 jack=FALSE,
 normalize=TRUE,
 scale.weights=1,
 offset=NULL,

 alpha=0.05,                            # confidence limit
 maxit=50,                              # max. iterations
 maxhs=5,                               # half steps
 epsilon=1e-6,                          #
 maxstep=2.5,                           #
 censcorr=FALSE,
 x=TRUE,
 ...                                    # N -> breslow; km -> prentice
)
### by MP und GH, 2007
### Surv-Objekt 'simple'  (e.g. Surv(time,event)) or
### Counting-Process-Style (e.g. Surv(start,stop,event))
### event: 1=dead/event, 0=censored
{
        ## evt. use of abbrev. parameter names
        l <- list(...)
        if("N" %in% names(l)) breslow <- l$N
        if("km" %in% names(l)) prentice <- l$km
        
        n <- nrow(data)
        pl<-FALSE
        robcov<-0+1*(robust==T)+2*(jack==T)
	
	obj <- decomposeSurv(formula, data, sort=FALSE, offset)
        obj$resp[, 1] <- obj$resp[, 1] + 0.00001
        obj$resp[, 2] <- obj$resp[, 2] + 0.00002
        NTDE <- obj$NTDE                # number time-dep effects
	mmm <- cbind(obj$mm1, obj$timedata)
  ind.offset <- sum((length(offset)!=0))      
	cov.name <- obj$covnames
	
        k <- ncol(obj$mm1)              # number covariates
        n <- nrow(data)
        
        if(is.null(id)) id<-1:n
        else {
         orig.id<-id
         unique.id<-unique(orig.id)
         new.id<-order(unique.id)
         id<-1:n
         for (i in 1:n) id[i]<-new.id[unique.id==orig.id[i]]
        }
        
	## ####### weights ################
        weights <- matrix(1, n, k+NTDE) # weight matrix
        w.raw<-weights[,1]
        w<-w.raw
        formulaAll <- as.formula(paste(as.character(formula)[2], "1", sep="~"))
        my.survfit<-getFromNamespace("survfit","survival") 
        fit <- my.survfit(Surv(obj$resp[,1],obj$resp[,2],obj$resp[,3]), data)
        event <- event1 <- obj$resp[, 3]
        event1[-1][diff(obj$resp[, 2]) == 0] <- 0  
        ties <- (diff(c(-999, obj$resp[, 2])) == 0)
  if(censcorr) {
   time.obs.km <- (-(1:max(id)))
   cens.obs.km <- ((1:max(id)))
   for(id.i in 1:max(id)){
    time.obs.km[id.i]<-max(obj$resp[id==id.i,2])
    cens.obs.km[id.i]<-1-max(obj$resp[id==id.i,3])
   }
   obskm<-my.survfit(Surv(time.obs.km,cens.obs.km))
   obskm$surv<-cbind(1,obskm$surv)
   obskm$time<-cbind(0,obskm$time)
   if(!is.call(breslow)&!is.call(prentice)&!is.call(taroneware)) {
     w.obskm<-w
     for(i in 1:n) {
      w.obskm[i] <- 1/min(obskm$surv[obskm$time<obj$resp[i,2]])
      w[i]<-w.raw[i]*w.obskm[i]
     }
    }
    if(normalize)    w <- w * sum(event1) / sum((w * event1))
    weights<-matrix(w,length(w),k+NTDE)
  }
	if(is.call(breslow)) {
		formulaB <- as.formula(paste(as.character(formula)[2], 
                                             as.character(breslow)[2], sep="~"))
		objB <- decomposeSurv(formulaB, data, sort=TRUE)
		cov.nameB <- objB$covnames[objB$ind]
    # n at risk with counting process style
    ## count n rows with time1<time2[i] and time2>=time2[i]
    w <- (n:1)   #fit$n.risk
    for (i in 1:n){
     if(event[i]==1){
       w[i] <- sum((obj$resp[,1] < obj$resp[i,2]) & (obj$resp[,2] >= obj$resp[i,2]))      
     }
     else w[i]<-0
    }
    w[!event] <- 0
    w[ties] <- w[c(ties[-1], FALSE)]
    w.raw<-w
    if(censcorr) {
     w.obskm<-w
     for(i in 1:n) {
      w.obskm[i] <- 1/min(obskm$surv[obskm$time<obj$resp[i,2]])
      w[i]<-w.raw[i]*w.obskm[i]
     }
    }
    if(normalize)    w <- w * sum(event1) / sum((w * event1))
    weights[, match(cov.nameB, cov.name)] <- w
	 }
   if(is.call(prentice)) {
     formulaP <- as.formula(paste(as.character(formula)[2],
                    as.character(prentice)[2], sep="~"))
     objP <- decomposeSurv(formulaP, data, sort=TRUE)
     cov.nameP <- objP$covnames[objP$ind]
     w.sorted <- c(1, fit$surv[-length(fit$surv)]) #S_t
     w<-obj$resp[,3]
     for(i in 1:n) {
      w[i]<-w[i]*w.sorted[fit$time==obj$resp[i,2]]
     }
     w[!event] <- 0
     w[ties] <- w[c(ties[-1], FALSE)]
    w.raw<-w
    if(censcorr) {
     w.obskm<-w
     for(i in 1:n) {
      w.obskm[i] <- 1/min(obskm$surv[obskm$time<obj$resp[i,2]])
      w[i]<-w.raw[i]*w.obskm[i]
     }
    }
    if(normalize)    w <- w * sum(event1) / sum((w * event1))
    weights[, match(cov.nameP, cov.name)] <- w
    }
   if(is.call(taroneware)) {
     formulaTW <- as.formula(paste(as.character(formula)[2],
                                              as.character(taroneware)[2], sep="~"))
     objTW <- decomposeSurv(formulaTW, data, sort=TRUE)
     cov.nameTW <- objTW$covnames[objTW$ind]
     # n at risk with counting process style
     w <- (n:1)   #fit$n.risk
     for (i in 1:n){
      if(event1[i]==1){
        w[i] <- sum((obj$resp[,1] < obj$resp[i,2]) & (obj$resp[,2] >= obj$resp[i,2]))      
      }
      else w[i]<-0
     }
     w <- w ^0.5

       w[!event] <- 0
       w[ties] <- w[c(ties[-1], FALSE)]
    w.raw<-w
    if(censcorr) {
     w.obskm<-w
     for(i in 1:n) {
      w.obskm[i] <- 1/min(obskm$surv[obskm$time<obj$resp[i,2]])
      w[i]<-w.raw[i]*w.obskm[i]
     }
    }
    if(normalize)    w <- w * sum(event1) / sum((w * event1))
    weights[, match(cov.nameTW, cov.name)] <- w
    }

        ## number of weighted variables
  NGV <- sum(weights[1, ] != 1 | weights[n, ] != 1)
        
	## standardisierung
	sd1 <- sd(obj$mm1)
	sd2 <- sd(obj$timedata)
	Z.sd <- c(sd1, sd2 * sd1[obj$timeind])
	weights<-weights*scale.weights
  mm1.orig <- obj$mm1      
	obj$mm1 <- scale(obj$mm1, FALSE, sd1)
	if (ind.offset) obj$mm1o <- cbind(obj$offset.values, obj$mm1)
	else obj$mm1o <- obj$mm1
	obj$timedata <- scale(obj$timedata, FALSE, sd2)
	mmm <- cbind(obj$mm1, obj$timedata)
    DFBETA <- matrix(0,max(id),k+NTDE)
    CARDS <- cbind(obj$mm1o, obj$resp, weights, obj$timedata, id)
    PARMS <- c(n, k, robcov, maxit, maxhs, maxstep, epsilon, 0, 0, 0, 0, 0, NGV, NTDE, max(id), ind.offset)
    IOARRAY <- rbind(rep(1, k+NTDE), matrix(0, 2+3*(k+NTDE), k+NTDE))
    if(NTDE>0)
          IOARRAY[4, (k+1):(k+NTDE)] <- obj$timeind
        storage.mode(CARDS) <- "double"
        storage.mode(PARMS) <- "double"
        storage.mode(IOARRAY) <- "double"
#        storage.mode(DFBETA) <- "double"
	
        ## --------------- call Fortran-Routine WEIGHTEDCOX ----------------------------------
        value <- .Fortran("weightedcox",
                          CARDS = CARDS,
                          outpar = PARMS,
                          outtab = IOARRAY, 
                          PACKAGE=coxphw)
        if(value$outpar[8])
          warning("Error in routine WEIGHTEDCOX; parms8 <> 0")
        outtab <- matrix(value$outtab, nrow=3+3*(k+NTDE))
        iter <- value$outpar[10]
        ## --------------- process output object <fit> of class coxphw -------
        coef.orig <- outtab[3,  ]
        coefs <- coef.orig / Z.sd
        cov.ls<-matrix(outtab[4:(k+3+NTDE), ], ncol=k+NTDE) / (Z.sd %*% t(Z.sd))
        cov.lw<-NULL
        cov.j<-NULL
        dfbeta.resid<-NULL
        if(robust==T) {
           cov.lw<-matrix(outtab[(k+NTDE+4):(3+2*(k+NTDE)), ], ncol=k+NTDE) / (Z.sd %*% t(Z.sd))
           dfbeta.resid <- value$CARDS[1:max(id),1:(k+NTDE)] / t(matrix(Z.sd,k+NTDE,max(id)))
        }
        if(jack==T) {
           cov.j<-matrix(outtab[(3+2*(k+NTDE)+1):(3+3*(k+NTDE)), ], ncol=k+NTDE) / (Z.sd %*% t(Z.sd))
           dfbeta.resid <- value$CARDS[1:max(id),1:(k+NTDE)] / t(matrix(Z.sd,k+NTDE,max(id)))
           }
        cov.method<-"Lin-Sasieni"
        covs<-cov.ls
        if (robust==T) {
         covs<-cov.lw
         cov.method<-"Lin-Wei"
        }
        if (jack==T) {
          covs<-cov.j
          cov.method<-"Jackknife"
        }
        dimnames(covs) <- list(cov.name, cov.name)
        vars <- diag(covs)
        names(coefs) <- cov.name
        if(value$outpar[10]>=maxit) {
          cat("No convergence attained in ", value$outpar[10], " iterations.\n", sep="")
          }
        if(!censcorr){
         w.obskm<-rep(1,n)
         }  
        w.matrix<-cbind(time=obj$resp[,2],w.raw,w.obskm,w)[obj$resp[,3]!=0,]
        orderw<-order(w.matrix[,1])
        w.matrix<-w.matrix[orderw,,drop=F]
        fit <- list(coefficients = coefs,  dfbeta.resid = dfbeta.resid,
                    alpha = alpha, var = covs, df = k+NTDE,
                    iter = iter,
                    method.ties = "breslow", n = n, ##terms = terms(formula),
                    y = obj$resp, formula = formula,  call = match.call(),
                    cov.j=cov.j, cov.lw=cov.lw, cov.ls=cov.ls, cov.method=cov.method, w.matrix=w.matrix)
#        if(all(weights - weights[, 1] == 0)) {
#          fit$loglik <- value$outpar[12:11]
#          fit$score <- value$outpar[7]
#          }
        
        fit$Wald <- t(coefs)%*%solve(fit$var)%*%coefs
        
        fit$means <- apply(mmm, 2, mean)
        fit$linear.predictors <- as.vector(scale(mmm, fit$means, scale=FALSE) %*% coefs)
        fit$method <- "Weighted Estimation"
        
        fit$method.ci <- "Wald"
        fit$ci.lower <- exp(coefs + qnorm(alpha/2) * vars^0.5)
        fit$ci.upper <- exp(coefs + qnorm(1 - alpha/2) * vars^0.5)
        fit$prob <- 1 - pchisq((coefs^2/vars), 1)
        if(x){
         fit$x <- mm1.orig
         }
        fit$offset.values <- obj$offset.values
        names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- cov.name
        attr(fit, "class") <- c("coxphw", "coxph")
        fit
}

"plotw" <- function
(
 object,    # object of class coobjectphw
 rank=F,
 log=F,
  ...            # dummy
 )
 {
   if(rank) {
    time<-order(object$w.matrix[,1])
    label<-"Ranked time"
    }
   else {
    time<-object$w.matrix[,1]
    label<-"Time"
    }
   if(log) {
    weights<-log(object$w.matrix[,2:4])
    wlabel<-"Log of weight"
   }
   else {
    weights<-object$w.matrix[,2:4]
    wlabel<-"Weight"
    }
   plot(time, weights[,1],type="l",lty=1,ylim=c(min(weights),max(weights)), xlab=label, ylab=wlabel)
   lines(time, weights[,2],lty=2)
   lines(time, weights[,3],lty=3)
   legend(min(time),0.95*max(weights),c("Raw weight","Censoring weight","Normalized total weight"),lty=1:3, lwd=1)
 }
 

"print.coxphw" <- function
(
  x,     # object of class coxphw
  ...            # dummy
 )
### overall test is Wald test
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
        cat("\nWald test=", x$Wald, " on ", x$df, "df, p=", 1 - pchisq(x$Wald, x$df), ", n=", x$n, "\n\n", sep = "")
        
        invisible(x)
}


"summary.coxphw" <- function
(
 object,              # object of class coxphf
 ...                  # dummy
 )
### MP and GH
### 2007-07
{
        print(object$call)
        cat("\nModel fitted by", object$method, "\n\n")
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
        
        cat("Wald test =", object$Wald, "on", object$df, " df, p =",
            1 - pchisq(object$Wald, object$df))
        cat("\n\nCovariance-Matrix:\n")
        print(object$var)
        
        invisible(object)
}
