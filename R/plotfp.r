plotfp<-function (obj, plot.x = NA, ref = NA, plot = TRUE, treatment=NULL, ref.type="value", variable=NULL, xlab=NULL, ylab=NULL,
    ...)
{
    ### version 2010-04-29: NOT COMPATIBLE TO OLDER VERSIONS!!!
    
    ### obj:   coxphf.beta object
    ### plot.x: range of x axis
    ### ref:    reference value (only for ref.type=="value" mode)
    ### plot:   if function should be plotted
    ### treatment: variable which is in interaction with fp variable (use "")
    ### ref.type:  "value" for simple fp term (no interaction),
    ###            "interaction.time" for interaction for treatment with fp(time)
    ###            "interaction.treat" for interaction of treatment with fp(variable)
    ###            any value to specify the level of the treatment variable for which the fp of "variable" should be plotted
    ### variable:  name of fp variable (use "")

    data<-obj$dataline      # a dataline from original data set of analysis
    if (is.null(variable)) stop("No variable specified.\n")
    if (ref.type=="interaction.time" & is.null(treatment)) stop("No interaction variable for time specified.\n")
    if (is.na(ref))
        ref <- min(plot.x)
    if (any(plot.x <= 0) & ref.type=="interaction.time")
        warning(paste("model.frame() might cancel lines with ", variable, "<=0..."))
    n.x <- length(plot.x)
    data <- data[1,]
    data[,variable]<-NA
    if(ref.type=="interaction.time" | ref.type=="interaction.treat") data[,treatment]<-1
    else if(ref.type != "value") data[,treatment]<-ref.type
    i.x <- which(is.na(data))
    plotData <- data[rep(1, n.x), , drop = FALSE]
    plotData[, i.x] <- plot.x
    rownames(plotData) <- 1:n.x
    objP <- decomposeSurv(obj$formula, plotData, sort = FALSE,
        obj$offset, PTpreset = obj$PTcoefs)
    kk <- ncol(objP$mm1)
    objP$mm1 <- objP$mm1[, obj$ind[1:kk], drop = FALSE]
    objP$covnames <- objP$covnames[obj$ind]
    objP$timedata <- objP$timedata[, obj$ind[-(1:kk)], drop = FALSE]
    plot.y <- (cbind(objP$mm1, objP$timedata)) %*% obj$coefficients

    if (ref.type != "interaction.time" & ref.type != "interaction.treat"){
        plotDataRef <- plotData
        plotDataRef[, i.x] <- ref
        objPRef <- decomposeSurv(obj$formula, plotDataRef, sort = FALSE,
        obj$offset, PTpreset = obj$PTcoefs)
        objPRef$mm1 <- objPRef$mm1[, obj$ind[1:kk], drop = FALSE]
        objPRef$covnames <- objPRef$covnames[obj$ind]
        objPRef$timedata <- objPRef$timedata[, obj$ind[-(1:kk)], drop = FALSE]
        plot.yRef <- (cbind(objPRef$mm1, objPRef$timedata)) %*% obj$coefficients
        plot.y <- plot.y - plot.yRef
        }
   if (ref.type == "interaction.treat"){
        plotDataRef <- plotData
        plotDataRef[, treatment] <- 0
        objPRef <- decomposeSurv(obj$formula, plotDataRef, sort = FALSE,
        obj$offset, PTpreset = obj$PTcoefs)
        objPRef$mm1 <- objPRef$mm1[, obj$ind[1:kk], drop = FALSE]
        objPRef$covnames <- objPRef$covnames[obj$ind]
        objPRef$timedata <- objPRef$timedata[, obj$ind[-(1:kk)], drop = FALSE]
        plot.yRef <- (cbind(objPRef$mm1, objPRef$timedata)) %*% obj$coefficients
        plot.y <- plot.y - plot.yRef
        }
    if (plot) {
        if(is.null(xlab)) xlab<-variable
        if(is.null(ylab)) ylab<-"log relative hazard"
        plot(plot.x, plot.y, lty = 1, xlab=xlab, ylab=ylab, ...)
        abline(0, 0, lty = 2)
    }
    return(plot.y)
}