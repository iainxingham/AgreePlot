## Functions for agreeplot package
## IRC May 2015

#' Bland Altman plot
#' 
#' \code{BAplot} draws a Bland Altman plot
#' 
#' <Need to fill in details here>
#' 
#' @param x,y Either two numeric vectors containing paired measurements made with
#'   different methods or a single numeric vector x with y=NULL and method defined.
#' @param limit Number of standard deviations each side of the mean difference
#'   to draw the limits of agreement. Defaults to 1.96.
#' @param subject A factor vector identifying the subject for replicates. 
#'   NULL if only one pair of measures per subject.
#' @param method A factor vector identifying the method used for x. Ignored 
#'   if y defined.
#' @param methodlevels An optional character vector specifying two levels of 
#'   the factor method to be included in the comparison   
#' @param data An optional dataframe containing x, y, subject and/or method
#' @param annotation Determines what annotation to add to the plot. Possible
#'   values include 
#'     \code{"zero"}    Plot a horizontal line at difference = zero
#'     \code{"bias"}    Plot a horizontal line at mean difference
#'     \code{"limits"}  Plot limits of agreement at +/- \code{limit} standard
#'                      deviations from the mean difference 
#' @param xlab x axis label
#' @param ylab y axis label
#' @param main Plot title
#' @param ... Other parameters to be passed to plot()
#'
#' @references Bland JM, Altman DG. Statistical Methods in Medical Research 1999; 8: 135-160
#'   Bland JM, Altman DG. Journal of Biopharmaceutical Statistics 2007; 17(4): 571-582
#'
#' @examples <to do> 
#' 
#' @export
BAplot <- function(x, y=NULL, limit=1.96,
                   subject=NULL,
                   method=NULL,
                   methodlevels=NULL,
                   data=NULL,
                   annotation=c("zero", "bias", "limits"),
                   xlab="Average of x & y",
                   ylab="Difference, x-y",
                   main="Bland Altman plot",
                   ...) {
  
  # If provided a dataframe then use this
  x <- eval(substitute(x), data, parent.frame())
  y <- eval(substitute(y), data, parent.frame())
  subject <- eval(substitute(subject), data, parent.frame())
  method <- eval(substitute(method), data, parent.frame())
  
  # If subject isn't defined then need pairs of measurements in x and y
  if(is.null(subject)) {
    if(is.null(y))
      stop("Need either paired x & y or for subject to be defined")
      
    # Paired x and y, no replicates    
    if(length(x) != length(y)) stop("Need paired measurements in x and y 
                                    if subject not defined")
    
    # Remove NAs
    nas <- is.na(x) | is.na(y)
    if(sum(nas) > 0) {
      warning("Excluding observations with missing data")
      x <- x[!nas]
      y <- y[!nas]
    }
    
    LA <- mean(x-y) + c(0-limit, limit)*sd(x-y)  
    
    plot(c(max((x+y)/2), min((x+y)/2)), c(max(LA[2], (x-y)), min(LA[1], (x-y))), 
         type="n", xlab=xlab, ylab=ylab, main=main, ...)
    points((x+y)/2, (x-y))
    
    if("bias" %in% annotation) abline(h=mean(x-y))
    if("zero" %in% annotation) abline(h=0, lty=2)
    if("limits" %in% annotation) {
      abline(h=LA[1])
      abline(h=LA[2])
    }
  }
    
  # Subject defined
  else {
    # Different methods in x and y - paired 
    if(!is.null(y)) {
      if(length(x) != length(y))
        stop("Should have paired measurements in x and y")
      if(length(x) != length(subject))
        stop("Subject and x vectors should be equal length")
      if(!is.factor(subject)) stop("Subject should be a factor")

      # Remove NAs
      nas <- is.na(x) | is.na(y) | is.na(subject)
      if(sum(nas) > 0) {
        warning("Excluding observations with missing data")
        x <- x[!nas]
        y <- y[!nas]
        subject <- subject[!nas]
      }
      
      tempdf <- data.frame(x, y, subject, diff=(x-y))
      anova_result <- anova(lm(diff~subject, data = tempdf))
      
      # Using dplyr for this - replace with aggregate()??
      #sumdf <- tempdf %>%
        #group_by(subject) %>%
        #summarise(avg=mean(diff), stddev=sd(diff), num=n())
      
      # Calculate sd for limits of agreement
      mi <- aggregate(tempdf["diff"], by=tempdf["subject"], FUN=length)
      divisor <- (sum(mi[,2])^2 - sum(mi[,2]^2)) / ((nrow(mi)-1) * sum(mi[,2]))
      stdev <- sqrt(anova_result[2,3] + 
                      ((anova_result[1,3] - anova_result[2,3])/divisor))
      LA <- mean(x-y) + c(0-limit, limit) * stdev 
      
      plot(c(max((x+y)/2), min((x+y)/2)), c(max(LA[2], (x-y)), min(LA[1], (x-y))), 
           type="n", xlab=xlab, ylab=ylab, main=main, ...)
      points((x+y)/2, x-y)
      
      if("bias" %in% annotation) abline(h=mean(x-y))
      if("zero" %in% annotation) abline(h=0, lty=2)
      if("limits" %in% annotation) {
        abline(h=LA[1])
        abline(h=LA[2])
      }
    }
    
    # x only, method defined
    else {
      if(length(x) != length(method)) 
        stop("x and method should be equal length")
      if(length(x) != length(subject))
        stop("x and subject should be the same length")
      if(!is.factor(method)) stop("method should be a factor")
      if(nlevels(method) < 2) stop("method should have at least two levels")
      
      # Remove NAs
      nas <- is.na(x) | is.na(subject) | is.na(method)
      if(sum(nas) > 0) {
        warning("Excluding observations with missing data")
        x <- x[!nas]
        subject <- subject[!nas]
        method <- method[!nas]
      }
      
      
      # Allow user to specify levels from method
      if(is.null(methodlevels)) methodlevels <- levels(method)[1:2]
      else {
        if(length(methodlevels) < 2) {
          warning("Invalid methodlevels - using first two levels of method instead")
          methodlevels <- levels(method)[1:2]
        }
        else if(length(methodlevels) > 2) {
          warning("methodlevels has more than two items - using first two")
          methodlevels <- methodlevels[1:2]
        }
        levelcheck <- methodlevels %in% levels(method)
        if(!(levelcheck[1] & levelcheck[2])) {
          warning("methodlevels contains items that are not levels of method.
                  Using first two levels of method instead")
          methodlevels <- levels(method)[1:2]
        }
      }
      
      tempdf <- data.frame(x, subject, method)

      tempx <- tempdf[tempdf$method==methodlevels[1],]
      tempy <- tempdf[tempdf$method==methodlevels[2],]
      mix <- aggregate(tempx["x"], by=tempx["subject"], FUN=length)
      miy <- aggregate(tempy["x"], by=tempy["subject"], FUN=length)

      # Check for subjects with no measurements by one method and remove
      noy <- mix$subject %in% miy$subject
      nox <- miy$subject %in% mix$subject
      if((!(sum(nox) == length(nox))) | (!(sum(noy) == length(noy)))) {
        warning("One or more subjects have measurements by one method only")
        tempx <- tempx[tempx$subject %in% mix$subject[noy],]
        tempy <- tempy[tempy$subject %in% miy$subject[nox],]
        mix <- aggregate(tempx["x"], by=tempx["subject"], FUN=length)
        miy <- aggregate(tempy["x"], by=tempy["subject"], FUN=length)
      }
      
      anovax <- anova(lm(x~subject, data=tempx))
      avgx <- aggregate(tempx["x"], by=tempx["subject"], FUN=mean)
      divx <- 1 - ((1/nrow(mix))*(sum(1/mix[,2])))
      
      anovay <- anova(lm(x~subject, data=tempy))
      avgy <- aggregate(tempy["x"], by=tempy["subject"], FUN=mean)
      divy <- 1 - ((1/nrow(miy))*(sum(1/miy[,2])))
      
      stdev <- sqrt(var(avgx[,2]-avgy[,2]) + divx*anovax[2,3] + divy*anovay[2,3])
      LA <- weighted.mean(avgx[,2]-avgy[,2], w=mix[,2]+miy[,2]) + c(0-limit, limit) * stdev
      
      plot(c(max((avgx[,2]+avgy[,2])/2), min((avgx[,2]+avgy[,2])/2)), 
           c(max(LA[2], (avgx[,2]-avgy[,2])), min(LA[1], (avgx[,2]-avgy[,2]))), 
           type="n", xlab=xlab, ylab=ylab, main=main, ...)
      points((avgx[,2]+avgy[,2])/2, (avgx[,2]-avgy[,2]))

      if("bias" %in% annotation) 
        abline(h=weighted.mean(avgx[,2]-avgy[,2], w=mix[,2]+miy[,2]))
      if("zero" %in% annotation) abline(h=0, lty=2)
      if("limits" %in% annotation) {
        abline(h=LA[1])
        abline(h=LA[2])
      }
      
      #<to do>
      # ?? Denote mi in averaged plot
    }  
  }
}
