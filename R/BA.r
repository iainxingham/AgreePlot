## Functions for agreeplot package
## IRC May 2015

#' Bland Altman plot
#' 
#' \code{BAplot} draws a Bland Altman plot given a set of paired measurements
#' 
#' <Need to fill in details here>
#' 
#' @param x,y Two numeric vectors containing paired measurements
#' @param limit Number of standard deviations each side of \code{mean(x-y)}
#'   to draw the limits of agreement. Defaults to 1.96
#' @param subject A vector identifying the subject for replicates. NULL if only 
#'   one pair of measures per subject
#' @param measure.method A vector identifying the method used for x. Ignored 
#'   if y defined
#' @param xlab x axis label
#' @param ylab y axis label
#' @param main Plot title
#' @param ... Other parameters to be passed to plot()
#'
#' @references Bland JM, Altman DG. Statistical Methods in Medical Research 1999; 8: 135-160
#'
#' @examples <to do> 
#' 
BAplot <- function(x, y=NULL, limit=1.96,
                   subject=NULL,
                   measure.method=NULL,
                   xlab="Average of x & y",
                   ylab="Difference, x-y",
                   main="Bland Altman plot",
                   ...) {
  
  # If subject isn't defined then need pairs of measurements in x and y
  if(subject==NULL) {
    # Need to add if measure.method is defined
    if(y==NULL) stop("Need either paired x & y or for subject to be defined")
    if(length(x) != length(y)) stop("Need paired measurements in x and y")
    
    # Paired x and y, no replicates
    LA <- mean(x-y, na.rm=TRUE) + c(0-limit, limit)*sd(x-y, na.rm=TRUE)  
    plot((x+y)/2, (x-y), type="p", xlab=xlab, ylab=ylab, 
         main=main, ...)
    abline(h=0, lty=2)
    abline(h=LA[1])
    abline(h=LA[2])
  }
  else {
    # Vector in x with subject 
  }

}
