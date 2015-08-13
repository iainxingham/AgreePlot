## Functions for agreeplot package
## IRC May 2015

#' Bland Altman plot
#' 
#' \code{BAplot} draws a Bland Altman plot
#' 
#' <Need to fill in details here>
#' 
#' @param x,y Two numeric vectors containing paired measurements
#' @param limit Number of standard deviations each side of \code{mean(x-y)}
#'   to draw the limits of agreement. Defaults to 1.96
#' @param subject A vector identifying the subject for replicates. NULL if only 
#'   one pair of measures per subject
#' @param measure.method A factor vector identifying the method used for x. Ignored 
#'   if y defined
#' @param data An optional dataframe containing x, y, subject and/or 
#'   measure.method  
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
                   data=NULL,
                   xlab="Average of x & y",
                   ylab="Difference, x-y",
                   main="Bland Altman plot",
                   ...) {
  
  # If provided a dataframe then use this
  x <- eval(substitute(x), data, parent.frame())
  y <- eval(substitute(y), data, parent.frame())
  subject <- eval(substitute(subject), data, parent.frame())
  measure.method <- eval(substitute(measure.method), data, parent.frame())
  
  # If subject isn't defined then need pairs of measurements in x and y
  if(is.null(subject)) {
    if(is.null(y))
      stop("Need either paired x & y or for subject to be defined")
      
    # Paired x and y, no replicates    
    if(length(x) != length(y)) stop("Need paired measurements in x and y 
                                    if subject not defined")
    LA <- mean(x-y, na.rm=TRUE) + c(0-limit, limit)*sd(x-y, na.rm=TRUE)  
    plot((x+y)/2, (x-y), type="p", xlab=xlab, ylab=ylab, 
         main=main, ...)
    abline(h=0, lty=2)
    abline(h=LA[1])
    abline(h=LA[2])
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
      
      tempdf <- data.frame(x, y, subject, diff=(x-y), avg=(x+y)/2)
      anova_result <- anova(lm(diff~subject, data = tempdf))
      
      # Using dplyr for this - replace with aggregate()??
      #sumdf <- tempdf %>%
        #group_by(subject) %>%
        #summarise(avg=mean(diff), stddev=sd(diff), num=n())
      
      # Calculate mean differences to plot
      avgdiff <- aggregate(tempdf["diff"], by=tempdf["subject"], FUN=mean)
      avgread <- aggregate(tempdf["avg"], by=tempdf["subject"], FUN=mean)
      
      # Calculate sd for limits of agreement
      # <to do> Check handling of zeros here
      mi <- aggregate(tempdf["diff"], by=tempdf["subject"], FUN=length)
      divisor <- (sum(mi[2])^2 - sum(mi[2]^2)) / ((nrow(mi)-1) * sum(mi[2]))
      stdev <- sqrt(anova_result[2,3] + 
                      ((anova_result[1,3] - anova_result[2,3])/divisor))
      LA <- mean(avgdiff) + c(0-limit, limit) * stdev 
      
      plot(avgread[2], avgdiff[2], type="p", xlab=xlab, ylab=ylab, 
           main=main, ...)
      abline(h=0, lty=2)
      abline(h=LA[1])
      abline(h=LA[2])
    }
    
    # x only, method defined
    else {
      if(length(x) != length(measure.method)) 
        stop("x and measure.method should be equal length")
      if(length(x) != length(subject))
        stop("x and subject should be the same lenght")
      if(!is.factor(measure.method)) stop("measure.method should be a factor")
      if(nlevels(measure.method) < 2) stop("measure.method should have 
                                           at least two levels")
      
      # Could have bit here allowing user to specify levels from measure.method
      methlevels <- levels(measure.method)
      tempdf <- data.frame(x, subject, measure.method)
      # <to do>
      # Now cycle through subject and calculate 
      #   mean x per subject
      #   mean x per subject & measure.method == methlevels[1]
      #   mean x per subject & measure.method == methlevels[2]
    }  
  }

}
