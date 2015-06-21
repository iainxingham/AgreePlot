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
#' @param xlab x axis label
#' @param ylab y axis label
#' @param main Plot title
#' @param ... Other parameters to be passed to plot()
#'
#' @references <bland altman - lancet?>
#'
#' @examples <to do> 
#' 
BAplot <- function(x, y, limit=1.96, 
                   xlab="Average of x & y",
                   ylab="Difference, x-y",
                   main="Bland Altman plot",
                   ...) {
  LA <- mean(x-y) + c(limit, 0-limit)*sd(x-y)  
  plot((x+y)/2, (x-y), type="p", xlab=xlab, ylab=ylab, 
       main=main, ...)
  abline(h=0, lty=2)
  abline(h=LA[1])
  abline(h=LA[2])
}
