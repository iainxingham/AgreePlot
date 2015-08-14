## Functions for agreeplot package
## IRC June 2015


### --- Internal functions ---

## Get quadrant
# x, y    Point co ordinates
# limit   Extent of exclusion zone
# centre  Centre of four quadrant plot
# zone    Use a exclusion zone?
#
# Returns 
#   0 Exclusion zone
#   1 Right upper
#   2 Right lower
#   3 Left lower
#   4 Left upper
get_quadrant <- function(x, y, limit, centre=0, zone=TRUE) {
  quad <- 1
  
  if(zone) {
    if((x > (centre-limit)) & (x < (centre+limit)) & 
         (y > (centre-limit)) & (y < (centre+limit)))
      quad <- 0 
  }
  else limit <- 0     # No zone so set limit to 0
  
  if(quad) {
    if(x < (centre-limit)) {
      if(y > centre) quad <- 4 
      else quad <- 3 
    }
    else if(x > centre+limit) {
      if(y > centre) quad <- 1 
      else quad <- 2   
    }
    else if(y > centre+limit) {
      if(x > centre) quad <- 1
      else quad <- 4     
    }
    else {
      if(x > centre) quad <- 2
      else quad <- 3    
    }
  }
  
  quad
}

## Concordance function for four quadrant plot
#   x1, y1  numeric vectors of paired measurements - first reading
#   x2, y2  numeric vectors of paired measurements - second reading
#   excl_zone     Create an exclusion zone?
#   excl_frac     Central fraction to exclude
#   exclude_mean  Fraction of what
concord_rate <- function(x1, y1, x2, y2, excl_zone=TRUE, excl_frac=0.15,
                         excl_mean=mean(x1+y1+x2+y2, na.rm=TRUE)/4) {
  
  if(excl_zone) excl_limit <- excl_mean * excl_frac
  else excl_limit <- 0
  
  delta_mon1 <- x2-x1
  delta_mon2 <- y2-y1
  quads <- sapply(1:length(delta_mon1), 
                  function(i, ...) 
                    get_quadrant(delta_mon1[i], delta_mon2[i], ...), 
                  limit=excl_limit, centre=0, zone=excl_zone)
  
  out_zone <- sum(quads>0)
  good_conc <- sum((quads==1) | (quads==3))
  # Concordance rate
  good_conc / out_zone
}

## Function to draw circles, for polarplot()
# x, y    centre of circle
# radius  circle radius
# nv      number of vertices for call to polygon
drawcircle <- function(x, y, radius, nv=100) {
  dtheta <- 2 * pi / nv
  thetas <- seq(0, (2*pi)-dtheta, by=dtheta)
  xv <- radius * cos(thetas) + x
  yv <- radius * sin(thetas) + y
  polygon(xv, yv)
}

### --- Functions to be visible ---

#' Four quadrant plot 
#'
#' \code{fourquadplot} draws a four quadrant plot from two sets of 
#' paired measurements
#' 
#' <details here>
#' 
#' @param x1,y1  numeric vectors of paired measurements - first reading
#' @param x2,y2  numeric vectors of paired measurements - second reading
#' @param data   Optional data frame containing x1, y1, x2, and/or y2
#' @param exclusion     Central fraction to exclude
#' @param exclude_mean  Fraction of what?
#' @param regress_line  Include summary or comparison lines? Takes a value of
#' \code{"none"} (draw no extra lines), \code{"x=y"} (the default - draw the
#' line of perfect concordance), \code{"lm"} (draw the linear regression line
#' y~x) or \code{"both"}
#' @param regress_lty,regress_col Line type and colour to draw regression line
#' @param xy_lty,xy_col  Line type and colour to draw the line x=y
#' @param zone_density,zone_col  Passed to polygon to fill exclusion zone
#' @param xlab,ylab,main,...  passed to plot()
#' 
#' @examples <to do>
#' 
#' @export
fourquadplot <- function(x1, y1, x2, y2, 
                         data=NULL,
                         exclusion=0.15,
                         exclude_mean=mean(x1+y1+x2+y2, na.rm=TRUE)/4,
                         regress_line="x=y",
                         regress_lty=1, regress_col="black",
                         xy_lty=2, xy_col="black",
                         zone_density=10, zone_col="red",
                         xlab="Change recorded with method x (x2-x1)", 
                         ylab="Change recorded with method y (y2-y1)", 
                         main="Four quadrant plot", ...) {
  
  # If provided a dataframe then use this
  x1 <- eval(substitute(x1), data, parent.frame())
  y1 <- eval(substitute(y1), data, parent.frame())
  x2 <- eval(substitute(x2), data, parent.frame())
  y2 <- eval(substitute(y2), data, parent.frame())
  
  x <- x2-x1
  y <- y2-y1
  plot(x, y, type="p", xlab=xlab, ylab=ylab, main=main, ...)
  
  abline(h=0)
  abline(v=0)
  
  excl_limit <- exclude_mean * exclusion
  polygon(c(-excl_limit, -excl_limit, excl_limit, excl_limit), 
          c(-excl_limit, excl_limit, excl_limit, -excl_limit),
          density=zone_density, col=zone_col)
  
  if(regress_line == "both") {
    abline(a=0, b=1, lty=xy_lty, col=xy_col)
    rl <- lm(y~x)
    abline(rl, lty=regress_lty, col=regress_col)
  }
  else if(regress_line == "x=y") {
    abline(a=0, b=1, lty=xy_lty, col=xy_col)
  }
  else if(regress_line == "lm") {
    rl <- lm(y~x)
    abline(rl, lty=regress_lty, col=regress_col)
  }
}

#' Polar plot
#' 
#' \code{polarplot} plots a polar plot of trend data
#' 
#' <details here>
#' 
#' @param x1,y1  numeric vectors of paired measurements - first reading
#' @param x2,y2  numeric vectors of paired measurements - second reading
#' @param axis_ticks  number of ticks to include on axis
#' @param limit_type  limit lines to be drawn. Takes a value of \code{"none"} 
#' (no limit lines drawn), \code{"absolute"} (limit lines are drawn at +/- 
#' limit_var) or \code{"fraction"} (limit lines are drawn +/- limit_var * mean
#' of all the measurements)
#' @param limit_var   parameter for drawing limit lines (See limit type above)
#' 
#' @references Critchley LA, Lee A & Ho A M-H. A Critical Review of the Ability
#' of Continuous Cardiac Output Monitors to Measure Trends in Cardiac Output.
#' Anesth Analg. 2010; 111: 1180-92
#' 
#' @examples <to do>
#'
#' @export 
polarplot <- function(x1, y1, x2, y2, axis_ticks=5, limit_type="fraction",
                      limit_var=0.1,
                      main = "Polar plot", ...) {
  deltax <- x2 - x1
  deltay <- y2 - y1
  hemis <- sapply(deltay, function(j) ifelse(j>0, 0, 1))
  meandelta <- abs((deltax-deltay)/2)
  theta <- atan(deltay/deltax)+((hemis-0.25)*pi)
  plotx <- cos(theta) * meandelta
  ploty <- sin(theta) * meandelta
  
  graphsize <- max(abs(meandelta))
  plot(c(-graphsize, graphsize), c(-graphsize, graphsize), type="n", 
       asp=1, xaxt="n", yaxt="n", main=main, xlab="", ylab="", ...)
  points(plotx, ploty)
  
  # As plot adds 4% to axis size, need to use this for calculations
  graphsize <- graphsize*1.04
  
  # Offset for text
  ## textoff <-
  
  # Outer circle and grid
  drawcircle(0, 0, graphsize)
  lines(c(0,0), c(-graphsize, graphsize), col="grey14")
  lines(c(-graphsize, graphsize), c(0,0), col="grey14")
  lines(c(-graphsize, graphsize)*sin(pi/4), 
        c(-graphsize, graphsize)*sin(pi/4), col="grey14")
  lines(c(-graphsize, graphsize)*sin(pi/4), 
        c(graphsize, -graphsize)*sin(pi/4), col="grey14")
  
  # Labels
  text(graphsize, 0, "0\xB0")
  text(graphsize*sin(pi/4), graphsize*sin(pi/4), "45\xB0")
  text(0, graphsize, "90\xB0")
  text(-graphsize*sin(pi/4), graphsize*sin(pi/4), "135\xB0")
  text(-graphsize, 0, "180\xB0")
  text(-graphsize*sin(pi/4), -graphsize*sin(pi/4), "225\xB0")
  text(0, -graphsize, "270\xB0")
  text(graphsize*sin(pi/4), -graphsize*sin(pi/4), "315\xB0")

  # Inner circles
  ticksize <- graphsize/(axis_ticks+1)
  for(i in 1:axis_ticks) {
    drawcircle(0, 0, i*ticksize)
    text(0, i*ticksize, signif(i*ticksize, 2))
  }
  
  # Limit lines
  if(limit_type == "fraction") {
    lineoff <- limit_var*mean(c(x1, y1, x2, y2))
  }
  else if(limit_type == "absolute") {
    lineoff <- limit_var
  }
  else lineoff <- 0
  
  if(lineoff) {
    abline(h=lineoff, lty=2)
    abline(h=-lineoff, lty=2)
  }
}
