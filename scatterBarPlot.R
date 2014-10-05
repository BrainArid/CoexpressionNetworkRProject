scatterBarPlot <- function(x, y, dcol="blue", lhist=20, num.dnorm=5*lhist, ...){
  ## check input
  ##stopifnot(ncol(x)==2)
  ## set up layout and graphical parameters
  layMat <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(layMat, widths=c(5/7, 2/7), heights=c(2/7, 5/7))
  ospc <- 0.5 # outer space
  pext <- 4 # par extension down and to the left
  bspc <- 1 # space between scatter plot and bar plots
  par. <- par(mar=c(pext, pext, bspc, bspc), oma=rep(ospc, 4)); # plot parameters
  ## scatter plot
  plot(x,y, xlim=range(is.finite(x)), ylim=range(is.finite(y)), ...)
  ## 3) determine barplot and height parameter
  ## histogram (for barplot-ting the density)
  hist(x, plot=FALSE, breaks=seq(from=min(x), to=max(x),length.out=lhist));
  hist(y, plot=FALSE, breaks=seq(from=min(y), to=max(y),length.out=lhist)); # note: this uses probability=TRUE
  ## determine the plot range and all the things needed for the barplots and lines
  ##xx <- seq(min(x, max(x), length.out=num.dnorm)) # evaluation points for the overlaid density
  ###xy <- dnorm(xx, mean=mean(x), sd=sd(x)) # density points
  ##yx <- seq(min(y, max(y), length.out=num.dnorm))
  ##yy <- dnorm(yx, mean=mean(y), sd=sd(y))
  ## barplot and line for x (top)
  par(mar=c(0, pext, 0, 0))
  barplot(xhist$density, axes=FALSE, ylim=c(0, xhist$density),
          space=0) # barplot
  lines(seq(from=0, to=lhist-1, length.out=num.dnorm), xy, col=dcol) # line
  ## barplot and line for y (right)
  par(mar=c(pext, 0, 0, 0))
  barplot(yhist$density, axes=FALSE, xlim=c(0, yhist$density),
          space=0, horiz=TRUE) # barplot
  lines(yy, seq(from=0, to=lhist-1, length.out=num.dnorm), col=dcol) # line
  ## restore parameters
  par(par.)
}