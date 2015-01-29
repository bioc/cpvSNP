##' Plots association test p-values, highlighting p-values within
#'  gene sets.
##'
##' Creates a density plot of all -log10(p-value), and overlaying 
#' a density plot of -log10(p-value) for SNPs within a gene set 
#' of interest.
##' @title Create a Density Plot of P-Values, Highlighting SNPs 
#' Within a Gene Set
##' @param pvals A numerical vector of p-values with names 
#' corresponding to elements listed in the geneIds slot in the 
#' set object.
##' @param set An object of class \code{GeneSet} where geneIds holds 
#' the items of each set corresponding to the pvals.
##' @param title A character string containing the title of the plot.
#' Default is "", a blank title.
##' @param geneSetColor A character vector holding the color the 
#' highlighted SNPs in the set object should be. Default is 
#' \code{red}.
##' @param xlab A character vector with the label for the x-axis on
#'  the plot. Default is \code{-log10(PValue)}.
##' @param ylab A character vector holding the label for the y-axis 
#' on the plot. Default is \code{Density}.
##' @param xlim A vector of two numeric values with the limits of 
#' plotting on the x-axis. By default, the range of the x-axis 
#' values will be determined from the data.
##' @param ylim A vector with two numeric values used as limits for 
#' plotting on the y-axis. By default, the range of the y-axis values
#'  will be determined from the data.
##' @param ... Further arguments to be passed to the plotting methods,
#' such as graphical parameters.
##' @return Creates a plot.
##' @author Jason Hackney, Jessica Larson, Caitlin McHugh
#' \email{mchughc@@uw.edu}
##' @export
assocPvalBySetPlot <- function(pvals, set, title="", 
    geneSetColor = "red", 
    xlab = expression(-log[10]~p-value), 
    ylab = "Density",
    xlim=NULL, ylim=NULL, ...)
{
    if(!is(set)[1]=="GeneSet"){
        stop("set must be of class GeneSet.") }
	
    geneSetColor <- rgb(t(col2rgb(geneSetColor)), alpha=90, 
        maxColorValue=255)
        
	setlocs <- names(pvals) %in% geneIds(set)
    setStats <- pvals[setlocs]	
	
    if(length(setlocs)<1){
        stop("No pvalues included for ids listed in geneIds(set).") }
	
    data <- list("main" = -log10(pvals+1e-09), 
        "geneSet" = -log10(setStats+1e-09))
    ## data is a list of numerical vectors: 'main', 'geneSet'
    dens <- bplapply(data , function(x) { if( length( x ) > 1 ){
        density( x, na.rm = TRUE )
	} else { NULL } } )
    if(is.null(xlim)) {
        xlim <- range(sapply(dens, function(x) {
	    if(!is.null(x)) { range(x$x, na.rm = TRUE) }
	    else { return(NA) }
            }), na.rm =TRUE)
    }
    if ( is.null(ylim) ) {
        ylim <- range(sapply(dens, function(x) {
	    if(!is.null(x)) {
	    range(x$y, na.rm = TRUE)
            } else { return(NA) }
	    }), na.rm =TRUE)
    }
    xlim <- c(0, max(abs(xlim)))
    for(n in setdiff(c("main", "geneSet"), names(dens))){ 
        dens[n]= NA }

    plot(dens$main, xlim = xlim, ylim = ylim, xlab = xlab, 
        main=title, ylab = ylab, col = "grey75", xaxt='n', ...)
	legend("topright", c("all assoc results","Set of interest"),
        col=c("grey75","red"), lty=c(1,1), lwd=c(3,3))

    polygon(dens$main, col="grey75", border="grey75")

    lines(dens$geneSet, xlim = xlim, ylim = ylim, xlab = xlab, 
        ylab = ylab, col = geneSetColor, xaxt='n', lwd=3)
    points(data$geneSet, y=jitter(rep(0,length(data$geneSet)), factor=0.1),
        cex=1, pch=19, col=rgb(t(col2rgb(geneSetColor)), alpha=100,
            maxColorValue=255))
}