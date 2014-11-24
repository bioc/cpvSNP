##' Plots association test p-values against number of SNPs per
#' gene set.
##'
##' Creates a plot of all p-values against the number of
#' SNPs in the gene set, for all sets within a collection.
##' @title Create a Plot of P-Values Against Number of SNPs
#' per Set for All Sets in a Collection
##' @param set An object of class \code{GeneSetResultCollection}. 
##' @param title A character string containing the title of the plot.
#' Default is "", a blank title.
##' @param xlab A character vector with the label for the x-axis on
#'  the plot. Default is \code{SNPs per Gene Set}.
##' @param ylab A character vector holding the label for the y-axis 
#' on the plot. Default is \code{PValue}.
##' @param ... Further arguments to be passed to the plotting methods,
#' such as graphical parameters.
##' @return Creates a plot.
##' @author Caitlin McHugh \email{mchughc@@uw.edu}
##' @export
plotPvals <- function(set, title="", xlab = "SNPs per Gene Set", 
ylab = "PValue", ...)
{
	if(!is.element(class(set),c("GeneSetResultCollection",
	    "GLOSSIResultCollection","VEGASResultCollection"))){
	        stop("set must be a GeneSetResultCollection object.")}

	result <- data.frame(df=unlist(degreesOfFreedom(set)),
	    pval=unlist(pValue(set)))
    ggplot(result, aes(x=df,y=pval, ...)) +
    geom_point() + xlab(xlab) + ylab(ylab) + ggtitle(title)
}