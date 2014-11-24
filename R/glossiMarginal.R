##' A helper function to calculate the chi-squared statistic 
#' corresponding to an observed set of p-values.
##'
##' This function calculates a chi-squared statistic from a set of
#' pvalues.
##' @title Calculate a Chi-squared Statistic and P-Value for Items
#' in a Set
##' @param pval A vector of p-values of length equal to the number
#' of geneIds in the \code{GeneSetCollection} object.
##' @param set An element of a \code{GeneSetCollection} object.
##' @return An object of type \code{GLOSSIResult}.
##' @author Jason Hackney, Jessica Larson
##' @export
glossiMarginal <- function(pval, set)
{
    ids <- intersect(unlist(geneIds(set)), names(pval))
    
    # need to add a bit so we don't have log(0)
    gamma.stat <- sum(-log(pval[ids]+.Machine$double.eps) )
    
    # find the corresponding gamma/chisq density pvalue
    gamma.p <- pgamma(gamma.stat, length(ids), lower.tail=FALSE)
    
    return( new("GLOSSIResult",
        geneSetName = setName(set),
        pValue = gamma.p,
        degreesOfFreedom = length(ids),
        statistic = gamma.stat) )
}