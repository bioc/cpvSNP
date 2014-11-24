##' Calculates the p-value representing the association of the 
#' set with the phenotype of interest, assuming all items in the
#' set are independent, using the GLOSSI method.
##'
##' This function calculates a p-value for sets of SNPs that
#' reside within a gene set collection. We calculate the chi-
#' square p-values and statistic by applying Fisher's 
#' transformation to the observed p-values. 
##' @title Calculate a Chi-squared Statistic and P-Value for
#' Independent Items of a Set Using the GLOSSI Method
##' @param pvals A numerical vector of p-values with names 
#' corresponding to elements listed in the \code{geneIds} slot in 
#' the snp.gsc object.
##' @param snp.gsc An object of class \code{GeneSetCollection} 
#' where geneIds holds the items of each set corresponding to the 
#' pvals.
##' @return An object with the corresponding GLOSSI results. If 
#' snp.gsc is a \code{GeneSetCollection} (i.e., multiple sets of 
#' interest), then the corresponding \code{GLOSSIResultCollection}
#' is returned. If snp.gsc is a \code{GeneSet}, a \code{GLOSSIResult}
#' object will be returned.
##' @references Chai, High-Seng and Sicotte, Hughes et al. 
#' GLOSSI: a method to assess the association of genetic loci-sets 
#' with complex diseases.
#' BMC Bioinformatics, 2009.
##' @author Jason Hackney, Jessica Larson, Caitlin McHugh 
#' \email{mchughc@@uw.edu}
##' @export
glossi <- function(pvals, snp.gsc)
{
    if(is.null(names(pvals))) stop("Names of p-values must 
        correspond to SNPs.")
    
    # subset the gene stats to only include genes in the gene set    
    pvals <- pvals[!is.na(pvals)]
    pathway.snps <- unique(unlist(geneIds(snp.gsc)))
    common.snps <- intersect(names(pvals), pathway.snps)
    pvals <- pvals[common.snps]
    set.names <- names(snp.gsc)
    
    if(is(snp.gsc)[1]=="GeneSetCollection"){
    	res <- sapply( set.names, function(x){glossiMarginal(pvals,
            snp.gsc[[x]])} )
    	return( new("GLOSSIResultCollection", res) )
    }
    
    if(is(snp.gsc)[1]=="GeneSet"){
    	res <- glossiMarginal(pvals,snp.gsc)
    	return(res)
    }
    
}