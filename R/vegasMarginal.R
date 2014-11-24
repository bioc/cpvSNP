##' Calculates the p-value representing the association of the 
#' set with the phenotype of interest.
##'
##' This is a helper function to calculate the p-value for a set 
#' of SNPs that reside within a gene set collection. The 
#' correlation among the SNPs is taken into account by the LD 
#' matrix. The resulting p-value is calculated from a null 
#' distribution that is simulated \code{num_sims} times based upon
#' the specified correlation structure.
##' @title Calculate the P-Value for a Set Using the VEGAS Method
##' @param pvals A vector of p-values corresponding to items in 
#' the set.
##' @param ld_matrix A square, symmetric matrix of LD values, with
#' each row and column corresponding to each of the items in the 
#' set. The diagonal entries should be 1, indicating the LD 
#' between an item in the set and itself is 1.
##' @param num_sims An integer value for the number of simulations
#' to be performed. 
##' @param correction A logical argument indicating whether a 
#' value of one should be added to the numerator when calculating 
#' the p-value based upon the simulated statistics. By default, 
#' the correction is added. An argument of FALSE will not add one
#' to the numerator.
##' @param seed An integer argument indicating what the random
#' seed should be set to. This allows for replication of results.
#' The default is NULL, and a random seed will be set internally. 
##' @param verbose A logical argument indicating whether periodic
#' output should be printed. Defaults to FALSE, indicating no 
#' output will be printed.
##' @return A \code{VEGASResult} object with the corresponding 
#' VEGAS results.
##' @author Caitlin McHugh \email{mchughc@@uw.edu}
##' @export
vegasMarginal <- function(pvals, ld_matrix, num_sims, 
    correction=TRUE, seed=NULL, verbose=FALSE)
{
    if(!length(pvals)==nrow(ld_matrix)){
        stop("pvals should correspond to the entries 
            in ldMatrix.")
    }
    if(nrow(ld_matrix)==0) stop("no SNPs in this gene set.")
    
    if(is.null(seed)) seed <- rbinom(1,1e6,0.5)
	set.seed(seed)
	
    ld_matrix <- make.positive.definite(ld_matrix)

    # take transpose of the the cholesky decomp of the ld_matrix
    chol_matrix <- t(chol(ld_matrix))
    
    # calculate the observed sum of chisqs for this gene set
    chisqVals <- qchisq(pvals,df=1,ncp=0,lower.tail=FALSE)
    chisq_obs <- sum(chisqVals)
    
    # simulate num_sims chisq stats with the same correlation
    chisq_sim <- simulate_chisq(chol_matrix,num_sims)
    
    output <- vector("list",3)
    names(output) <- c("pvalue","statistic","simulated.stats")
    numerP <- sum(chisq_sim>chisq_obs)
	if(correction){
        output[["pvalue"]] <- (numerP+1)/(num_sims+1)
	}else{
        output[["pvalue"]] <- numerP/num_sims
	}
    output[["statistic"]] <- chisq_obs
    output[["simulated.stats"]] <- chisq_sim
    
    return(output)
}