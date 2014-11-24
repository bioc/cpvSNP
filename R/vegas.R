##' Calculates the p-value representing the association of the set
#' with the phenotype of interest.
##'
##' This function calculates a p-value for sets of SNPs that 
#' reside within a gene set collection. We calculate the null 
#' distribution by taking into account the observed correlation 
#' among the SNPs and simulating a specified number of statistics
#' from which the resulting p-value is calculated.
##' @title Calculate the P-Value for a Set Using the VEGAS Method
##' @param set A \code{GeneSet} object containing a set of genes 
#' of interest or a \code{GeneSetCollection} object containing a 
#' collection of GeneSets.
##' @param assoc_table An object of class \code{GRanges}. This 
#' object should at least contain columns SNP and P which hold SNP
#' rsIDs and their corresponding association test p-values, 
#' respectively.
##' @param ldMatrix A square, symmetric matrix of LD values, with
#' each row and column corresponding to each of the items in the 
#' set. The diagonal entries should be 1, indicating the LD 
#' between an item in the set and itself is 1.
##' @param num_sims A positive integer value for the number of
#' simulations to be performed. Default is 1000.
##' @param correction A logical argument indicating whether a 
#' value of one should be added to the numerator when calculating 
#' the p-value based upon the simulated statistics. By default, 
#' the correction is added and this argument is TRUE.
##' @param seed An integer argument indicating what the random
#' seed should be set to. This allows for replication of results.
#' The default is NULL, and a random seed will be set internally. 
##' @param verbose A logical argument indicating whether output 
#' should be printed. The default is FALSE.
##' @return An object with the corresponding VEGAS results. If set
#' is a \code{GeneSetCollection} (i.e., multiple sets of 
#' interest), then the corresponding \code{VEGASResultCollection}
#' is returned. If set is a \code{GeneSet}, a \code{VEGASResult} 
#' object will be returned.
##' @references Liu, Jimmy Z. and Mcrae, Allan F. et al. 
#' A Versatile Gene-Based Test for Genome-Wide Association 
#' Studies.
#' The American Journal of Human Genetics, 2010.
##' @author Caitlin McHugh \email{mchughc@@uw.edu}
##' @export
vegas <- function(set, assoc_table, ldMatrix, num_sims=1000, 
    correction=TRUE, seed=NULL, verbose=FALSE)
{
    if(!is.logical(verbose)) stop("verbose must be either TRUE 
        or FALSE.")
    if(!is.logical(correction)) stop("correction must be either 
        TRUE or FALSE.")
    if(!all.equal(num_sims,as.integer(num_sims))) stop("num_sims 
        must be an integer.")  
    if(num_sims<=0){stop("num_sims must be positive.")}  
    
    if(num_sims<100 & verbose){
    	print("Results from running less than 100 simulations should be
            interpreted carefully.")
    }
    
    if(!(is.null(seed)|is.numeric(seed))) stop("seed must be a numeric value.")
    
    if(!is(ldMatrix, "data.frame")){
    	ldMatrix <- as.matrix(ldMatrix)
    	message("Warning: coercing ldMatrix from data.frame to matrix.")
    }
    
    # first, remove SNPs with all NA values in the LD matrix
    naSNPS <- apply(ldMatrix,1,function(x){sum(is.na(x))})
    torm <- names(naSNPS)[naSNPS>(nrow(ldMatrix)-3)]
    if(length(torm)>0){
	    colsKp <- !is.element(rownames(ldMatrix),torm)
    	ldMatrix <- ldMatrix[colsKp,colsKp]
    }
    
    # then remove SNPs with just one NA value in the LD matrix
    naSNPS <- apply(ldMatrix,1,function(x){sum(is.na(x))})
    torm <- names(naSNPS)[naSNPS>0]
    colsKp <- !is.element(rownames(ldMatrix),torm)
    ldMatrix <- ldMatrix[colsKp,colsKp]

    if(nrow(ldMatrix)!=ncol(ldMatrix)) stop("ldMatrix needs to 
        be square.")
    if(!isSymmetric(ldMatrix,tol=1e-04)) stop("ldMatrix needs to 
        be symmetric.")

    di <- diag(ldMatrix)
    if(!all(round(di,digits=5)==1)) stop("ldMatrix should have 
        all 1's on the diagonal.")
    
    if(is(set)[1]=="GeneSetCollection"){
    # loop through the gene sets, get the analysis results
        res <- bplapply(set,vegasPrep,assoc_table=assoc_table,
            ldMatrix=ldMatrix,num_sims=num_sims,
	        correction=correction,seed=seed,verbose=verbose)
	    return(new("VEGASResultCollection", res))
    }else{
    	res <- vegasPrep(set,assoc_table=assoc_table,
    	    ldMatrix=ldMatrix,num_sims=num_sims,
            correction=correction,seed=seed,verbose=verbose)
    	return(res)
    }
}