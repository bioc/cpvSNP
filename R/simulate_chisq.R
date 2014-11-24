##' Simulates a value from a chi-squared distribution with a 
#' specified correlation structure.
##'
##' This function simulates a value from a chi-squared 
#' distribution with a specified correlation structure and degrees
#' of freedom corresponding to the number of rows in the 
#' correlation matrix.
##' @title Simulate a Chi-squared Statistic from a Distribution 
#' with a Specified Correlation Structure
##' @param corr_matrix A square matrix of correlations.
##' @param num_sims An integer value for the number of simulations
#' to be performed. 
##' @return A numeric value corresponding to a chi-squared 
#' statistic simulated from a chi-squared-like distribution with 
#' specified correlation structure and degrees of freedom 
#' corresponding to the number of rows in the correlation matrix.
##' @author Caitlin McHugh \email{mchughc@@uw.edu}
##' @export
simulate_chisq <- function(corr_matrix, num_sims)
{
    # sample n random standard normal values
    z <- rnorm(nrow(corr_matrix)*num_sims,mean=0,sd=1)
    zmat <- matrix(z,nrow=num_sims)
    # multiply them by the correlation structure
    z_corr <- zmat%*%corr_matrix
    
    # square and sum them
    zsq <- z_corr**2
    z_chisq <- rowSums(zsq)
    
    return(z_chisq)
}