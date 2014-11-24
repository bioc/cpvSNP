\name{glossiMarginal}
\alias{glossiMarginal}
\title{Calculate a Chi-squared Statistic and P-Value for Items
in a Set}
\usage{
glossiMarginal(pval, set)
}
\arguments{
\item{pval}{A vector of p-values of length equal to the number
of geneIds in the \code{GeneSetCollection} object.}

\item{set}{An element of a \code{GeneSetCollection} object.}
}
\value{
An object of type \code{GLOSSIResult}.
}
\description{
A helper function to calculate the chi-squared statistic
corresponding to an observed set of p-values.
}
\details{
This function calculates a chi-squared statistic from a set of
pvalues.
}
\author{
Jason Hackney, Jessica Larson
}
\examples{
    set.seed(30)
    pvals <- runif(100)
    names(pvals) <- paste0("rs", 1:100)
    snpGS5 <- GeneSet(geneIds=names(pvals)[1:5], setName="set5")
    res <- glossiMarginal (pvals, snpGS5)
}
