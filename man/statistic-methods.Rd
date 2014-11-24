\name{statistic-methods}
\docType{methods}
\alias{statistic}
\alias{statistic-methods}
\alias{statistic,GeneSetResult-method}
\alias{statistic,GeneSetResultCollection-method}
\alias{statistic,VEGASResult-method}
\alias{statistic,VEGASResultCollection-method}
\alias{statistic,GLOSSIResult-method}
\alias{statistic,GLOSSIResultCollection-method}
\title{ ~~ Methods for Function \code{statistic}  ~~}
\description{
This function returns the calculated statistic for a specified \code{GeneSet} 
or a list of calculated statistics for a \code{GeneSetCollection}.  
}
\usage{
statistic(object)
}
\arguments{
  \item{object}{An object of type \code{GeneSetResult}, \code{VEGASResult}, \code{GLOSSIResult},
  	 \code{GeneSetResultCollection}, \code{VEGASResultCollection} or \code{GLOSSIResultCollection} }
}
\section{Methods}{
\describe{
 Defined methods include:
 \item{
      \code{signature(object = "GeneSetResult")}
    }{
      Returns the calculated statistic for a specified \code{GeneSetResult} object
    }
      \item{
      \code{signature(object = "GeneSetResultCollection")}
    }{
      Returns a list of calculated statistics for a \code{GeneSetResultCollection} object (1 for each set)
    }
 \item{
      \code{signature(object = "VEGASResult")}
    }{
      Returns the calculated statistic for a specified \code{VEGASResult} object
    }
      \item{
      \code{signature(object = "VEGASResultCollection")}
    }{
      Returns a list of calculated statistics for a \code{VEGASResultCollection} object (1 for each set)
    }
 \item{
      \code{signature(object = "GLOSSIResult")}
    }{
      Returns the calculated statistic for a specified \code{GLOSSIResult} object
    }
      \item{
      \code{signature(object = "GLOSSIResultCollection")}
    }{
      Returns a list of calculated statistics for a \code{GLOSSIResultCollection} object (1 for each set)
    }

	}
}
\section{Value}{
\describe{
Returns decimal values or a list of decimal values corresponding to calculated statistics from a GLOSSI function call.
}
}
\seealso{
  \code{\linkS4class{GLOSSIResult}}-class, \code{pValue}
}
\examples{
	\dontrun{
	statistic( glossi(geneSet, assoc_table, ldMatrix) ) 
	}
}
\author{
Caitlin McHugh
}
\keyword{methods}
