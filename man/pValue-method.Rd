\name{pValue-methods}
\docType{methods}
\alias{pValue}
\alias{pValue-methods}
\alias{pValue,GeneSetResult-method}
\alias{pValue,VEGASResult-method}
\alias{pValue,GLOSSIResult-method}
\alias{pValue,GeneSetResultCollection-method}
\alias{pValue,VEGASResultCollection-method}
\alias{pValue,GLOSSIResultCollection-method}
\title{ ~~ Methods for Function \code{pValue}  ~~}
\description{
This function returns the calculated p-value for a specified \code{GeneSet} 
or a list of p-values for a \code{GeneSetCollection}.  
}
\usage{
pValue(object)
}
\arguments{
  \item{object}{An object of type \code{GeneSetResult}, \code{VEGASResult}, \code{GLOSSIResult}, 
  	\code{GeneSetResultCollection}, \code{VEGASResultCollection}, or \code{GLOSSIResultCollection} }
}
\section{Methods}{
\describe{
 Defined methods include:
 \item{
      \code{signature(object = "VEGASResult")}
    }{
      Returns the calculated p-value for a specified \code{VEGASResult} object
    }
 \item{
      \code{signature(object = "GeneSetResult")}
    }{
      Returns the calculated p-value for a specified \code{GeneSetResult} object
    }
 \item{
      \code{signature(object = "GLOSSIResult")}
    }{
      Returns the calculated p-value for a specified \code{GLOSSIResult} object
    }
    
     \item{
      \code{signature(object = "VEGASResultCollection")}
    }{
      Returns a list of calculated p-values for a \code{VEGASResultCollection} object (1 for each set)
    }
     \item{
      \code{signature(object = "GeneSetResultCollection")}
    }{
      Returns a list of calculated p-values for a \code{GeneSetResultCollection} object (1 for each set)
    }
      \item{
      \code{signature(object = "GLOSSIResultCollection")}
    }{
      Returns a list of calculated p-values for a \code{GLOSSIResultCollection} object (1 for each set)
    }

	}
}
\section{Value}{
\describe{
Returns a decimal p-value corresponding to the gene set or a list of p-values corresponding to the gene set collection.
}
}
\seealso{
  \code{\linkS4class{VEGASResult}}-class, \code{geneSetName}
}
\examples{
	\dontrun{
	pValue( vegas(geneSet, assoc_table, ldMatrix) ) 
	}
}
\author{
Caitlin McHugh
}
\keyword{methods}
