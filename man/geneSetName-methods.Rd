\name{geneSetName-methods}
\docType{methods}
\alias{geneSetName}
\alias{geneSetName-methods}
\alias{geneSetName,VEGASResult-method}
\alias{geneSetName,GLOSSIResult-method}
\alias{geneSetName,GeneSetResult-method}
\alias{geneSetName,VEGASResultCollection-method}
\alias{geneSetName,GLOSSIResultCollection-method}
\alias{geneSetName,GeneSetResultCollection-method}
\title{ ~~ Methods for Function \code{geneSetName}  ~~}
\description{
This function returns the gene set name for a specified \code{GeneSet} 
or a list of gene set names for a \code{GeneSetCollection}.  
}
\usage{
geneSetName(object)
}
\arguments{
  \item{object}{An object of type \code{GeneSetResult}, \code{VEGASResult}, \code{GLOSSIResult}, 
  	\code{GeneSetResultCollection}, \code{VEGASResultCollection}, or \code{GLOSSIResultCollection} }
}
\section{Methods}{
\describe{
 Defined methods include:
 \item{
      \code{signature(object = "GeneSetResult")}
    }{
      Returns the gene set name for a specified \code{GeneSetResult} object
    }
 \item{
      \code{signature(object = "VEGASResult")}
    }{
      Returns the gene set name for a specified \code{VEGASResult} object
    }
 \item{
      \code{signature(object = "GLOSSIResult")}
    }{
      Returns the gene set name for a specified \code{GLOSSIResult} object
    }
     \item{
      \code{signature(object = "GeneSetResultCollection")}
    }{
      Returns a list of gene set names for a \code{GeneSetResultCollection} object (1 for each set)
    }    
     \item{
      \code{signature(object = "VEGASResultCollection")}
    }{
      Returns a list of gene set names for a \code{VEGASResultCollection} object (1 for each set)
    }
      \item{
      \code{signature(object = "GLOSSIResultCollection")}
    }{
      Returns a list of gene set names for a \code{GLOSSIResultCollection} object (1 for each set)
    }

	}
}
\section{Value}{
\describe{
 Returns the name or a list of names of the gene sets in the corresponding object.
}
}

\seealso{
  \code{\linkS4class{VEGASResult}}-class, \code{pValue}
}
\examples{
	\dontrun{
	geneSetName( vegas(geneSet, assoc_table, ldMatrix) ) 
	}
}
\author{
Caitlin McHugh
}
\keyword{methods}
