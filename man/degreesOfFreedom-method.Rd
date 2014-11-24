\name{degreesOfFreedom-methods}
\docType{methods}
\alias{degreesOfFreedom}
\alias{degreesOfFreedom-methods}
\alias{degreesOfFreedom,VEGASResult-method}
\alias{degreesOfFreedom,GLOSSIResult-method}
\alias{degreesOfFreedom,GeneSetResult-method}
\alias{degreesOfFreedom,VEGASResultCollection-method}
\alias{degreesOfFreedom,GLOSSIResultCollection-method}
\alias{degreesOfFreedom,GeneSetResultCollection-method}
\title{ ~~ Methods for Function \code{degreesOfFreedom}  ~~}
\description{
This function returns the degrees of freedom used to calculate p-values from the corresponding \code{GeneSet} 
or a list of degrees of freedom for a \code{GeneSetCollection}.  
}
\usage{
degreesOfFreedom(object)
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
      Returns the degrees of freedom for a \code{GeneSetResult} object
    }
 \item{
      \code{signature(object = "VEGASResult")}
    }{
      Returns the degrees of freedom for a \code{VEGASResult} object
    }
 \item{
      \code{signature(object = "GLOSSIResult")}
    }{
      Returns the degrees of freedom for a \code{GLOSSIResult} object
    }
    
    \item{
      \code{signature(object = "GeneSetResultCollection")}
    }{
      Returns a list of degrees of freedom for a \code{GeneSetResultCollection} objects (1 for each set)
    }
     \item{
      \code{signature(object = "VEGASResultCollection")}
    }{
      Returns a list of degrees of freedom for a \code{VEGASResultCollection} objects (1 for each set)
    }
      \item{
      \code{signature(object = "GLOSSIResultCollection")}
    }{
      Returns a list of degrees of freedom for a \code{GLOSSIResultCollection} objects (1 for each set)
    }

	}
}
\section{Value}{
    \describe{
      Returns an integer or a list of integers indicating the degrees of freedom in the corresponding object.
    }
}
\seealso{
  \code{\linkS4class{VEGASResult}}-class, \code{pValue}
}
\examples{
	\dontrun{
	degreesOfFreedom( vegas(geneSet, assoc_table, ldMatrix) ) 
	}
}
\author{
Caitlin McHugh
}
\keyword{methods}
