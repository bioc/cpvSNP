\name{simulatedStats-methods}
\docType{methods}
\alias{simulatedStats}
\alias{simulatedStats-methods}
\alias{simulatedStats,VEGASResult-method}
\alias{simulatedStats,VEGASResultCollection-method}
\title{ ~~ Methods for Function \code{simulatedStats}  ~~}
\description{
This function returns the simulated statistics for a specified \code{GeneSet} 
or a list of simulated statistics for a \code{GeneSetCollection}.  
}
\usage{
simulatedStats(object)
}
\arguments{
  \item{object}{An object of type \code{VEGASResult} or 
  	\code{VEGASResultCollection} }
}
\section{Methods}{
\describe{
 Defined methods include:
 \item{
      \code{signature(object = "VEGASResult")}
    }{
      Returns the simulated statistics for a specified \code{VEGASResult} object
    }    
     \item{
      \code{signature(object = "VEGASResultCollection")}
    }{
      Returns a list of simulated statistics for a \code{VEGASResultCollection} object (1 for each set)
    }
	}
}
\section{Value}{
\describe{
  Returns decimal values of simulated statistics from a VEGAS function call or a list of simulated statistics from a VEGAS function call on a gene set collection.
}
}
\seealso{
  \code{\linkS4class{VEGASResult}}-class, \code{observedStat}
}
\examples{
	\dontrun{
	simulatedStats( vegas(geneSet, assoc_table, ldMatrix) ) 
	}
}
\author{
Caitlin McHugh
}
\keyword{methods}
