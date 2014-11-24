\name{VEGASResult-class}
\Rdversion{1.1}
\docType{class}
\alias{VEGASResult-class}

\title{Class \code{"VEGASResult"}}
\description{
Objects of this class store results from running VEGAS. 

}
\section{Objects from the Class}{
Objects can be created by calls of \code{vegas}.
}
\section{Slots}{
  \describe{
  	\item{\code{geneSetName}:}{Object of class \code{"character"}, the name of the geneSet }
    \item{\code{pValue}:}{Object of class \code{"numeric"}, the p-value }
    \item{\code{degreesOfFreedom}:}{Object of class \code{"integer"}, the degrees of freedom used to calculate the p-value }
    \item{\code{statistic}:}{Object of class \code{"numeric"}, the test statistic for the set }
    \item{\code{simulatedStats}:}{Object of class \code{"vector"}, the simulated statistics for the set used to calculate the p-value }
  }
}
\section{Methods}{
No methods defined with class "VEGASResult" in the signature.
}
\section{Accessors}{
\describe{
In the code snippets below, \code{object} is a VEGASResult object.
    \item{\code{geneSetName(object)}:}{The name of the gene set.}
    \item{\code{pValue(object)}:}{A numeric p-value.}
    \item{\code{degreesOfFreedom(object)}:}{Integer value of the degrees of freedom used in the statistical test.}
    \item{\code{statistic(object)}:}{A numeric value corresponding to the statistic calculated from the test performed.}
}
}
\author{
Caitlin McHugh
}

\seealso{
\code{\linkS4class{VEGASResultCollection}}, \code{\linkS4class{GLOSSIResult}} 
}
\examples{
showClass("VEGASResult")
}
\keyword{classes}
