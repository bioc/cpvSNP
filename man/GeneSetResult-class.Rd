\name{GeneSetResult-class}
\Rdversion{1.1}
\docType{class}
\alias{GeneSetResult-class}

\title{Class \code{"GeneSetResult"}}
\description{
Objects of this class store results from running GeneSet methods. 

}
\section{Objects from the Class}{
Objects can be created by calls of \code{glossi} or \code{vegas}.
}
\section{Slots}{
  \describe{
  	\item{\code{geneSetName}:}{Object of class \code{"character"}, the name of the geneSet }
    \item{\code{pValue}:}{Object of class \code{"numeric"}, the p-value }
    \item{\code{degreesOfFreedom}:}{Object of class \code{"integer"}, the degrees of freedom used to calculate the p-value }
    \item{\code{statistic}:}{Object of class \code{"numeric"}, the test statistic for the set }
  }
}
\section{Methods}{
No methods defined with class "GeneSetResult" in the signature.
}
\section{Accessors}{
\describe{
In the code snippets below, \code{object} is a GeneSetResult object.
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
\code{\linkS4class{GeneSetResultCollection}}, \code{\linkS4class{GLOSSIResult}} 
}
\examples{
showClass("GeneSetResult")
}
\keyword{classes}
