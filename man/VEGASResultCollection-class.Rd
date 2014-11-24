\name{VEGASResultCollection-class}
\Rdversion{1.1}
\docType{class}
\alias{VEGASResultCollection-class}

\title{Class \code{"VEGASResultCollection"}}
\description{
Objects of this class store results from running VEGAS. VEGASResultCollection objects contain a list of VEGASResult objects (one result for each GeneSet)
}
\section{Objects from the Class}{
Objects can be created by calls of \code{vegas}.
}
\section{Methods}{
No methods defined with class "VEGASResult" in the signature.
}
\section{Accessors}{
\describe{
In the code snippets below, \code{object} is a VEGASResultCollection object.
    \item{\code{geneSetName(object)}:}{A list of the names of the gene sets.}
    \item{\code{pValue(object)}:}{A list of numeric p-values.}
    \item{\code{degreesOfFreedom(object)}:}{List of integer values of the degrees of freedom used in the statistical tests performed on each gene set.}
    \item{\code{statistic(object)}:}{A list of numeric values corresponding to the statistics calculated from the tests performed on each gene set.}
}
}
\author{
Caitlin McHugh
}

\seealso{
\code{\linkS4class{VEGASResult}}, \code{\linkS4class{GLOSSIResultCollection}} 
}
\examples{
showClass("VEGASResultCollection")
}
\keyword{classes}
