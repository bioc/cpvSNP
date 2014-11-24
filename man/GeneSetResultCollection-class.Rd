\name{GeneSetResultCollection-class}
\Rdversion{1.1}
\docType{class}
\alias{GeneSetResultCollection-class}

\title{Class \code{"GeneSetResultCollection"}}
\description{
Objects of this class store results from running gene set methods such as GLOSSI or VEGAS. GeneSetResultCollection objects contain a list of GeneSetResult objects (one result for each GeneSet)
}
\section{Objects from the Class}{
Objects can be created by calls of \code{glossi} or \code{vegas}.
}
\section{Methods}{
No methods defined with class "GeneSetResultCollection" in the signature.
}
\section{Accessors}{
\describe{
In the code snippets below, \code{object} is a GeneSetResultCollection object.
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
\code{\linkS4class{GeneSetResult}}, \code{\linkS4class{GLOSSIResultCollection}} 
}
\examples{
showClass("GeneSetResultCollection")
}
\keyword{classes}
