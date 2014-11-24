\name{GLOSSIResultCollection-class}
\Rdversion{1.1}
\docType{class}
\alias{GLOSSIResultCollection-class}

\title{Class \code{"GLOSSIResultCollection"}}
\description{
Objects of this class store results from running GLOSSI. GLOSSIResultCollection objects contain a list of GLOSSIResult objects (one result for each GeneSet)
}
\section{Objects from the Class}{
Objects can be created by calls of \code{glossi}.
}
\section{Methods}{
No methods defined with class "GLOSSIResult" in the signature.
}
\section{Accessors}{
\describe{
In the code snippets below, \code{object} is a GLOSSIResultCollection object.
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
\code{\linkS4class{GLOSSIResult}}, \code{\linkS4class{VEGASResultCollection}} 
}
\examples{
showClass("GLOSSIResultCollection")
}
\keyword{classes}
