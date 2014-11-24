\name{geneSetAnalysis}
\docType{data}
\alias{geneSetAnalysis}
\title{Data to run gene set analysis methods}
\description{
	A list including the following components:

	\code{arrayData} -- a data.frame of dimension 119487 rows and 4 columns including P-values, SNP ids, position and chromosome for a GWAS experiment;

	\code{geneSets} -- a GeneSetCollection containing gene Entrez ids;

	\code{ldMat} -- a matrix of dimension 871 x 871 holding linkage disequilibrium values;

	\code{indepSNPs} -- a list containing 302 independent SNP ids;
}
\usage{geneSetAnalysis}
\examples{
    data(geneSetAnalysis)
    head(geneSetAnalysis$arrayData)
}
\keyword{datasets}
