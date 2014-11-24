
\name{geneToSNPList}
\alias{geneToSNPList}
\title{Translate a List of Gene Ids to Their Corresponding
SNP Ids}
\usage{
geneToSNPList(geneList, arrayData, genes, maxgap = 20000, quiet = TRUE)
}
\arguments{
\item{geneList}{An object of class \code{GeneSetCollection}
or \code{GeneSet} where geneIds holds the gene Entrez ids.}

\item{arrayData}{A \code{GRanges} object that corresponds to
the probes for which data are present.}

\item{genes}{A \code{GRanges} object holding the gene transcripts
for which the mapping will occur. The ids of the genes should match
the ids that are listed in the \code{geneList} object.}

\item{maxgap}{An integer value indicating the flanking region
around the gene for which SNPs are considered part of the gene.
Default is 20000 (20kb).}

\item{quiet}{A logical variable indicating whether warnings should
be output from the SNP to gene mapping step.}
}
\value{
A \code{GeneSetCollection} object containing all SNP ids that
lie within the genes listed in geneList.
}
\description{
Translates a set of gene ids to their corresponding SNPs.
}
\details{
Translates a GeneSet or GeneSetCollection of gene
Entrez ids to the corresponding SNPs that lie within a
prespecified region of the gene.
}
\author{
Jason Hackney, Jessica Larson,
Caitlin McHugh \email{mchughc@uw.edu}
}
\examples{
    data(geneSetAnalysis)
    head(geneSetAnalysis$arrayData)
    arrayDataGR <- createArrayData(geneSetAnalysis[["arrayData"]], positionName="Position")
    
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    genesHg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
	geneSets <- geneSetAnalysis[["geneSets"]]
	
	snpsGSC <- geneToSNPList(geneSets, arrayDataGR, genesHg19)
}
\keyword{methods}

