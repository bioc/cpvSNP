##' Translates a set of gene ids to their corresponding SNPs.
##'
##' Translates a GeneSet or GeneSetCollection object of gene 
#' Entrez ids to the corresponding SNPs that lie within a 
#' prespecified region of the gene.
##' @title Translate a List of Gene Ids to Their Corresponding 
#' SNP Ids
##' @param geneList An object of class \code{GeneSetCollection}
#'  or \code{GeneSet} where geneIds holds the gene Entrez ids.
##' @param arrayData A \code{GRanges} object that corresponds to 
#' the probes for which data are present. 
##' @param genes A \code{GRanges} object holding the gene transcripts
#' fow which the mapping will occur. The ids of the genes should match
#' the ids that are listed in the \code{geneList} object.
##' @param maxgap An integer value indicating the flanking region 
#' around the gene for which SNPs are considered part of the gene. 
#' Default is 20000 (20kb). 
##' @param quiet A logical variable indicating whether warnings should
#' be output from the SNP to gene mapping step.
##' @return A \code{GeneSetCollection} object containing all SNP ids that 
#' lie within the genes listed in geneList.
##' @author Jason Hackney, Jessica Larson, 
#' Caitlin McHugh \email{mchughc@@uw.edu}
##' @export
geneToSNPList <- function(geneList, arrayData, genes,
    maxgap=20000, quiet=TRUE)
{
    if(!is(genes, "GRanges")){stop("genes must be a GRanges object.")}
	
	if(!is.logical(quiet)){stop("quiet must be a logical argument.")}
	
    if(quiet){
        hitsRes <- suppressWarnings(findOverlaps(genes, arrayData, 
        maxgap = maxgap))
    }else{hitsRes <- findOverlaps(genes, arrayData, maxgap = maxgap)}
    
    hitsMat <- as.data.frame(hitsRes)
    spl <- split(hitsMat$subjectHits,names(genes)[hitsMat$queryHits])

    # this transforms the elements of the split from indices to the snp at 
        # that index
    spl2 <- bplapply(spl,function(x) arrayData$SNP[x])

    # this creates a gene set collection object that maps identifiers from 
    # the geneList to the AnnotationIdentifier() from the spl2 list 
    snp.gsc <- GeneSetCollection(bplapply(geneList,mapIdentifiers, 
        AnnotationIdentifier(), as.environment(spl2)))
    
    return(snp.gsc)
}