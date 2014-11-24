##' Creates a GRanges object used for SNP set analysis.
##'
##' This function takes a data.frame and creates a GRanges 
#' object used for SNP set analysis.
##' @title Create a GRanges Object 
##' @param arrayData A data.frame object containing array data from 
#' which the \code{GRanges} object will be created. This object is expected
#' to include the position and chromosome, along with a SNP 
#' identifier and corresponding p-value.
##' @param positionName The name of the column in the filepathData 
#' object that holds position data for each probe. By default, this 
#' value is NULL and there are columns Start and End which hold 
#' this information.
##' @param chromosomeName The name of the column in the filepathData 
#' object that holds chromosome data for each probe. By default, 
#' this value is \code{chromosome}.
##' @param chromosomeNameConvention The naming convention used for the 
#' chromosomes, either \code{NCBI} (default), \code{UCSC}, or a user
#' specified mapping. By default, this is \code{NCBI}, indicating the autosomes
#' are stored as integers and the mitochondrial and sex chromosomes are stored
#' as \code{MT} and \code{X, Y}, respectively. The \code{UCSC} naming 
#' convention appends "chr" to each chromosome, and the mitochondrial coding is
#' "chrM." Alternatively, a user can supply the mapping from each chromosome
#' to the UCSC naming convention in the form of a character vector where the
#' elements are the UCSC names and the corresponding names of the elements are
#' the current chromosome names. An `NA' chromosome will be mapped to unknown.
#' For more information, please see the \code{GenomeInfoDb} BioConductor R package.
##' @param verbose A logical argument indicating whether output 
#' should be printed. The default is FALSE.
##' @return A GRanges object.
##' @author Jason Hackney, Jessica Larson, Caitlin McHugh
#' \email{mchughc@@uw.edu}
##' @export
createArrayData <- function(arrayData, 
    positionName=NULL, chromosomeName="chromosome", 
    chromosomeNameConvention="NCBI", verbose=TRUE){

    if(!is(arrayData, "data.frame")){ 
        stop("arrayData must be a data.frame object.")
    }

    names(arrayData)[is.element(names(arrayData),
        chromosomeName)] <- "chromosome"
    if(is.null(arrayData$chromosome)){
        stop("Data must have a `chromosome' variable.")    	
    }

    if(!is.element(chromosomeNameConvention, c("NCBI", "UCSC"))){
        # check to be sure its a vector with old/new names mapped
        chrs <- unique(arrayData$chromosome)
        if(!all(is.element(names(chromosomeNameConvention)),chrs)){
            stop("There is an error in the chromosomeNameConvention argument.")
        }
        if(length(chromosomeNameConvention)!=length(chrs)){
        	stop("The chromosomeNameConvention argument does not include all
        	    chromosomes listed in the data.")
        }
    }   

    if(!is.null(arrayData$FILTER)){
        arrayData$FILTER <- factor(arrayData$FILTER)
    }
    if(!is.null(arrayData$REF)){arrayData$REF <- factor(arrayData$REF)}
    if(!is.null(arrayData$ALT)){arrayData$ALT <- factor(arrayData$ALT)}
    if(is.null(positionName)&(is.null(arrayData$Start)|
        is.null(arrayData$End))){
        stop("Data must have a `position' variable.")
    }
    if(!is.null(positionName)){
        thisCol <- which(names(arrayData)==positionName)
        arrayData$Start <- arrayData[,thisCol]
        arrayData$End <- arrayData[,thisCol] 
    }
            
    # subset to only those that have calculated P-values
    arrayData <- arrayData[!is.na(arrayData$P),]
            
    arrayData$Start[is.na(arrayData$chromosome)] <- 1
    arrayData$End[is.na(arrayData$chromosome)] <- 1
                
    arrayData$chromosome[is.na(arrayData$chromosome)] <- 'U'
    
    if(chromosomeNameConvention=="NCBI"){
        arrayData$chromosome[is.element(arrayData$chromosome,0)] <- 'U'
        arrayData$chromosome[is.element(arrayData$chromosome,27)] <- 'U'
        arrayData$chromosome[is.element(arrayData$chromosome,"MT")] <- 'M'    
           
        arrayData$chromosome <- factor(arrayData$chromosome)
        levels(arrayData$chromosome) <- paste('chr',
            levels(arrayData$chromosome), sep = '')
    }
    
    if(!is.element(chromosomeNameConvention,c("NCBI","UCSC"))){
    	chrFact <- as.factor(arrayData$chromosome)
    	arrayData$chromosome <- chromosomeNameConvention[chrFact]
    }
        
    arrayData$chromosome <- factor(arrayData$chromosome,
        levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
        "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
        "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
        "chr22", "chrX", "chrY", "chrM", "chrXY", "chrU"))
                
    if(verbose){ message("making GRanges object") }
     
    arrayData <- arrayData[order(arrayData$chromosome, arrayData$Start),] 
            
    array.ranges <- GRanges( Rle(arrayData$chromosome), 
        IRanges(start=arrayData$Start, end=arrayData$End) )
    elementMetadata(array.ranges) <- arrayData

    return(array.ranges)
}