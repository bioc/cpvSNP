vegasPrep <- function(set, assoc_table, ldMatrix, num_sims,
    correction, seed, verbose){
	
    if(is(set)[1]=="GeneSetCollection"){gsn <- names(set)}
    if(is(set)[1]=="GeneSet"){gsn <- setName(set)}

    if(verbose){ message(paste("now analyzing gene set",gsn)) }
    snpsKp <- is.element(rownames(ldMatrix),
        unlist(geneIds(set)))&is.element(rownames(ldMatrix),
        assoc_table$SNP)
    this_ld <- ldMatrix[snpsKp,snpsKp]
	
    assoc_thisSet <- assoc_table[is.element(assoc_table$SNP,
        rownames(this_ld)),]
    this_vegas_result <- vegasMarginal(assoc_thisSet$P,this_ld,
        num_sims,correction,seed,verbose)
    if(is.null(this_vegas_result)){ return(NULL) }
	
    return( new("VEGASResult",
        geneSetName = gsn,
	    pValue = this_vegas_result[["pvalue"]],
	    degreesOfFreedom = dim(elementMetadata(assoc_thisSet))[1],
	    statistic = this_vegas_result[["statistic"]],
	    simulatedStats = this_vegas_result[["simulated.stats"]]) )
}