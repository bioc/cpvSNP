
# unit tests for the cpvSNP package

test_createArrayData <- function(){
	data(geneSetAnalysis)
	arrayData <- geneSetAnalysis[["arrayData"]]
	checkTrue(class(createArrayData(arrayData,positionName="Position"))=="GRanges")
}

test_geneToSNPList <- function(){
	data(geneSetAnalysis)
	data(test_genesHg19)
	
	geneList <- geneSetAnalysis[["geneSets"]]
	arrayData <- geneSetAnalysis[["arrayData"]]
	arrayDataGR <- createArrayData(geneSetAnalysis[["arrayData"]],
	    position="Position", verbose=FALSE)

	checkException(geneToSNPList(geneList,arrayDataGR,test_genesHg19,maxgap=NA))
	checkException(geneToSNPList(geneList,arrayDataGR,test_genesHg19,maxgap=TRUE))
	checkException(geneToSNPList(geneList,arrayDataGR,test_genesHg19,maxgap="yes"))	
	
	checkException(geneToSNPList(geneList,arrayDataGR,test_genesHg19,quiet=NA))
	checkException(geneToSNPList(geneList,arrayDataGR,test_genesHg19,quiet=12))
	checkException(geneToSNPList(geneList,arrayDataGR,test_genesHg19,quiet="yes"))	
		
	checkTrue(class(geneToSNPList(geneList,arrayDataGR,test_genesHg19))=="GeneSetCollection")
}

test_glossi <- function(){
	data(geneSetAnalysis)
    data(test_snpsGSC)

	pvals <- geneSetAnalysis[["arrayData"]]$P
	names(pvals) <- geneSetAnalysis[["arrayData"]]$SNP
	
	gRes <- glossi(pvals, snpsGSC[[1]])
	checkTrue(class(pValue(gRes))=="numeric")
	checkTrue(class(degreesOfFreedom(gRes))=="integer")
	checkTrue(class(statistic(gRes))=="numeric")
	
	checkEqualsNumeric(degreesOfFreedom(gRes),607)
	
	gRes <- glossi(pvals, snpsGSC)
	checkTrue(length(pValue(gRes))==2)
	checkTrue(length(degreesOfFreedom(gRes))==2)
	checkTrue(length(statistic(gRes))==2)
}

test_vegas <- function(){
	data(geneSetAnalysis)
    data(test_snpsGSC)
	
	arrayDataGR <- createArrayData(geneSetAnalysis[["arrayData"]],
	    position="Position", verbose=FALSE)
	vRes <- vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]])
	
	# check output
	checkTrue(class(vRes)=="VEGASResult")
	checkTrue(class(pValue(vRes))=="numeric")
	checkTrue(class(statistic(vRes))=="numeric")
	checkTrue(class(degreesOfFreedom(vRes))=="integer")
	checkTrue(length(simulatedStats(vRes))==1000)	
	
	# checking input errors/exceptions
	
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], num_sims=0))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], num_sims=NA))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]],
	    num_sims=TRUE))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]],
	    num_sims="yes"))

	checkException(vegas(NA, arrayDataGR, geneSetAnalysis[["ldMat"]]))
	checkException(vegas("yes", arrayDataGR, geneSetAnalysis[["ldMat"]]))
	checkException(vegas(12, arrayDataGR, geneSetAnalysis[["ldMat"]]))
	checkException(vegas(TRUE, arrayDataGR, geneSetAnalysis[["ldMat"]]))
	checkException(vegas(geneSet(c(0,0)), arrayDataGR, geneSetAnalysis[["ldMat"]]))
 
	checkException(vegas(snpsGSC[[1]], "yes", geneSetAnalysis[["ldMat"]]))
	checkException(vegas(snpsGSC[[1]], 12, geneSetAnalysis[["ldMat"]]))
	checkException(vegas(snpsGSC[[1]], TRUE, geneSetAnalysis[["ldMat"]]))
	checkException(vegas(snpsGSC[[1]], GRanges(c(0,0)), geneSetAnalysis[["ldMat"]]))
	
	checkException(vegas(snpsGSC[[1]], arrayDataGR, 1))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, NA))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, TRUE))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, "yes"))	
		
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], 
	    correction=NA))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], 
	    correction="yes"))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], 
	    correction=12))
	
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], 
	    seed=NA))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], 
	    seed="yes"))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], 
	    seed=TRUE))

	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], 
	    verbose=NA))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], 
	    verbose=12))
	checkException(vegas(snpsGSC[[1]], arrayDataGR, geneSetAnalysis[["ldMat"]], 
	    verbose="yes"))
}
	
	
	