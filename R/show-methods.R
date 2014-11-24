## show
setMethod("show", "VEGASResult",
function(object){
    cat(paste0( "VEGAS results for ", object@geneSetName, "\n"))
    cat(paste0("p-value = ", signif(object@pValue,3), "\n") )
    cat(paste0("observed statistic = ", 
        signif(object@statistic,3), "\n") )
    cat(paste0("mean of simulated statistics = ", 
        signif(mean(object@simulatedStats),3), "\n"))
    cat(paste0("variance of simulated statistics = ", 
        signif(var(object@simulatedStats),3), "\n"))
}
)

setMethod("show", "GLOSSIResult",
function(object){
    cat(paste0( "GLOSSI results for ", object@geneSetName, "\n"))
    cat(paste0("p-value = ", signif(object@pValue,3), "\n"))
    cat(paste0("observed statistic = ", 
        signif(object@statistic,3), "\n"))
    cat(paste0("degrees of freedom = ", object@degreesOfFreedom), "\n")
}
)

setMethod("show", "GeneSetResult",
function(object){
    cat(paste0( "Gene Set results for ", object@geneSetName, "\n"))
    cat(paste0("p-value = ", signif(object@pValue,3), "\n"))
    cat(paste0("observed statistic = ", 
        signif(object@statistic,3), "\n"))
    cat(paste0("degrees of freedom = ", object@degreesOfFreedom), "\n")
}
)