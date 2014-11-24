setMethod("geneSetName", "GeneSetResult",
function(object){object@geneSetName}
)
setMethod("geneSetName", "GeneSetResultCollection",
function(object){lapply(object, function(y) y@geneSetName)}
)

setMethod("pValue", "GeneSetResult",
function(object){object@pValue}
)
setMethod("pValue", "GeneSetResultCollection",
function(object){lapply(object, function(y) y@pValue)}
)

setMethod("degreesOfFreedom", "GeneSetResult",
function(object){object@degreesOfFreedom}
)
setMethod("degreesOfFreedom", "GeneSetResultCollection",
function(object){lapply(object, function(y) y@degreesOfFreedom)}
)

setMethod("statistic", "GeneSetResult",
function(object){object@statistic}
)
setMethod("statistic", "GeneSetResultCollection",
function(object){lapply(object, function(y) y@statistic)}
)