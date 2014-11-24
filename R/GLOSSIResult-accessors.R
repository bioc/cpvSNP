setMethod("geneSetName", "GLOSSIResult",
function(object){object@geneSetName}
)
setMethod("geneSetName", "GLOSSIResultCollection",
function(object){lapply(object, function(y) y@geneSetName)}
)

setMethod("pValue", "GLOSSIResult",
function(object){object@pValue}
)
setMethod("pValue", "GLOSSIResultCollection",
function(object){lapply(object, function(y) y@pValue)}
)

setMethod("degreesOfFreedom", "GLOSSIResult",
function(object){object@degreesOfFreedom}
)
setMethod("degreesOfFreedom", "GLOSSIResultCollection",
function(object){lapply(object, function(y) y@degreesOfFreedom)}
)

setMethod("statistic", "GLOSSIResult",
function(object){object@statistic}
)
setMethod("statistic", "GLOSSIResultCollection",
function(object){lapply(object, function(y) y@statistic)}
)