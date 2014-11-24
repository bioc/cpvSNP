setMethod("geneSetName", "VEGASResult",
function(object){object@geneSetName}
)
setMethod("geneSetName", "VEGASResultCollection",
function(object){lapply(object, function(y) y@geneSetName)}
)

setMethod("pValue", "VEGASResult",
function(object){object@pValue}
)
setMethod("pValue", "VEGASResultCollection",
function(object){lapply(object, function(y) y@pValue)}
)

setMethod("degreesOfFreedom","VEGASResult",
function(object){object@degreesOfFreedom}
)
setMethod("degreesOfFreedom", "VEGASResultCollection",
function(object){lapply(object,function(y) y@degreesOfFreedom)}
)

setMethod("statistic", "VEGASResult",
function(object){object@statistic}
)
setMethod("statistic", "VEGASResultCollection",
function(object){lapply(object, function(y) y@statistic)}
)

setMethod("simulatedStats", "VEGASResult",
function(object){object@simulatedStats}
)
setMethod("simulatedStats", "VEGASResultCollection",
function(object){lapply(object, function(y) y@simulatedStats)}
)