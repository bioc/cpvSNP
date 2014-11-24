setClass("GeneSetResult",
    representation(
        geneSetName = "character",
        pValue = "numeric",
        degreesOfFreedom = "integer",
        statistic = "numeric")
)
setClass("GeneSetResultCollection",
    contains="list",
    validity = function(object) {
        msg <- NULL
        if ( !all(sapply(object, is, "GeneSetResult")) )
            msg <- c(msg, "members must all be 'GeneSetResult' classes")
        if ( !is.null(msg) )
            msg
        else
            TRUE
    }
)

setClass("VEGASResult",
	contains="GeneSetResult",
    representation( 
        simulatedStats = "vector")
)

setClass("VEGASResultCollection",
    contains="list",
    validity = function(object) {
        msg <- NULL
        if ( !all(sapply(object, is, "VEGASResult")) )
            msg <- c(msg, "members must all be 'VEGASResult' classes")
        if ( !is.null(msg) )
            msg
        else
            TRUE
    }
)

setClass("GLOSSIResult",
    contains="GeneSetResult"
)

setClass("GLOSSIResultCollection",
    contains="list",
    validity = function(object) {
        msg <- NULL
        if ( !all(sapply(object, is, "GLOSSIResult")) )
            msg <- c(msg, "members must all be 'GLOSSIResult' classes")
        if ( !is.null(msg) )
            msg
        else
            TRUE
    }
)