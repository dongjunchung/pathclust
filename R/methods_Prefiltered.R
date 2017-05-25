#' @title Prefiltering summary
#'
#' @description Summary of prefiltering results.
#'
#' @aliases show,Prefiltered-method
#'
#' @param object output of prefilter function

setMethod(
  f="show",
  signature="Prefiltered",
  definition=function( object ) {

    gene.before <- colnames(object@inputdata$data)
    gene.after <- unique(unlist(sapply( object@xlist, colnames )))

    cat( "Summary: Pre-filtering results (class: Prefiltered)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Number of genes before prefiltering: ", length(gene.before), "\n", sep="" )
    cat( "Number of genes after prefiltering: ", length(gene.after), "\n", sep="" )
    cat( "--------------------------------------------------\n" )
  }
)
