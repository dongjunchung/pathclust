#' @title Summary of selectGene results
#'
#' @description Summary of selectGene results.
#'
#' @aliases show,FitGene-method
#'
#' @param object output of selectGene function

setMethod(
  f="show",
  signature="FitGene",
  definition=function( object ) {

    gene.before <- unique(unlist(sapply( object@dataPrefiltered, colnames )))
    gene.after <- unique(unlist(object@geneSelected))

    cat( "Summary: Gene-level analysis results (class: FitGene)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Number of prefiltered genes: ", length(gene.before), "\n", sep="" )
    cat( "Number of selected genes: ", length(gene.after), "\n", sep="" )
    cat( "--------------------------------------------------\n" )
  }
)

#' @title SPLS coefficient estimates for genes
#'
#' @description SPLS coefficient estimates for genes.
#'
#' @aliases coef,FitGene-method
#'
#' @param object output of selectGene function

setMethod(
  f="coef",
  signature="FitGene",
  definition=function( object ) {
    return(object@fit$coef)
  }
)
