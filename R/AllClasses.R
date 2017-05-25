
#' An S4 class to represent pre-filtered data.
#'
#' @slot xlist list of genes in each pathway
#' @slot inputdata input data for prefilter function

setClass( Class="Prefiltered",
    representation=representation(
        xlist="list",
        inputdata="list"
    )
)

#' An S4 class to represent gene-level analysis results
#'
#' @slot score list of SPLS latent components for each pathway
#' @slot geneSelected list of selected genes in each pathway
#' @slot fit list of SPLS fits
#' @slot dataPrefiltered list of prefiltered data for each pathway
#' @slot inputdata input data for prefilter function
#'
setClass( Class="FitGene",
    representation=representation(
        score="list",
        geneSelected="list",
        fit="list",
        dataPrefiltered="list",
        inputdata="list"
    )
)

#' An S4 class for pathway-level analysis results
#'
#' @slot PathAll all pathways
#' @slot PathSelected pathways selected by LASSO
#' @slot coef LASSO coefficient estimates for pathway latent components
#' @slot fitGene gene-level analysis results

setClass( Class="FitPath",
    representation=representation(
        pathAll="character",
        pathSelected="character",
        coef="data.frame",
        fitGene="FitGene"
    )
)
