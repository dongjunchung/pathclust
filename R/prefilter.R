#' @title prefilter function
#'
#' @description a supervised prefiltering by fitting a cox regression model on the expression measures of each mRNA.  Only genes with p values smaller than a pre-selected cut-off point were included in further analysis. Here we chose p=0.5 as default cut-off point.
#'
#' @export
#' @import foreach survival methods
#' @param data A N by P dataframe, where N is number of observations, P is number of genes
#' @param time time passed to selectGene's time argument.
#' @param status status selectGene's status argument.
#' @param p.cut p.cut is the pvalue cutoff point for prefiltering
#' @param plist user provided gene lists organized in pathways
#'
#'
#' @examples
#' data(TCGA)
#' prefilter.results=prefilter( data=TCGA$geneexpr, time=TCGA$t, status=TCGA$d, plist=TCGA$pathList )


prefilter <- function( data, time, status, p.cut=0.5, plist=plist ){
  pvals <- foreach(i=1:ncol(data), .combine='c') %do% {
    summary(coxph( Surv(time, status)~data[,i] ))$coef[5]
  }

  xlist<- lapply( plist,function(x){
    idx<- which( (names(data)%in%x)==T & (pvals<=p.cut)==T )
    data[,idx]
  } )

  methods::new( "Prefiltered",
    xlist = xlist,
    inputdata = list( data = data, time = time, status = status, pathway = names(xlist) )
  )
}
