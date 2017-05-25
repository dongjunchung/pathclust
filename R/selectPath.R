#' @title pathway selection using LASSO
#'
#' @description  LASSO step for pathway selection. To identify patient subgroups based on pathway level risk profile, we fit a LASSO-Penalized Cox regression on scores derived from all the pathways.
#'
#' @param object results of selectGene step
#' @param seed random seed that was set,  default = 123.
#'
#' @import stats glmnet foreach methods
#' @export
#' @docType methods
#' @rdname selectPath-methods
#' @aliases selectPath
#' @aliases selectPath,FitGene-method
#'
#' @examples
#' data(TCGA)
#' prefilter.results=prefilter( data=TCGA$geneexpr, time=TCGA$t, status=TCGA$d, plist=TCGA$pathList )
#' gene.results=selectGene( object=prefilter.results, fold=5, K=5, etas=c(0.1,0.9),seed=123)
#' path.results=selectPath( object=gene.results, seed=123)

setMethod(
  f="selectPath",
  signature="FitGene",
  definition=function( object, seed=123){

    time<-object@inputdata$time
    status<-object@inputdata$status
    pathways<-object@inputdata$pathway
    score<-object@score

    xx<-foreach(j=1:length(score), .combine='cbind')%do%{score[[j]]}

    set.seed(seed)
    cvlas<-cv.glmnet( x=as.matrix(xx),y=Surv(time, status),
                     family="cox",type.measure="deviance",alpha=1)

    cols<-sapply(score,ncol)
    ##this lambda should also be optional
    path.beta<-as.vector(coef(cvlas,s=cvlas$lambda.1se))
    path.beta<-data.frame( rep(pathways,cols), path.beta,
      stringsAsFactors=FALSE )

    methods::new( "FitPath",
      pathAll=pathways,
      pathSelected=path.beta[which(path.beta[,2]>0),1],
      coef=path.beta,
      fitGene=object
    )
  }
)
