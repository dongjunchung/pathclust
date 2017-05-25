#' @title Gene selection using SPLS model
#'
#' @description Sparse Partial Least Square (SPLS) step for gene selection and dimension reduction.By applying SPLS to each pathway, we achieve the goal of gene selection and dimension reduction at the same time.
#'
#' @param object output list of prefilter step
#' @param fold The number of folds to use to perform the cross-validation process.
#' @param K the maximum number of hidden features in spls.
#' @param etas Thresholding parameter. eta should be between 0 and 1.
#' @param seed random seed that was set,  default = 123.
#'
#' @import stats foreach plsRcox survival methods
#' @export
#' @docType methods
#' @rdname selectGene-methods
#' @aliases selectGene
#' @aliases selectGene,Prefiltered-method
#'
#' @examples
#' data(TCGA)
#' prefilter.results=prefilter( data=TCGA$geneexpr, time=TCGA$t, status=TCGA$d, plist=TCGA$pathList )
#' gene.results=selectGene( object=prefilter.results, fold=5, K=5, etas=c(0.1,0.9),seed=123)

setMethod(
  f="selectGene",
  signature="Prefiltered",
  definition=function( object, fold=5, K=5, etas=seq(0.1,0.9,0.1), seed=123 ) {

    time<-object@inputdata$time
    status<-object@inputdata$status
    data<-object@xlist
    pathways<-object@inputdata$pathway

    set.seed(seed)
    n<-length(time)
    cvfolds <- split(sample(n), rep(1:fold, length=n))
    dimx<-unlist( lapply(data,function(x){ncol(as.matrix(x))}) )

    k.opt<-eta.opt<-NULL
    score<-genes<-beta<-spls.beta<-w<-list()

    for(j in 1:length(pathways)){
      xx<-as.matrix( data[[j]],nrow=n,ncol=dimx[j] )
      kmax<-min( K, ncol(xx) )

      if(kmax>1){
        aucs <- foreach(i=1:length(etas),.combine='rbind') %do% {
          suppressWarnings(cvi <- cv.coxsplsDR2(
            data=list(x=xx,time=time,status=status),
            nfold=fold, nt=kmax, eta=etas[i],
            se=TRUE, givefold=cvfolds,
            scaleX=TRUE, scaleY=FALSE ))
            #plot.it=F, sclaleY=F ))

          cbind( cvi$cv.error10, cvi$cv.se10	)
        }

        h<-which.max(aucs[,1])
        se1<-aucs[h,1]-aucs[h,2]

        ###find aucs within 1se of the maximum
        mat<-cbind( aucs, rep(0:kmax,length(etas)),
                   rep(etas,each=(kmax+1)) )[aucs[,1]>se1,]

        tmp<-aucs[,1][aucs[,1]>se1]

        if(length(tmp)>1){
          ##choose the most parsimonious model in terms of genes
          ##always choose the smallest k among those with largest eta
          q<-which.max(mat[,4])
          k.opt[j] <- mat[q,3]
          eta.opt[j] <- mat[q,4]
        }

        if(length(tmp)==1){
          k.opt[j]<-mat[3]
          eta.opt[j]<-mat[4]
        }

        cox<-coxph(Surv(time, status) ~ 1)
        devres<-residuals(cox,type="deviance")

        spls.mod<-spls.cox (x=xx, y=devres, K=k.opt[j], eta=eta.opt[j],
                           kappa=0.5, select="pls2", scale.x=T, scale.y=F, trace=F)

        score[[j]]<- data.frame(spls.mod$plsmod$variates$X)
        spls.beta[[j]]<-data.frame( colnames(xx),spls.mod$betahat )
        rownames(spls.beta[[j]])<-NULL

        ##objects saved for prediction
        xA<-spls.mod$x[,spls.mod$A]
        genes[[j]]<-colnames(xx)[spls.mod$A]
        w[[j]]<-spls.mod$pred$w

      }else{

        score[[j]]<-xx
        genes[[j]]<-names(xx)
        spls.beta[[j]]<-NA
        k.opt[j]<-1
        eta.opt[j]<-w[[j]]<-NA
      }
    }

    names(genes)<-pathways
    names(spls.beta)<-pathways

    methods::new( "FitGene",
      score=score,
      geneSelected=genes,
      fit = list( coef=spls.beta, direction=w ),
      dataPrefiltered=data,
      inputdata = object@inputdata
    )
  }
)
