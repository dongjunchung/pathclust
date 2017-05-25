
#' @import stats graphics ggplot2 survivalROC
#' @exportMethod coef predict plot

#' @title Summary of selectPath results
#'
#' @description Summary of selectPath results.
#'
#' @aliases show,FitPath-method
#'
#' @param object output of selectPath function

setMethod(
  f="show",
  signature="FitPath",
  definition=function( object ) {

    path.before <- object@pathAll
    path.after <- unique(object@pathSelected)

    cat( "Summary: Pathway-level analysis results (class: FitPath)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Number of all pathways: ", length(path.before), "\n", sep="" )
    cat( "Number of selected pathways: ", length(path.after), "\n", sep="" )
    cat( "\n" )
    cat( "List of selected pathways:\n" )
		for ( i in 1:length(path.after) ) {
			cat( "\t",path.after[i],":\n", sep="" )
		}
    cat( "--------------------------------------------------\n" )
  }
)

#' @title LASSO coefficient estimates for pathways
#'
#' @description LASSO coefficient estimates for pathways.
#'
#' @aliases coef,FitPath-method
#'
#' @param object output of selectPath function

setMethod(
  f="coef",
  signature="FitPath",
  definition=function( object ) {
    return(object@coef)
  }
)

#' @title Risk group prediction
#'
#' @description Risk group prediction.
#'
#' @rdname FitPath-class
#' @aliases predict,FitPath-method
#'
#' @param object output of selectPath function
#' @param newx test gene expression data
#' @param cuts cut points to determine risk groups

setMethod(
  f="predict",
  signature="FitPath",
  definition=function( object, newx=NULL, cuts=NULL){

    ##number of pathways
    if ( is.null(newx) ) {
      newxlist <- object@fitGene@dataPrefiltered
    } else {
      newxlist <- lapply( object@fitGene@dataPrefiltered, function(x){
        idx<- which( (colnames(newx)%in%colnames(x))==T )
        newx[,idx]
      } )
    }
    beta<-object@coef
    w<-object@fitGene@fit$direction
    genes<-object@fitGene@geneSelected

    pathway<-unique(beta[,1])
    np<-length(pathway)
    psrp<-matrix(0,nrow=nrow(newxlist[[1]]),ncol=np)

    for(j in 1:np){
      betaj<-x<-NULL
      idx<-which(beta[,1]==pathway[j])
      betaj<-beta[idx,2]
      x<-as.matrix(newxlist[[j]])
      x<-scale(x,T,T)

      ##predict score
      A<-which(colnames(x)%in%genes[[j]]==T)
      xA<-as.matrix(x[,A])

      if(is.na(w[j])==T){sc<-xA}
      if(is.na(w[j])==F){sc<-xA %*% w[[j]]}

      if(ncol(sc)>1){psrp[,j]<-sc %*% betaj}
      if(ncol(sc)==1){psrp[,j]<-sc*betaj}
    }

    risk<-apply(psrp,1,function(x){length(which(x>0))})

    if(is.null(cuts)){
      low<-as.numeric(summary(risk)[2])
      high<-as.numeric(summary(risk)[5])
    }
    if(length(cuts)==2){
      low<-cuts[1]
      high<-cuts[2]
    }
    riskcat<-ifelse(risk<=low,"low","med")
    riskcat<-ifelse(risk>=high,"high",riskcat)

    return(list(
      risk.index = risk,
      riskcat = riskcat,
      cuts=c(low,high),
      time = object@fitGene@inputdata$time,
      status = object@fitGene@inputdata$status
    ))
  }
)

#' @title Plotting function
#'
#' @description Plot Kaplan-Meier curves (\code{type="KM"}), ROC (\code{type="ROC"}), or hazard ratio (\code{type="HR"}).
#'
#' @rdname FitPath-class
#' @aliases plot,FitPath-method
#'
#' @param x output of selectPath function
#' @param y (unused argument)
#' @param type type of plot. one of Kaplan-Meier curves (\code{type="KM"}), ROC (\code{type="ROC"}), or hazard ratio (\code{type="HR"}).

setMethod(
  f="plot",
  signature=c("FitPath","missing"),
  definition=function( x, y, type="KM" ) {
    if ( type == "KM" ) {

      # plot Kaplan-Meier curves

      predicted <- predict(x)

      time<-predicted$time
      status<-predicted$status
      reg<-survfit(Surv(time,status)~factor(predicted$riskcat))
      cuts<-predicted$cuts
      tab<-table(predicted$riskcat)
      plot( reg, col=c("red","black","green"),xlim=c(0,60),
            xlab="Months From Diagnosis to Death",ylab="Survival Fraction (KM)" )
      legend("bottomleft",lty=c(1,1,1),col=c("red","black","green"),
             c( paste("High",">",cuts[2],"n=",tab[1]),paste("Medium", "n=",tab[2]),paste("Low","<",cuts[1], "n=",tab[3])))

    } else if ( type == "ROC" ) {

      # plot time-dependent ROC

      predicted <- predict(x)

      time<-predicted$time
      status<-predicted$status
      roc<-survivalROC::survivalROC(Stime=time,status=status,
        marker<-predicted$risk.index, predict.time=60, method="KM")


      plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="red",
           xlab="1-Specificity", ylab="Sensitivity")
      text(x=0.8, y=0.1, labels = paste("AUC =", round((roc$AUC),3)))
      abline(0,1)

    } else if ( type == "HR" ) {

      # plot hazard ratios of selected pathways

      path.results<-x@coef
      betas<-path.results[,2]
      paths<-path.results[,1]
      index<-coefs<-list()
      mat<-path.results[betas!=0,]
      lab2=paths
      #lab1=gsub( "KEGG_", "", paths )
      #lab2=gsub( "_PATHWAY", "", lab1 )

      for(i in 1:length(paths)){
        j<-which( mat[,1]==paths[i] )
        coefs[[i]]<-mat[j,2]
      }

      names(coefs)<-lab2
      o<-rev(order(-sapply(coefs,max)))
      coefs<-coefs[o]
      lab2<-names(coefs)

      coefs<-unique(coefs)
      lab2<- unique(lab2)
      for(i in 1:length(coefs)){
        if(length(coefs[[i]])>1){
          index[[i]]<-rev(seq( i-0.1*length(coefs[[i]])+0.1, i , 0.1 ))
        }
        if(length(coefs[[i]])==1){index[[i]]<-i}
      }

      index.plot<- unlist(index)
      coefs.plot<- exp(unlist(coefs))
      path.plot<-lab2
      data.plot<- as.data.frame(cbind(index.plot,coefs.plot))

      ggplot(data.plot,aes(x=index.plot,xend=index.plot,y=1,yend=coefs.plot))+
        geom_segment()+coord_flip()+
        scale_x_discrete(limits=path.plot)+
        ggtitle("Hazard ratio") +
        theme_bw()+scale_y_continuous(expand = c(0, 0),limits= c(1 , max(data.plot$coefs.plot)+0.05))+
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              plot.title =element_text(size=15, vjust=1,hjust=0.5))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))

    } else {
      stop("Inappropriate 'type' argument!")
    }
  }
)

