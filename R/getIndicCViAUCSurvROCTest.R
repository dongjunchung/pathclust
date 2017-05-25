#' @import survcomp survival

getIndicCViAUCSurvROCTest = function(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),times.prederr=1:500,train.fit,plot.it=FALSE,tmax.train=max(Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ]),tmax.test=max(Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ])){
  #try(attachNamespace("survival"),silent=TRUE)
  #on.exit(try(unloadNamespace("survival"),silent=TRUE))
  #try(attachNamespace("survcomp"),silent=TRUE)
  #on.exit(try(unloadNamespace("survcomp"),silent=TRUE))
  object <- NULL

  object$lp <- lp
  object$lpnew <- lpnew
  object$Surv.rsp <- Surv.rsp
  object$Surv.rsp.new <- Surv.rsp.new
  object$times.auc <- times.auc
  object$train.fit <- train.fit
  object$test.fit <- survival::coxph(object$Surv.rsp.new~object$lpnew, iter.max=0, init=1)

  object$nulltrain.fit <- survival::coxph(object$Surv.rsp~1)
  object$lp0 <- predict(object$nulltrain.fit)
  object$nulltest.fit <- survival::coxph(object$Surv.rsp.new~1)
  object$lp0new <- predict(object$nulltest.fit)

  object$tmax.train <- tmax.train
  object$tmax.test <- tmax.test


  object$utimes.test <- unique( object$Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ] )
  object$utimes.test <- object$utimes.test[ order(object$utimes.test) ]

  #  library(survcomp)
  ##test
  mytdroc.test <- NULL
  object$AUC_survivalROC.test <- NULL
  object$AUC_survivalROC.test$auc <- rep(NA,length(object$utimes.test))
  object$AUC_survivalROC.test$iauc <- NA
  object$AUC_survivalROC.test$times <- object$utimes.test
  class(object$AUC_survivalROC.test) <- "survAUC"
  test.cc.ix <- complete.cases(object$lpnew, object$Surv.rsp.new[,"time"], object$Surv.rsp.new[,"status"], NULL)
  test.surv.event.cc.ix <- object$Surv.rsp.new[,"status"][test.cc.ix]
  if (all(sort(unique(test.surv.event.cc.ix)) == c(0, 1))) {
    for(i in 1:length(object$utimes.test)) {
      rr.test <- survcomp::tdrocc(x=object$lpnew, surv.time=object$Surv.rsp.new[,"time"], surv.event=object$Surv.rsp.new[,"status"], time=object$utimes.test[i], na.rm=TRUE, verbose=FALSE)
      mytdroc.test <- c(mytdroc.test, list(rr.test))
    }
    object$AUC_survivalROC.test$auc <- unlist(lapply(mytdroc.test, function(x) { return(x$AUC) }))
    cc.ix.test <- complete.cases(object$AUC_survivalROC.test$auc)
    auc.survivalROC.test.cc <- object$AUC_survivalROC.test$auc[cc.ix.test]
    time.test.cc <- object$utimes.test[cc.ix.test]
    if(length(time.test.cc)>0){
      diffs.test.cc <- c(time.test.cc[1], time.test.cc[2:length(time.test.cc)] - time.test.cc[1:(length(time.test.cc) - 1)])
      object$AUC_survivalROC.test$iauc <- sum(diffs.test.cc * auc.survivalROC.test.cc)/max(time.test.cc)
      if(plot.it){
        plot(object$AUC_survivalROC.test)
        abline(h = 0.5)
      }
    }
  }

  return(object)
}
