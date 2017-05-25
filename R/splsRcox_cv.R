#' @import rms plsRcox

cv.coxsplsDR2=
function (data, method = c("efron", "breslow"), nfold = 5, nt = 10,
    eta = 0.5, se = TRUE, givefold, scaleX = TRUE, scaleY = FALSE )
{
    cv.error10 <- NULL
    cv.se10 <- NULL
    pred10<-NULL

    method <- match.arg(method)
    x <- data$x
    time <- data$time
    status <- data$status
    n <- length(time)

    if (missing(givefold)) {
        folds <- split(sample(seq(n)), rep(1:nfold, length = n))
    }else {
        folds <- givefold
    }

    show_nbr_var <-TRUE
    errormat10<-matrix(NA, nt+1, nfold)

    for (i in seq(nfold)) {
        pred10=rep(NA, nt+1)
        o <- folds[[i]]
        tr <- list(x=x[-o, ], time=time[-o], status=status[-o])
        ts <- list( x=x[o, ], time=time[o], status=status[o])

        cox<-coxph(Surv(tr$time, tr$status) ~ 1)
        devres<-residuals(cox,type="deviance")
     	spls.mod<-spls.cox(x=tr$x, y=devres, K=nt, eta=eta,
       kappa=0.5, select="pls2", scale.x=T, scale.y=F, trace=F)

		sci<-data.frame(spls.mod$plsmod$variates$X)
		wi<-spls.mod$pred$w
		coeffit<-matrix(NA,nrow=nt,ncol=nt)
		newAs<-list()
		for(s in 1:nt){
			coxfit=coxph(Surv(tr$time, tr$status)~as.matrix(sci)[,1:s])
		    if(s==1){coeffit[1:s,s]=summary(coxfit)$coef[1]}
		    if(s>1){coeffit[1:s,s]=summary(coxfit)$coef[,1]}
		}

        nzb <- cumsum(c(0, sapply(spls.mod$new2As, length)))

        for (jj in 1:nt) {
            Aval <- spls.mod$A
            newxdata <- predict.pls.cox( spls.mod$plsmod, newdata=scale(ts$x[,Aval], 				spls.mod$meanx[Aval], spls.mod$normx[Aval]),scale.X=F, scale.Y=F)$variates[,1:jj, drop=F]

            oldxdata <- sci[, 1:jj, drop=F]
            allxdata <- predict.pls.cox( spls.mod$plsmod, newdata=scale(x[, Aval], 						spls.mod$meanx[Aval], spls.mod$normx[Aval]),scale.X=F, scale.Y=F)$variates[, 1:jj, drop=F]

            if (jj == 1) {pred10[1]=0.5}

            predict.trainvectjj <- as.matrix(oldxdata) %*% coeffit[1:jj, jj, drop=F]
            predictvectjj <- as.matrix(newxdata) %*% coeffit[1:jj, jj, drop=F]
            Xlp <- rep(NA, length(time))
            Xlp[-o] <- predict.trainvectjj
            Xlp[o] <- predictvectjj
            tmp=as.data.frame(cbind(time=time, status=status, Xlp=Xlp))
            TR <- tmp[-o, ]
            TE <- tmp[o, ]

            survival.time <- tmp[, "time"]
            survival.status <- tmp[, "status"]

            Surv.rsp <- Surv(tr$time, tr$status)
            Surv.rsp.new <- Surv(ts$time, ts$status)
            train.fit <- coxph(Surv(time, status) ~ Xlp, x = TRUE,
                y = TRUE, method = method, data = TR, iter.max = 0,
                init = 1)
            train.fit.cph <- rms::cph(Surv(time, status) ~ Xlp, x = TRUE,
                y = TRUE, method = method, data = TR, iter.max = 0,
                init = 1)
            lp <- predict(train.fit)
            lpnew <- predict(train.fit, newdata = TE)

            AUCs <- getIndicCViAUCSurvROCTest(lp, lpnew,
                  Surv.rsp, Surv.rsp.new, times.auc = seq(0,
                    max(time), length.out = 1000), times.prederr = seq(0,
                    max(time), length.out = 1000)[-(990:1000)],
                  train.fit, plot.it=F)

            pred10[jj + 1] = AUCs$AUC_survivalROC.test$iauc
        }
        if (is.na(pred10[1])) {
            pred10[1] <- 0.5
        }
        if (length(o) == 1) {
            for (ind in 1:number_ind) {
                assign(paste("pred", ind, sep = ""), matrix(get(paste("pred",
                  ind, sep = "")), nrow = 1))
            }
        }

        errormat10[, i] <- ifelse(is.finite(pred10), pred10, NA)
    }

        cv.error10<- apply(errormat10, 1, mean, na.rm = TRUE)
        cv.se10<-sqrt(apply(errormat10, 1, var, na.rm = TRUE))/nfold

        object <- list(nt = nt, cv.error10 = cv.error10, cv.se10 = cv.se10,
            folds = folds, nzb = nzb)
}
