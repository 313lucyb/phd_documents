library(foreign)
library(rstpm2)
library(caret)
library(mfp)
library(rms)
library(ROCR)
library(survival)
library(Hmisc)
library(MASS)
library(reshape2)
library(flexsurv)
library(pROC)
library(DescTools)
library(boot)
library(randomForestSRC)
library(party)
library(ggRandomForests)
library(risksetROC)
library(akima)
library(ipred)
library(pec)
library(colorRamps)
library(grDevices)
library(RColorBrewer)
cols <-brewer.pal(11, "RdBu")

## implementation challenge calibration pre-post optimisation

SRF <- read.dta("baseline_model2.dta")
SRF$smoke<-factor(SRF$smoke)
SRF$firsttreat<-factor(SRF$firsttreat)
SRF$pgen<-factor(SRF$pgen)
SRF$renal<-factor(SRF$renal)
SRF$lung<-factor(SRF$lung)
SRF$diabetes<-factor(SRF$diabetes)
SRF$steroids<-factor(SRF$steroids)

set.seed(16645)
si.obj.before <- rfsrc(Surv(timevent, event) ~ age + pgen + firsttreat + 
                    disdur + previous_dmards + lung + 
                    diabetes + bmi + renal + steroids + ovmean + 
                    smoke + dascore, data = SRF, tree.err=TRUE, 
                  ntree = 500, importance = TRUE)
test0<-SRF
test0$pred0<-predictSurvProb(object = si.obj.before, newdata=test0, 
                             times = 1)

set.seed(76783)
val_ests0 <- val.surv(est.surv = as.numeric(test0$pred0), 
                      S = Surv(test0$timevent, test0$event),
                      u=1, fun=function(p)log(-log(p)),
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests0,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
groupkm(x = as.numeric(test0$pred0), S = Surv(test0$timevent,test0$event), 
        g=10, u=1, pl=T, add=T, lty=0, lwd=1.5, col=cols[10], 
        cex.subtitle=FALSE, cex=1.2, pch=21, bg=cols[10])
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[10], "black","black"), cex = 1.1)

min(test0$pred0)
# # CPM Development

SRF <- read.dta("baseline_model2.dta")
SRF$smoke<-factor(SRF$smoke)
SRF$firsttreat<-factor(SRF$firsttreat)
SRF$pgen<-factor(SRF$pgen)
SRF$renal<-factor(SRF$renal)
SRF$lung<-factor(SRF$lung)
SRF$diabetes<-factor(SRF$diabetes)
SRF$steroids<-factor(SRF$steroids)

set.seed(56655)

# Hyperparameter optimisation

rsf.tune<-tune(formula = Surv(timevent, event) ~ age + pgen + firsttreat + 
                 disdur + previous_dmards + lung + 
                 diabetes + bmi + renal + steroids + ovmean + 
                 smoke + dascore, data = SRF, 
               ntreeTry = 500, nodesizeTry = seq(10, 300, by = 10), 
               nsplit = 5, trace = TRUE, doBest = TRUE)

plot.tune <- function(o, linear = TRUE) {
  x <- o$results[,1]
  y <- o$results[,2]
  z <- o$results[,3]
  so <- interp(x=x, y=y, z=z, linear = linear)
  #data <- data.frame(x=x, y=y, z=z)
  #d <- ggplot(data, aes(x, y, z=z))
  #d + geom_tile(aes(fill=z))  + scale_fill_gradient2(low="blue", high="red")
  
  idx <- which.min(z)
  x0 <- x[idx]
  y0 <- y[idx]
  ++
  filled.contour(x = so$x,
                 y = so$y,
                 z = so$z,
                 xlim = range(so$x, finite = TRUE) + c(-2, 2),
                 ylim = range(so$y, finite = TRUE) + c(-2, 2),
                 color.palette =
                   colorRampPalette(c("lightblue", "red")),
                 xlab = "Terminal node size",
                 ylab = "No. variables tried at each split",
                 #main = "error rate for nodesize and mtry",
                 key.title = {title(main = "OOB error", cex.main=1, cex.lab=2)},
                 plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1.2,font=2);
                   points(x,y,pch=16,cex=.25)})
}
plot.tune(rsf.tune)

# Optimised ntry/nodesize: 3, 280
# 300 trees 
start_time <- Sys.time()
si.obj <- rfsrc(Surv(timevent, event) ~ age + pgen + firsttreat + 
                  disdur + previous_dmards + lung + 
                  diabetes + bmi + renal + steroids + ovmean + 
                  smoke + dascore, 
                data = SRF, tree.err=TRUE, ntree = 1000, 
                nodesize = 280, mtry = 3, nsplit = 5,
                importance = TRUE, seed = 334)
end_time <- Sys.time()
fit_time <- end_time - start_time


plot(gg_vimp(si.obj))
plot(si.obj)

set.seed(333)

# 300 tree model

start_time <- Sys.time()
si.obj2 <- rfsrc(Surv(timevent, event) ~ age + pgen + firsttreat + 
                   disdur + previous_dmards + lung + 
                   diabetes + bmi + renal + steroids + ovmean + 
                   smoke + dascore, 
                 data = SRF, mtry = 3, nodesize =280, tree.err=TRUE, 
                 ntree = 300, importance = TRUE, seed = 333)
end_time <- Sys.time()
fit_time <- end_time - start_time
plot(si.obj2)


# Out-of-bag (OOB) predictive accuracy

spred.oob<-si.obj2$survival.oob[,291]

set.seed(76783)
par(pty="s")
val_ests0.oob <- val.surv(est.surv = spred.oob, 
                      S = Surv(SRF$timevent, SRF$event),
                      u=1, fun=function(p)log(-log(p)),
                      pred = sort(runif(100, 0.8, 1)))
plot(val_ests0.oob,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(x = spred.oob, S = Surv(SRF$timevent,SRF$event), 
        g=10,u=1, pl=T, add=T,lty=0, col="plum", 
        cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c("plum", "black", "black"), cex = 1)


###########################################################
##                  Temporal Assessment                  ##
###########################################################

# Performance measure vectors

LT<-c(0,0.5,1,1.5,2,2.5,3)
Cind<-0
CindL<-0
CindU<-0
BS<-0

# # LT = 0 

# Discrimination

#test0<-SRF[sample(nrow(SRF), nrow(SRF), replace = TRUE),]
test0<-SRF
start_time0 <- Sys.time()
test0$pred0<-predictSurvProb(object = si.obj2, newdata=test0, 
                             times = 1)
end_time0 <- Sys.time()
pred_time0 <- end_time0 - start_time0

Cind[1]<- rcorr.cens(test0$pred0, 
                     Surv(test0$timevent, test0$event))["C Index"]

CindL[1]<-rcorr.cens(test0$pred0, Surv(test0$timevent, 
                                 test0$event))["C Index"]-
  1.96*rcorr.cens(test0$pred0, Surv(test0$timevent, 
                              test0$event))["S.D."]
CindU[1]<-rcorr.cens(test0$pred0, Surv(test0$timevent, 
                                 test0$event))["C Index"]+
  1.96*rcorr.cens(test0$pred0, Surv(test0$timevent, 
                              test0$event))["S.D."]

# Calibration (Brier score)

BS[1]<-sbrier(obj = Surv(test0$timevent, test0$event), 
              pred = as.numeric(test0$pred0), btime = 1)

# Calibration plot
par(pty="s")
set.seed(76783)
val_ests0 <- val.surv(est.surv = as.numeric(test0$pred0), 
                      S = Surv(test0$timevent, test0$event),
                      u=1, fun=function(p)log(-log(p)),
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests0,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(x = as.numeric(test0$pred0), S = Surv(test0$timevent,test0$event), 
        g=10, u=1, pl=T, add=T, lty=0, col=cols[4], 
        cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[4], "black", "black"), cex = 1)


min(test0$pred0)

# # LT = 0.5

test1<-read.dta("validation_6months.dta")
test1$smoke<-factor(test1$smoke)
test1$firsttreat<-factor(test1$firsttreat)
test1$pgen<-factor(test1$pgen)
test1$renal<-factor(test1$renal)
test1$lung<-factor(test1$lung)
test1$diabetes<-factor(test1$diabetes)
test1$steroids<-factor(test1$steroids)

start_time1 <- Sys.time()
test1$pred1<-predictSurvProb(object = si.obj2, newdata=test1, 
                             times = 1.5)
end_time1 <- Sys.time()
pred_time1 <- end_time1 - start_time1

Cind[2]<- rcorr.cens(test1$pred1, 
                     Surv(test1$timevent, test1$event))["C Index"]

CindL[2]<-rcorr.cens(test1$pred1, Surv(test1$timevent, 
                                       test1$event))["C Index"]-
  1.96*rcorr.cens(test1$pred1, Surv(test1$timevent, 
                                    test1$event))["S.D."]
CindU[2]<-rcorr.cens(test1$pred1, Surv(test1$timevent, 
                                       test1$event))["C Index"]+
  1.96*rcorr.cens(test1$pred1, Surv(test1$timevent, 
                                    test1$event))["S.D."]

# Calibration (Brier score)

BS[2]<-sbrier(obj = Surv(test1$timevent, test1$event), 
              pred = as.numeric(test1$pred1), btime = 1.5)

# Calibration plot
set.seed(7678)
val_ests1 <- val.surv(est.surv = test1$pred1, 
                      S = Surv(test1$timevent, test1$event),
                      u=1.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.7, 1)))
plot(val_ests1,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test1$pred1, S = Surv(test1$timevent,test1$event), 
        g=10,u=1.5, pl=T, add=T,lty=0, col=cols[4], 
        cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[4], "black", "black"), cex = 0.8)


# # LT = 1

test2<-read.dta("validation_12months.dta")
test2$smoke<-factor(test2$smoke)
test2$firsttreat<-factor(test2$firsttreat)
test2$pgen<-factor(test2$pgen)
test2$renal<-factor(test2$renal)
test2$lung<-factor(test2$lung)
test2$diabetes<-factor(test2$diabetes)
test2$steroids<-factor(test2$steroids)

start_time2 <- Sys.time()
test2$pred2<-predictSurvProb(object = si.obj2, newdata=test2, 
                             times = 2)
end_time2 <- Sys.time()
pred_time2 <- end_time2 - start_time2

Cind[3]<- rcorr.cens(test2$pred2, 
                     Surv(test2$timevent, test2$event))["C Index"]

CindL[3]<-rcorr.cens(test2$pred2, Surv(test2$timevent, 
                                       test2$event))["C Index"]-
  1.96*rcorr.cens(test2$pred2, Surv(test2$timevent, 
                                    test2$event))["S.D."]
CindU[3]<-rcorr.cens(test2$pred2, Surv(test2$timevent, 
                                       test2$event))["C Index"]+
  1.96*rcorr.cens(test2$pred2, Surv(test2$timevent, 
                                    test2$event))["S.D."]

# Calibration (Brier score)

BS[3]<-sbrier(obj = Surv(test2$timevent, test2$event), 
              pred = as.numeric(test2$pred2), btime = 2)

# Calibration plot
set.seed(767832)
val_ests2 <- val.surv(est.surv = test2$pred2, 
                      S = Surv(test2$timevent, test2$event),
                      u=2, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests2,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test2$pred2, S = Surv(test2$timevent,test2$event), 
        g=10,u=2, pl=T, add=T,lty=0, col=cols[4], 
        cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[4], "black", "black"), cex = 0.8)


# # LT = 1.5

test3<-read.dta("validation_18months.dta")
test3$smoke<-factor(test3$smoke)
test3$firsttreat<-factor(test3$firsttreat)
test3$pgen<-factor(test3$pgen)
test3$renal<-factor(test3$renal)
test3$lung<-factor(test3$lung)
test3$diabetes<-factor(test3$diabetes)
test3$steroids<-factor(test3$steroids)

start_time3 <- Sys.time()
test3$pred3<-predictSurvProb(object = si.obj2, newdata=test3, 
                             times = 2.5)
end_time3 <- Sys.time()
pred_time3 <- end_time3 - start_time3

Cind[4]<- rcorr.cens(test3$pred3, 
                     Surv(test3$timevent, test3$event))["C Index"]

CindL[4]<-rcorr.cens(test3$pred3, Surv(test3$timevent, 
                                       test3$event))["C Index"]-
  1.96*rcorr.cens(test3$pred3, Surv(test3$timevent, 
                                    test3$event))["S.D."]
CindU[4]<-rcorr.cens(test3$pred3, Surv(test3$timevent, 
                                       test3$event))["C Index"]+
  1.96*rcorr.cens(test3$pred3, Surv(test3$timevent, 
                                    test3$event))["S.D."]

# Calibration (Brier score)

BS[4]<-sbrier(obj = Surv(test3$timevent, test3$event), 
              pred = as.numeric(test3$pred3), btime = 2.5)

# Calibration plot
set.seed(7678324)
val_ests3 <- val.surv(est.surv = test3$pred3, 
                      S = Surv(test3$timevent, test3$event),
                      u=2.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.5, 1)))
plot(val_ests3,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test3$pred3, S = Surv(test3$timevent,test3$event), 
        g=10,u=2, pl=T, add=T,lty=0, col=cols[4], 
        cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[4], "black", "black"), cex = 0.8)


# # LT = 2

test4<-read.dta("validation_24months.dta")
test4$smoke<-factor(test4$smoke)
test4$firsttreat<-factor(test4$firsttreat)
test4$pgen<-factor(test4$pgen)
test4$renal<-factor(test4$renal)
test4$lung<-factor(test4$lung)
test4$diabetes<-factor(test4$diabetes)
test4$steroids<-factor(test4$steroids)

start_time3 <- Sys.time()
test4$pred4<-predictSurvProb(object = si.obj2, newdata=test4, 
                             times = 3)

end_time3 <- Sys.time()
pred_time3 <- end_time3 - start_time3

Cind[5]<- rcorr.cens(test4$pred4, 
                     Surv(test4$timevent, test4$event))["C Index"]

CindL[5]<-rcorr.cens(test4$pred4, Surv(test4$timevent, 
                                       test4$event))["C Index"]-
  1.96*rcorr.cens(test4$pred4, Surv(test4$timevent, 
                                    test4$event))["S.D."]
CindU[5]<-rcorr.cens(test4$pred4, Surv(test4$timevent, 
                                       test4$event))["C Index"]+
  1.96*rcorr.cens(test4$pred4, Surv(test4$timevent, 
                                    test4$event))["S.D."]

# Calibration (Brier score)

BS[5]<-sbrier(obj = Surv(test4$timevent, test4$event), 
              pred = as.numeric(test4$pred4), btime = 3)

# Calibration plot
set.seed(76783244)
val_ests4 <- val.surv(est.surv = test4$pred4, 
                      S = Surv(test4$timevent, test4$event),
                      u=3, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.4, 1)))
plot(val_ests4,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test4$pred4, S = Surv(test4$timevent,test4$event), 
        g=10,u=3, pl=T, add=T,lty=0, col=cols[4], 
        cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[4], "black", "black"), cex = 0.8)


# # LT = 2.5

test5<-read.dta("validation_30months.dta")
test5$smoke<-factor(test5$smoke)
test5$firsttreat<-factor(test5$firsttreat)
test5$pgen<-factor(test5$pgen)
test5$renal<-factor(test5$renal)
test5$lung<-factor(test5$lung)
test5$diabetes<-factor(test5$diabetes)
test5$steroids<-factor(test5$steroids)

start_time5 <- Sys.time()
test5$pred5<-predictSurvProb(object = si.obj2, newdata=test5, 
                             times = 3.5)
end_time5 <- Sys.time()
pred_time5 <- end_time5 - start_time5


Cind[6]<- rcorr.cens(test5$pred5, 
                     Surv(test5$timevent, test5$event))["C Index"]

CindL[6]<-rcorr.cens(test5$pred5, Surv(test5$timevent, 
                                       test5$event))["C Index"]-
  1.96*rcorr.cens(test5$pred5, Surv(test5$timevent, 
                                    test5$event))["S.D."]
CindU[6]<-rcorr.cens(test5$pred5, Surv(test5$timevent, 
                                       test5$event))["C Index"]+
  1.96*rcorr.cens(test5$pred5, Surv(test5$timevent, 
                                    test5$event))["S.D."]

# Calibration (Brier score)

BS[6]<-sbrier(obj = Surv(test5$timevent, test5$event), 
              pred = as.numeric(test5$pred5), btime = 3.5)

# Calibration plot
set.seed(767832445)
val_ests5 <- val.surv(est.surv = test5$pred5, 
                      S = Surv(test5$timevent, test5$event),
                      u=3.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.4, 1)))
plot(val_ests5,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test5$pred5, S = Surv(test5$timevent,test5$event), 
        g=10,u=3.5, pl=T, add=T,lty=0, col=cols[4], 
        cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[4], "black", "black"), cex = 0.8)


# # LT = 3

test6<-read.dta("validation_36months.dta")
test6$smoke<-factor(test6$smoke)
test6$firsttreat<-factor(test6$firsttreat)
test6$pgen<-factor(test6$pgen)
test6$renal<-factor(test6$renal)
test6$lung<-factor(test6$lung)
test6$diabetes<-factor(test6$diabetes)
test6$steroids<-factor(test6$steroids)
start_time6 <- Sys.time()
test6$pred6<-predictSurvProb(object = si.obj2, newdata=test6, 
                             times = 4)
end_time6 <- Sys.time()
pred_time6 <- end_time6 - start_time6

Cind[7]<- rcorr.cens(test6$pred6, 
                     Surv(test6$timevent, test6$event))["C Index"]

CindL[7]<-rcorr.cens(test6$pred6, Surv(test6$timevent, 
                                       test6$event))["C Index"]-
  1.96*rcorr.cens(test6$pred6, Surv(test6$timevent, 
                                    test6$event))["S.D."]
CindU[7]<-rcorr.cens(test6$pred6, Surv(test6$timevent, 
                                       test6$event))["C Index"]+
  1.96*rcorr.cens(test6$pred6, Surv(test6$timevent, 
                                    test6$event))["S.D."]

# Calibration (Brier score)

BS[7]<-sbrier(obj = Surv(test6$timevent, test6$event), pred = as.numeric(test6$pred6), btime = 4)

# Calibration plot
set.seed(767832445)

val_ests6 <- val.surv(est.surv = test6$pred6, 
                      S = Surv(test6$timevent, test6$event),
                      u=4, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.4, 1)))
plot(val_ests6,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test6$pred6, S = Surv(test6$timevent,test6$event), 
        g=10,u=4, pl=T, add=T,lty=0, col=cols[4], 
        cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[4], "black", "black"), cex = 0.8)


pred_times<-c(pred_time0, pred_time1, pred_time2, pred_time3, pred_time4, pred_time5, pred_time6)
fit_times<-rep(fit_time, times = 7)
SRF.TA<-cbind(LT, Cind, CindL, CindU, BS, fit_times, pred_times)
write.table(SRF.TA, file = "srfTA.csv", sep = ",", col.names = NA,
            qmethod = "double")



# # Bootstrap function

bootvadSRF<-function(data, B, seed){
  dat <- data
  indiv <- unique(dat$studyno)
  
  c<-0
  c.od<-0
  C<-matrix(0, nrow = B, ncol = 7)
  bs<-matrix(0,nrow = B, ncol = 7)
  C.od<-matrix(0, nrow = B, ncol = 7)
  bs.od<-matrix(0, nrow = B, ncol = 7)
  opt.d<-matrix(0, nrow = B, ncol = 7)
  opt.c<-matrix(0, nrow = B, ncol = 7)
  results<-matrix(0, nrow = 2, ncol = 7)
  
  # original data - temporal assessment
  
  od0 <- test0[,-19]
  od1 <- test1
  od2 <- test2
  od3 <- test3
  od4 <- test4
  od5 <- test5
  od6 <- test6
  
  for (j in 1:B){
    seed<-seed+j
    set.seed(seed)
    smp <- sort(sample(indiv, length(indiv), replace=TRUE))
    smp.df <- data.frame(studyno=smp)
    boot <- merge(smp.df, dat, by = "studyno", all.x=TRUE)
    
    # fitting SRFs
    
    rsf <- rfsrc(Surv(timevent, event) ~ age+pgen+dascore+ovmean
                     +firsttreat+bmi+smoke+steroids+lung+renal+
                       diabetes+disdur+previous_dmards, 
                     data = boot, mtry = rsf.tune$optimal[1], 
                 nodesize = rsf.tune$optimal[2], tree.err=TRUE, 
                 ntree = 300, importance = TRUE)

    # Temporal assessment sets
    
    val0<-boot
    val1<-na.omit(merge(smp.df, od1, by = "studyno", all.x=TRUE))
    val2<-na.omit(merge(smp.df, od2, by = "studyno", all.x=TRUE))
    val3<-na.omit(merge(smp.df, od3, by = "studyno", all.x=TRUE))
    val4<-na.omit(merge(smp.df, od4, by = "studyno", all.x=TRUE))
    val5<-na.omit(merge(smp.df, od5, by = "studyno", all.x=TRUE))
    val6<-na.omit(merge(smp.df, od6, by = "studyno", all.x=TRUE))

    # Survival predictions at different time points (FPM 1)
    
    val0$pred0<-predictSurvProb(object = rsf, newdata=val0, 
                                 times = 1)
    val1$pred1<-predictSurvProb(object = rsf, newdata=val1, 
                                 times = 1.5)
    val2$pred2<-predictSurvProb(object = rsf, newdata=val2, 
                                 times = 2)
    val3$pred3<-predictSurvProb(object = rsf, newdata=val3, 
                                 times = 2.5)
    val4$pred4<-predictSurvProb(object = rsf, newdata=val4, 
                                 times = 3)
    val5$pred5<-predictSurvProb(object = rsf, newdata=val5, 
                                 times = 3.5)
    val6$pred6<-predictSurvProb(object = rsf, newdata=val6, 
                                 times = 4)

    # Discrimination
    
    c[1]<- rcorr.cens(val0$pred0, 
                         Surv(val0$timevent, val0$event))["C Index"]
    c[2]<- rcorr.cens(val1$pred1, 
                      Surv(val1$timevent, val1$event))["C Index"]
    c[3]<- rcorr.cens(val2$pred2, 
                      Surv(val2$timevent, val2$event))["C Index"]
    c[4]<- rcorr.cens(val3$pred3, 
                      Surv(val3$timevent, val3$event))["C Index"]
    c[5]<- rcorr.cens(val4$pred4, 
                      Surv(val4$timevent, val4$event))["C Index"]
    c[6]<- rcorr.cens(val5$pred5, 
                      Surv(val5$timevent, val5$event))["C Index"]
    c[7]<- rcorr.cens(val6$pred6, 
                      Surv(val6$timevent, val6$event))["C Index"]

    C[j,]<-c

    # Calibration measure
    
    bs0<-sbrier(obj = Surv(val0$timevent,val0$event), 
                pred = as.numeric(val0$pred0), btime = 1)
    bs1<-sbrier(obj = Surv(val1$timevent,val1$event), 
                pred = as.numeric(val1$pred1), btime = 1.5)
    bs2<-sbrier(obj = Surv(val2$timevent,val2$event), 
                pred = as.numeric(val2$pred2), btime = 2)
    bs3<-sbrier(obj = Surv(val3$timevent,val3$event), 
                pred = as.numeric(val3$pred3), btime = 2.5)
    bs4<-sbrier(obj = Surv(val4$timevent,val4$event), 
                pred = as.numeric(val4$pred4), btime = 3)
    bs5<-sbrier(obj = Surv(val5$timevent,val5$event), 
                pred = as.numeric(val5$pred5), btime = 3.5)
    bs6<-sbrier(obj = Surv(val6$timevent,val6$event), 
                pred = as.numeric(val6$pred6), btime = 4)
    
    bs[j,]<-c(bs0, bs1, bs2, bs3, bs4, bs5, bs6)

    # Survival predictions - original data
    
    od0$pred0<-predictSurvProb(object = rsf, newdata=od0, 
                                times = 1)
    od1$pred1<-predictSurvProb(object = rsf, newdata=od1, 
                               times = 1.5)
    od2$pred2<-predictSurvProb(object = rsf, newdata=od2, 
                               times = 2)
    od3$pred3<-predictSurvProb(object = rsf, newdata=od3, 
                               times = 2.5)
    od4$pred4<-predictSurvProb(object = rsf, newdata=od4, 
                               times = 3)
    od5$pred5<-predictSurvProb(object = rsf, newdata=od5, 
                               times = 3.5)
    od6$pred6<-predictSurvProb(object = rsf, newdata=od6, 
                               times = 4)

    # Discrimination
    
    c.od[1]<- rcorr.cens(od0$pred0, 
                      Surv(od0$timevent, od0$event))["C Index"]
    c.od[2]<- rcorr.cens(od1$pred1, 
                      Surv(od1$timevent, od1$event))["C Index"]
    c.od[3]<- rcorr.cens(od2$pred2, 
                      Surv(od2$timevent, od2$event))["C Index"]
    c.od[4]<- rcorr.cens(od3$pred3, 
                      Surv(od3$timevent, od3$event))["C Index"]
    c.od[5]<- rcorr.cens(od4$pred4, 
                      Surv(od4$timevent, od4$event))["C Index"]
    c.od[6]<- rcorr.cens(od5$pred5, 
                      Surv(od5$timevent, od5$event))["C Index"]
    c.od[7]<- rcorr.cens(od6$pred6, 
                      Surv(od6$timevent, od6$event))["C Index"]
    
    C.od[j,]<-c.od

    # Calibration measure
    
    bs0.od<-sbrier(obj = Surv(od0$timevent,od0$event), 
                pred = as.numeric(od0$pred0), btime = 1)
    bs1.od<-sbrier(obj = Surv(od1$timevent,od1$event), 
                   pred = as.numeric(od1$pred1), btime = 1.5)
    bs2.od<-sbrier(obj = Surv(od2$timevent,od2$event), 
                   pred = as.numeric(od2$pred2), btime = 2)
    bs3.od<-sbrier(obj = Surv(od3$timevent,od3$event), 
                   pred = as.numeric(od3$pred3), btime = 2.5)
    bs4.od<-sbrier(obj = Surv(od4$timevent,od4$event), 
                   pred = as.numeric(od4$pred4), btime = 3)
    bs5.od<-sbrier(obj = Surv(od5$timevent,od5$event), 
                   pred = as.numeric(od5$pred5), btime = 3.5)
    bs6.od<-sbrier(obj = Surv(od6$timevent,od6$event), 
                   pred = as.numeric(od6$pred6), btime = 4)
    
    bs.od[j,]<-c(bs0.od, bs1.od, bs2.od, bs3.od, bs4.od, bs5.od, bs6.od)

    opt.d[j,]<-as.numeric(C[j,]-C.od[j,])
    opt.c[j,]<-as.numeric(bs.od[j,]-bs[j,])
  }
  
  results[1,]<-c(mean(opt.d[,1]), mean(opt.d[,2]), mean(opt.d[,3]), mean(opt.d[,4]), 
                  mean(opt.d[,5]), mean(opt.d[,6]), mean(opt.d[,7]))

  results[2,]<-c(mean(opt.c[,1]), mean(opt.c[,2]), mean(opt.c[,3]), mean(opt.c[,4]), 
                  mean(opt.c[,5]), mean(opt.c[,6]), mean(opt.c[,7]))

  results<-data.frame(LT0 = results[,1], LT0.5 = results[,2], LT1 = results[,3], 
                       LT1.5 = results[,4], LT2 = results[,5], LT2.5 = results[,6], 
                       LT3 = results[,7], row.names = c("C statistic (OE)", 
                                                         "Brier Score (OE)"))

  print(results)
}

####################################################################

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-2) #always useful to keep one core free
start_time<-Sys.time()
registerDoParallel(cl)
SRF.opt <- foreach(i=1:200, .combine=rbind, .packages = c("ipred", "MASS", "randomForestSRC", 
                                                          "pec", "Hmisc", "dplyr")) %dopar% {
                                                            seed<-199429+i
                                                            bootvadSRF(SRF, 1, seed = seed)
                                                          }
stopCluster(cl)
end_time<-Sys.time()
bootcomptimeSRF <- end_time - start_time

SRF.opt.C=SRF.opt[seq(1,199,2),]
SRF.opt.B=SRF.opt[seq(2,200,2),]
SRF.opt2<-matrix(0, nrow=2, ncol=7)
SRF.opt2[1,]=c(mean(SRF.opt.C[,1]), mean(SRF.opt.C[,2]),
                mean(SRF.opt.C[,3]), mean(SRF.opt.C[,4]),
                mean(SRF.opt.C[,5]), mean(SRF.opt.C[,6]),
                mean(SRF.opt.C[,7]))
SRF.opt2[2,]=c(mean(SRF.opt.B[,1]), mean(SRF.opt.B[,2]),
                mean(SRF.opt.B[,3]), mean(SRF.opt.B[,4]),
                mean(SRF.opt.B[,5]), mean(SRF.opt.B[,6]),
                mean(SRF.opt.B[,7]))

SRF.OPT<-data.frame(LT0 = SRF.opt2[,1], LT0.5 = SRF.opt2[,2], LT1 = SRF.opt2[,3], 
                     LT1.5 = SRF.opt2[,4], LT2 = SRF.opt2[,5], LT2.5 = SRF.opt2[,6], 
                     LT3 = SRF.opt2[,7], row.names = c("C statistic (OE)", 
                                                        "Brier Score (OE)"))

write.table(SRF.OPT, file = "srfOPT.csv", sep = ",", col.names = NA,
            qmethod = "double")

