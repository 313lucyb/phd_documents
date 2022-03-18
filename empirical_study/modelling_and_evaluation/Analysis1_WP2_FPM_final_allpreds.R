library(foreign)
library(rstpm2)
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
#library(plyr)
library(boot)
library(risksetROC)
library(colorRamps)
library(grDevices)
library(ipred)
library(dplyr)
library(RColorBrewer)
cols <-brewer.pal(11, "RdBu")

# CPM development (With only event data in first 12 months)

FPM <- read.dta("baseline_model.dta")
FPM$smoke<-factor(FPM$smoke)
FPM$firsttreat<-factor(FPM$firsttreat)

# Final FPM

set.seed(76735)
start_time_fit1 <- Sys.time()
m2<-stpm2(Surv(timevent, event)~ age + pgen + firsttreat + 
            disdur + previous_dmards + lung + 
            diabetes + bmi + renal + steroids + ovmean + 
            smoke + dascore, data = FPM)
end_time_fit1 <- Sys.time()
fit_time1<- end_time - start_time

# Compare coefficients with Cox model

m2.cox<-coxph(Surv(timevent, event)~age + pgen + firsttreat + 
                disdur + previous_dmards + lung + 
                diabetes + bmi + renal + steroids + ovmean + 
                smoke + dascore, data = FPM)
m2.cox$coef
m2.coef<-as.numeric(summary(m2)@coef[2:18,1])
cor(m2.coef,m2.cox$coef)

# Plot baseline hazard and survival functions

base<- data.frame(age=0, pgen=0, firsttreat="1", 
                    disdur=as.integer(0), previous_dmards =as.integer(0), lung=0,
                    diabetes=0, bmi=0, renal=0, steroids=as.integer(0), ovmean=0, 
                    smoke="0", dascore=0)

par(pty="s")

plot(m2,type="hazard",newdata=base, xlab="Time since baseline (years)")
plot(m2,type="surv",newdata=base, xlab = "Time since baseline (years)")

# CPM development (With event data - 4 years)

FPM2 <- read.dta("baseline_model2.dta")
FPM2$smoke<-factor(FPM2$smoke)
FPM2$firsttreat<-factor(FPM2$firsttreat)

# Final FPM
start_time_fit2 <- Sys.time()
m3<-stpm2(Surv(timevent, event)~age + pgen + firsttreat +
            disdur + previous_dmards + lung + 
            diabetes + bmi + renal + steroids + ovmean + 
            smoke + dascore, data = FPM2)
end_time_fit2 <- Sys.time()
fit_time2 <- end_time_fit2 - start_time_fit2

# Compare coefficients with Cox model
m3.cox<-coxph(Surv(timevent, event)~age + pgen + firsttreat + 
                disdur + previous_dmards + lung + 
                diabetes + bmi + renal + steroids + ovmean + 
                smoke + dascore, data = FPM2)
m3.cox$coef
m3.coef<-as.numeric(summary(m3)@coef[2:18,1])
cor(m3.coef,m3.cox$coef)

#Plot baseline hazard and survival functions

plot(m3,type="hazard",newdata=base, xlab = "Time since baseline (years)")
plot(m3,type="surv",newdata=base, xlab = "Time since baseline (years)")

###########################################################
##              Temporal Assessment                     ##
###########################################################

# CPM with 12 month imposed censoring

## PA measure vectors

LT<-c(0,0.5,1,1.5,2,2.5, 3)
Cslope<-0
CslopeL<-0
CslopeU<-0
Cstat<-0
CstatL<-0
CstatU<-0
BS<-0

## LT = 0 ##

test0<- FPM
test0$timevent<-rep(1,nrow(test0))
testFPM<-test0

t0sub<-test0 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
  disdur, previous_dmards, lung, diabetes, bmi, renal, steroids, ovmean,
  smoke1, smoke2, dascore)

lpmat0<-as.matrix(t0sub)
lp0<-m2.coef %*% t(lpmat0)
lp0<-as.vector(t(lp0))

FPM$lp0<-lp0
start_time_pred10 <- Sys.time()
pred0<-predict(m2, newdata = test0, type = "surv", 
               se.fit =TRUE)
end_time_pred10 <- Sys.time()
pred_time10 <- end_time_pred10 - start_time_pred10

FPM$survpred<- as.numeric(pred0$Estimate)
FPM$pred<- 1-as.numeric(pred0$Estimate)


#Calibration slope and 95% CI

cslope0<-coxph(Surv(timevent, event)~lp0, data = FPM)
Cslope[1]<-as.numeric(cslope0$coefficients)
CslopeL[1]<-confint(cslope0)[1]
CslopeU[1]<-confint(cslope0)[2]

#Harrell's C-statistic and 95% CI

Cstat[1]<-as.numeric(summary(cslope0)$concordance[1])
CstatL[1]<-as.numeric(summary(cslope0)$concordance[1]
                      -(1.96*summary(cslope0)$concordance[2]))
CstatU[1]<-as.numeric(summary(cslope0)$concordance[1]
                      +(1.96*summary(cslope0)$concordance[2]))

# Brier score 

BS[1]<-sbrier(obj = Surv(FPM$timevent,FPM$event), 
              FPM$survpred, btime = 1)

# Visual assessment of calibration by risk groups

# Calibration plot


val_ests0 <- val.surv(est.surv = FPM$survpred, 
                      S = Surv(FPM$timevent, FPM$event),
                     u=1, fun=function(p)log(-log(p)), 
                     pred = sort(runif(100, 0.6, 1)))
plot(val_ests0,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE, title = "LT = 0", cex.lab=1.2) 
groupkm(FPM$survpred, S = Surv(FPM$timevent, FPM$event), 
        g=10,u=1, pl=T, add=T,lty=0, col=cols[2], 
        cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[2], "black", "black"))


## LT = 0.5 ##

test1<-read.dta("validation_6months.dta")
test1$smoke<-factor(test1$smoke, levels=c("0","1","2"), 
                    ordered=F)
test1$firsttreat<-factor(test1$firsttreat, ordered = F)
testFPM1<-test1
testFPM1$timevent<-rep(1.5,nrow(testFPM1))

t1sub<-test1 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat1<-as.matrix(t1sub)
lp1<-m2.coef %*% t(lpmat1)
lp1<-as.vector(t(lp1))

test1$lp1<-lp1
start_time_pred11 <- Sys.time()
pred1<-predict(m2, newdata = testFPM1, type = "surv", 
               se.fit =TRUE)
end_time_pred11 <- Sys.time()
pred_time11 <- end_time_pred11 - start_time_pred11
test1$survpred<- as.numeric(pred1$Estimate)
test1$pred<- 1-as.numeric(pred1$Estimate)

#Calibration slope and 95% CI

cslope1<-coxph(Surv(timevent, event)~lp1, data = test1)
Cslope[2]<-as.numeric(cslope1$coefficients)
CslopeL[2]<-confint(cslope1)[1]
CslopeU[2]<-confint(cslope1)[2]

#Harrell's C-statistic and 95% CI

Cstat[2]<-as.numeric(summary(cslope1)$concordance[1])
CstatL[2]<-as.numeric(summary(cslope1)$concordance[1]
                      -(1.96*summary(cslope1)$concordance[2]))
CstatU[2]<-as.numeric(summary(cslope1)$concordance[1]
                      +(1.96*summary(cslope1)$concordance[2]))

# Brier score 

BS[2]<-sbrier(obj = Surv(test1$timevent,test1$event), 
              test1$survpred, btime = 1.5)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests1 <- val.surv(est.surv = test1$survpred, 
                      S = Surv(test1$timevent, test1$event),
                      u=1.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests1,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE, main = "LT = 0.5") 
groupkm(test1$survpred, S = Surv(test1$timevent, test1$event), 
        g=10,u=1.5, pl=T, add=T,lty=0, col=cols[2], 
        cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[2], "black", "black"))


## LT = 1 ##

test2<-read.dta("validation_12months.dta")
test2$smoke<-factor(test2$smoke, levels=c("0","1","2"), 
                    ordered=F)
test2$firsttreat<-factor(test2$firsttreat, ordered = F)
testFPM2<-test2
testFPM2$timevent<-rep(2,nrow(testFPM2))

t2sub<-test2 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat2<-as.matrix(t2sub)
lp2<-m2.coef %*% t(lpmat2)
lp2<-as.vector(t(lp2))

test2$lp2<-lp2
start_time_pred12 <- Sys.time()
pred2<-predict(m2, newdata = testFPM2, type = "surv", 
               se.fit =TRUE)
end_time_pred12 <- Sys.time()
pred_time12 <- end_time_pred12 - start_time_pred12
test2$survpred<- as.numeric(pred2$Estimate)
test2$pred<- 1-as.numeric(pred2$Estimate)

#Calibration slope and 95% CI

cslope2<-coxph(Surv(timevent, event)~lp2, data = test2)
Cslope[3]<-as.numeric(cslope2$coefficients)
CslopeL[3]<-confint(cslope2)[1]
CslopeU[3]<-confint(cslope2)[2]

#Harrell's C-statistic and 95% CI

Cstat[3]<-as.numeric(summary(cslope2)$concordance[1])
CstatL[3]<-as.numeric(summary(cslope2)$concordance[1]
                      -(1.96*summary(cslope2)$concordance[2]))
CstatU[3]<-as.numeric(summary(cslope2)$concordance[1]
                      +(1.96*summary(cslope2)$concordance[2]))

# Brier score 

BS[3]<-sbrier(obj = Surv(test2$timevent,test2$event), 
              test2$survpred, btime = 2)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests2 <- val.surv(est.surv = test2$survpred, 
                      S = Surv(test2$timevent, test2$event),
                      u=2, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests2,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE, main = "LT = 1") 
groupkm(test2$survpred, S = Surv(test2$timevent, test2$event), 
        g=10,u=2, pl=T, add=T,lty=0, col=cols[2], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[2], "black", "black"))


## LT = 1.5 ##

test3<-read.dta("validation_18months.dta")
test3$smoke<-factor(test3$smoke, levels=c("0","1","2"), ordered=F)
test3$firsttreat<-factor(test3$firsttreat, ordered = F)
testFPM3<-test3
testFPM3$timevent<-rep(2.5,nrow(testFPM3))

t3sub<-test3 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat3<-as.matrix(t3sub)
lp3<-m2.coef %*% t(lpmat3)
lp3<-as.vector(t(lp3))

test3$lp3<-lp3
start_time_pred13 <- Sys.time()
pred3<-predict(m2, newdata = testFPM3, type = "surv", se.fit =TRUE)
end_time_pred13 <- Sys.time()
pred_time13 <- end_time_pred13 - start_time_pred13
test3$survpred<- as.numeric(pred3$Estimate)
test3$pred<- 1-as.numeric(pred3$Estimate)

#Calibration slope and 95% CI

cslope3<-coxph(Surv(timevent, event)~lp3, data = test3)
Cslope[4]<-as.numeric(cslope3$coefficients)
CslopeL[4]<-confint(cslope3)[1]
CslopeU[4]<-confint(cslope3)[2]

#Harrell's C-statistic and 95% CI

Cstat[4]<-as.numeric(summary(cslope3)$concordance[1])
CstatL[4]<-as.numeric(summary(cslope3)$concordance[1]
                      -(1.96*summary(cslope3)$concordance[2]))
CstatU[4]<-as.numeric(summary(cslope3)$concordance[1]
                      +(1.96*summary(cslope3)$concordance[2]))

# Brier score 

BS[4]<-sbrier(obj = Surv(test3$timevent,test3$event), 
              test3$survpred, btime = 2.5)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests3 <- val.surv(est.surv = test3$survpred, 
                      S = Surv(test3$timevent, test3$event),
                      u=2.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests3,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE, main = "LT = 1.5") 
groupkm(test3$survpred, S = Surv(test3$timevent, test3$event), 
        g=10,u=2.5, pl=T, add=T,lty=0, col=cols[2], cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),lty=c(0,2,1),
       pch=c(19,NA,NA), col=c(cols[2], "black", "black"), bty="n")


## LT = 2 ##

test4<-read.dta("validation_24months.dta")
test4$smoke<-factor(test4$smoke, levels=c("0","1","2"), ordered=F)
test4$firsttreat<-factor(test4$firsttreat, ordered = F)
testFPM4<-test4
testFPM4$timevent<-rep(3,nrow(testFPM4))

t4sub<-test4 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat4<-as.matrix(t4sub)
lp4<-m2.coef %*% t(lpmat4)
lp4<-as.vector(t(lp4))

test4$lp4<-lp4
start_time_pred14 <- Sys.time()
pred4<-predict(m2, newdata = testFPM4, type = "surv", se.fit =TRUE)
end_time_pred14 <- Sys.time()
pred_time14 <- end_time_pred14 - start_time_pred14
test4$survpred<- as.numeric(pred4$Estimate)
test4$pred<- 1-as.numeric(pred4$Estimate)

#Calibration slope and 95% CI

cslope4<-coxph(Surv(timevent, event)~lp4, data = test4)
Cslope[5]<-as.numeric(cslope4$coefficients)
CslopeL[5]<-confint(cslope4)[1]
CslopeU[5]<-confint(cslope4)[2]

#Harrell's C-statistic and 95% CI

Cstat[5]<-as.numeric(summary(cslope4)$concordance[1])
CstatL[5]<-as.numeric(summary(cslope4)$concordance[1]
                      -(1.96*summary(cslope4)$concordance[2]))
CstatU[5]<-as.numeric(summary(cslope4)$concordance[1]
                      +(1.96*summary(cslope4)$concordance[2]))

# Brier score 

BS[5]<-sbrier(obj = Surv(test4$timevent,test4$event), 
              test4$survpred, btime = 3)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests4 <- val.surv(est.surv = test4$survpred, 
                      S = Surv(test4$timevent, test4$event),
                      u=3, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests4,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE, main = "LT = 2") 
groupkm(test4$survpred, S = Surv(test4$timevent, test4$event), 
        g=10,u=3, pl=T, add=T,lty=0, col=cols[2], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA), col=c(cols[2], "black", "black"), 
       bty="n")

## LT = 2.5 ##

test5<-read.dta("validation_30months.dta")
test5$smoke<-factor(test5$smoke, levels=c("0","1","2"), ordered=F)
test5$firsttreat<-factor(test5$firsttreat, ordered = F)
testFPM5<-test5
testFPM5$timevent<-rep(3.5,nrow(testFPM5))

t5sub<-test5 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat5<-as.matrix(t5sub)
lp5<-m2.coef %*% t(lpmat5)
lp5<-as.vector(t(lp5))

test5$lp5<-lp5
start_time_pred15 <- Sys.time()
pred5<-predict(m2, newdata = testFPM5, type = "surv", se.fit =TRUE)
end_time_pred15 <- Sys.time()
pred_time15 <- end_time_pred15 - start_time_pred15
test5$survpred<- as.numeric(pred5$Estimate)
test5$pred<- 1-as.numeric(pred5$Estimate)

#Calibration slope and 95% CI

cslope5<-coxph(Surv(timevent, event)~lp5, data = test5)
Cslope[6]<-as.numeric(cslope5$coefficients)
CslopeL[6]<-confint(cslope5)[1]
CslopeU[6]<-confint(cslope5)[2]

#Harrell's C-statistic and 95% CI

Cstat[6]<-as.numeric(summary(cslope5)$concordance[1])
CstatL[6]<-as.numeric(summary(cslope5)$concordance[1]
                      -(1.96*summary(cslope5)$concordance[2]))
CstatU[6]<-as.numeric(summary(cslope5)$concordance[1]
                      +(1.96*summary(cslope5)$concordance[2]))

# Brier score 

BS[6]<-sbrier(obj = Surv(test5$timevent,test5$event), 
              test5$survpred, btime = 3.5)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests5 <- val.surv(est.surv = test5$survpred, 
                      S = Surv(test5$timevent, test5$event),
                      u=3.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests5,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE, main = "LT = 2.5") 
groupkm(test5$survpred, S = Surv(test5$timevent, test5$event), 
        g=10,u=3.5, pl=T, add=T,lty=0, col=cols[2], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA), col=c(cols[2], "black", "black"), 
       bty="n")

## LT = 3 ##

test6<-read.dta("validation_36months.dta")
test6$smoke<-factor(test6$smoke, levels=c("0","1","2"), ordered=F)
test6$firsttreat<-factor(test6$firsttreat, ordered = F)
testFPM6<-test6
testFPM6$timevent<-rep(4,nrow(testFPM6))

t6sub<-test6 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat6<-as.matrix(t6sub)
lp6<-m2.coef %*% t(lpmat6)
lp6<-as.vector(t(lp6))

test6$lp6<-lp6
start_time_pred16 <- Sys.time()
pred6<-predict(m2, newdata = testFPM6, type = "surv", se.fit =TRUE)
end_time_pred16 <- Sys.time()
pred_time16 <- end_time_pred16 - start_time_pred16
test6$survpred<- as.numeric(pred6$Estimate)
test6$pred<- 1-as.numeric(pred6$Estimate)

#Calibration slope and 95% CI

cslope6<-coxph(Surv(timevent, event)~lp6, data = test6)
Cslope[7]<-as.numeric(cslope6$coefficients)
CslopeL[7]<-confint(cslope6)[1]
CslopeU[7]<-confint(cslope6)[2]

#Harrell's C-statistic and 95% CI

Cstat[7]<-as.numeric(summary(cslope6)$concordance[1])
CstatL[7]<-as.numeric(summary(cslope6)$concordance[1]
                      -(1.96*summary(cslope6)$concordance[2]))
CstatU[7]<-as.numeric(summary(cslope6)$concordance[1]
                      +(1.96*summary(cslope6)$concordance[2]))

# Brier score 

BS[7]<-sbrier(obj = Surv(test6$timevent,test6$event), 
              test6$survpred, btime = 4)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests6 <- val.surv(est.surv = test6$survpred, 
                      S = Surv(test6$timevent, test6$event),
                      u=3.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests6,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE, main = "LT = 3") 
groupkm(test6$survpred, S = Surv(test6$timevent, test6$event), 
        g=10,u=3.5, pl=T, add=T,lty=0, col=cols[2], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA), col=c(cols[2], "black", "black"), 
       bty="n")

pred_times_one <- c(pred_time10, pred_time11, pred_time12, pred_time13, pred_time14, pred_time15, pred_time16)
fit_times_one <- rep(fit_time1, times=7)
FPM.one<-data.frame(LT,Cslope,CslopeL,CslopeU,Cstat,CstatL, CstatU, BS, fit_times_one, pred_times_one)
FPM.one

# CPM with censoring imposed after 4 years

## PA measure vectors

LT<-c(0,0.5,1,1.5,2, 2.5, 3)
Cslope<-0
CslopeL<-0
CslopeU<-0
Cstat<-0
CstatL<-0
CstatU<-0
BS<-0

## LT = 0 ##

test0<- FPM
test0FPM<-test0
test0FPM$timevent<-rep(1,nrow(test0FPM))

t0sub<-test0 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat0<-as.matrix(t0sub)
lp0<-m3.coef %*% t(lpmat0)
lp0<-as.vector(t(lp0))

FPM$lp0<-lp0
start_time_pred20 <- Sys.time()
pred0<-predict(m3, newdata = test0FPM, type = "surv", 
               se.fit =TRUE)
end_time_pred20 <- Sys.time()
pred_time20 <- end_time_pred20 - start_time_pred20
test0$survpred<- as.numeric(pred0$Estimate)
test0$pred<- 1-as.numeric(pred0$Estimate)

#Calibration slope and 95% CI

cslope0<-coxph(Surv(timevent, event)~lp0, data = FPM)
Cslope[1]<-as.numeric(cslope0$coefficients)
CslopeL[1]<-confint(cslope0)[1]
CslopeU[1]<-confint(cslope0)[2]

#Harrell's C-statistic and 95% CI

Cstat[1]<-as.numeric(summary(cslope0)$concordance[1])
CstatL[1]<-as.numeric(summary(cslope0)$concordance[1]
                      -(1.96*summary(cslope0)$concordance[2]))
CstatU[1]<-as.numeric(summary(cslope0)$concordance[1]
                      +(1.96*summary(cslope0)$concordance[2]))

# Brier score 

BS[1]<-sbrier(obj = Surv(FPM$timevent,FPM$event), 
              pred = test0$survpred, btime = 1)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests0 <- val.surv(est.surv = test0$survpred, 
                      S = Surv(FPM$timevent, FPM$event),
                      u=1, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests0,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test0$survpred, S = Surv(FPM$timevent,FPM$event), 
        g=10,u=1, pl=T, add=T,lty=0, col=cols[3], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n",
       col=c(cols[3], "black", "black"))


## LT = 0.5 ##

test1<-read.dta("validation_6months.dta")
test1$smoke<-factor(test1$smoke)
test1$firsttreat<-factor(test1$firsttreat)
testFPM1<-test1
testFPM1$timevent<-rep(1.5,nrow(testFPM1))

t1sub<-test1 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat1<-as.matrix(t1sub)
lp1<-m3.coef %*% t(lpmat1)
lp1<-as.vector(t(lp1)) 

test1$lp1<-lp1
start_time_pred21 <- Sys.time()
pred1<-predict(m3, newdata = testFPM1, type = "surv", 
               se.fit =TRUE)
end_time_pred21 <- Sys.time()
pred_time21 <- end_time_pred21 - start_time_pred21
test1$survpred<- as.numeric(pred1$Estimate)
test1$pred<- 1-as.numeric(pred1$Estimate)

#Calibration slope and 95% CI

cslope1<-coxph(Surv(timevent, event)~lp1, data = test1)
Cslope[2]<-as.numeric(cslope1$coefficients)
CslopeL[2]<-confint(cslope1)[1]
CslopeU[2]<-confint(cslope1)[2]

#Harrell's C-statistic and 95% CI

Cstat[2]<-as.numeric(summary(cslope1)$concordance[1])
CstatL[2]<-as.numeric(summary(cslope1)$concordance[1]
                      -(1.96*summary(cslope1)$concordance[2]))
CstatU[2]<-as.numeric(summary(cslope1)$concordance[1]
                      +(1.96*summary(cslope1)$concordance[2]))

# Brier score 

BS[2]<-sbrier(obj = Surv(test1$timevent,test1$event), 
              pred = test1$survpred, btime = 1.5)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests1 <- val.surv(est.surv = test1$survpred, 
                      S = Surv(test1$timevent, test1$event),
                      u=1.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests1,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test1$survpred, S = Surv(test1$timevent, test1$event), 
        g=10,u=1.5, pl=T, add=T,lty=0, col=cols[3], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[3], "black", "black"))


## LT = 1 ##

test2<-read.dta("validation_12months.dta")
test2$smoke<-factor(test2$smoke)
test2$firsttreat<-factor(test2$firsttreat)
testFPM2<-test2
testFPM2$timevent<-rep(2,nrow(testFPM2))

t2sub<-test2 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat2<-as.matrix(t2sub)
lp2<-m3.coef %*% t(lpmat2)
lp2<-as.vector(t(lp2))

test2$lp2<-lp2
start_time_pred22 <- Sys.time()
pred2<-predict(m3, newdata = testFPM2, type = "surv", 
               se.fit =TRUE)
end_time_pred22 <- Sys.time()
pred_time22 <- end_time_pred22 - start_time_pred22
test2$survpred<- as.numeric(pred2$Estimate)
test2$pred<- 1-as.numeric(pred2$Estimate)

#Calibration slope and 95% CI

cslope2<-coxph(Surv(timevent, event)~lp2, data = test2)
Cslope[3]<-as.numeric(cslope2$coefficients)
CslopeL[3]<-confint(cslope2)[1]
CslopeU[3]<-confint(cslope2)[2]

#Harrell's C-statistic and 95% CI

Cstat[3]<-as.numeric(summary(cslope2)$concordance[1])
CstatL[3]<-as.numeric(summary(cslope2)$concordance[1]
                      -(1.96*summary(cslope2)$concordance[2]))
CstatU[3]<-as.numeric(summary(cslope2)$concordance[1]
                      +(1.96*summary(cslope2)$concordance[2]))

# Brier score 

BS[3]<-sbrier(obj = Surv(test2$timevent,test2$event), 
              test2$survpred, btime = 2)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests2 <- val.surv(est.surv = test2$survpred, 
                      S = Surv(test2$timevent, test2$event),
                      u=2, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests2,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test2$survpred, S = Surv(test2$timevent, test2$event), 
        g=10,u=2, pl=T, add=T,lty=0, col=cols[3], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA),bty="n", 
       col=c(cols[3], "black", "black"))


## LT = 1.5 ##

test3<-read.dta("validation_18months.dta")
test3$smoke<-factor(test3$smoke)
test3$firsttreat<-factor(test3$firsttreat)
testFPM3<-test3
testFPM3$timevent<-rep(2.5,nrow(testFPM3))

t3sub<-test3 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat3<-as.matrix(t3sub)
lp3<-m3.coef %*% t(lpmat3)
lp3<-as.vector(t(lp3))

test3$lp3<-lp3
start_time_pred23<- Sys.time()
pred3<-predict(m3, newdata = testFPM3, type = "surv", se.fit =TRUE)
end_time_pred23 <- Sys.time()
pred_time23 <- end_time_pred23 - start_time_pred23
test3$survpred<- as.numeric(pred3$Estimate)
test3$pred<- 1-as.numeric(pred3$Estimate)

#Calibration slope and 95% CI

cslope3<-coxph(Surv(timevent, event)~lp3, data = test3)
Cslope[4]<-as.numeric(cslope3$coefficients)
CslopeL[4]<-confint(cslope3)[1]
CslopeU[4]<-confint(cslope3)[2]

#Harrell's C-statistic and 95% CI

Cstat[4]<-as.numeric(summary(cslope3)$concordance[1])
CstatL[4]<-as.numeric(summary(cslope3)$concordance[1]
                      -(1.96*summary(cslope3)$concordance[2]))
CstatU[4]<-as.numeric(summary(cslope3)$concordance[1]
                      +(1.96*summary(cslope3)$concordance[2]))

# Brier score 

BS[4]<-sbrier(obj = Surv(test3$timevent,test3$event), 
              test3$survpred, btime = 2.5)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests3 <- val.surv(est.surv = test3$survpred, 
                      S = Surv(test3$timevent, test3$event),
                      u=2.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests3,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test3$survpred, S = Surv(test3$timevent, test3$event), 
        g=10,u=2.5, pl=T, add=T,lty=0, col=cols[3], cex.subtitle=FALSE)
legend(0.85,0.7,c("Risk groups","Reference line","95% CI"),lty=c(0,2,1),
       pch=c(19,NA,NA), col=c(cols[3], "black", "black"), bty="n")


## LT = 2 ##

test4<-read.dta("validation_24months.dta")
test4$smoke<-factor(test4$smoke)
test4$firsttreat<-factor(test4$firsttreat)
testFPM4<-test4
testFPM4$timevent<-rep(3,nrow(testFPM4))

t4sub<-test4 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat4<-as.matrix(t4sub)
lp4<-m3.coef %*% t(lpmat4)
lp4<-as.vector(t(lp4))

test4$lp4<-lp4
start_time_pred24 <- Sys.time()
pred4<-predict(m3, newdata = testFPM4, type = "surv", se.fit =TRUE)
end_time_pred24 <- Sys.time()
pred_time24 <- end_time_pred24 - start_time_pred24
test4$survpred<- as.numeric(pred4$Estimate)
test4$pred<- 1-as.numeric(pred4$Estimate)

#Calibration slope and 95% CI

cslope4<-coxph(Surv(timevent, event)~lp4, data = test4)
Cslope[5]<-as.numeric(cslope4$coefficients)
CslopeL[5]<-confint(cslope4)[1]
CslopeU[5]<-confint(cslope4)[2]

#Harrell's C-statistic and 95% CI

Cstat[5]<-as.numeric(summary(cslope4)$concordance[1])
CstatL[5]<-as.numeric(summary(cslope4)$concordance[1]
                      -(1.96*summary(cslope4)$concordance[2]))
CstatU[5]<-as.numeric(summary(cslope4)$concordance[1]
                      +(1.96*summary(cslope4)$concordance[2]))

# Brier score 

BS[5]<-sbrier(obj = Surv(test4$timevent,test4$event), 
              test4$survpred, btime = 3)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests4 <- val.surv(est.surv = test4$survpred, 
                      S = Surv(test4$timevent, test4$event),
                      u=3, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests4,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test4$survpred, S = Surv(test4$timevent, test4$event), 
        g=10,u=3, pl=T, add=T,lty=0, col=cols[3], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA), col=c(cols[3], "black", "black"), 
       bty="n")

## LT = 2.5 ##

test5<-read.dta("validation_30months.dta")
test5$smoke<-factor(test5$smoke)
test5$firsttreat<-factor(test5$firsttreat)
testFPM5<-test5
testFPM5$timevent<-rep(3.5,nrow(testFPM5))

t5sub<-test5 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat5<-as.matrix(t5sub)
lp5<-m3.coef %*% t(lpmat5)
lp5<-as.vector(t(lp5)) 

test5$lp5<-lp5
start_time_pred25 <- Sys.time()
pred5<-predict(m3, newdata = testFPM5, type = "surv", se.fit =TRUE)
end_time_pred25 <- Sys.time()
pred_time25 <- end_time_pred25 - start_time_pred25
test5$survpred<- as.numeric(pred5$Estimate)
test5$pred<- 1-as.numeric(pred5$Estimate)

#Calibration slope and 95% CI

cslope5<-coxph(Surv(timevent, event)~lp5, data = test5)
Cslope[6]<-as.numeric(cslope5$coefficients)
CslopeL[6]<-confint(cslope5)[1]
CslopeU[6]<-confint(cslope5)[2]

#Harrell's C-statistic and 95% CI

Cstat[6]<-as.numeric(summary(cslope5)$concordance[1])
CstatL[6]<-as.numeric(summary(cslope5)$concordance[1]
                      -(1.96*summary(cslope5)$concordance[2]))
CstatU[6]<-as.numeric(summary(cslope5)$concordance[1]
                      +(1.96*summary(cslope5)$concordance[2]))

# Brier score 

BS[6]<-sbrier(obj = Surv(test5$timevent,test5$event), 
              test5$survpred, btime = 3.5)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests5 <- val.surv(est.surv = test5$survpred, 
                      S = Surv(test5$timevent, test5$event),
                      u=3.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests5,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test5$survpred, S = Surv(test5$timevent, test5$event), 
        g=10,u=3.5, pl=T, add=T,lty=0, col=cols[3], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA), col=c(cols[3], "black", "black"), 
       bty="n")

## LT = 3 ##

test6<-read.dta("validation_36months.dta")
test6$smoke<-factor(test6$smoke)
test6$firsttreat<-factor(test6$firsttreat)
testFPM6<-test6
testFPM6$timevent<-rep(4,nrow(testFPM6))

t6sub<-test6 %>%
  dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                disdur, previous_dmards, lung, 
                diabetes, bmi, renal, steroids, ovmean,
                smoke1, smoke2, dascore)

lpmat6<-as.matrix(t6sub)
lp6<-m3.coef %*% t(lpmat6)
lp6<-as.vector(t(lp6))

test6$lp6<-lp6
start_time_pred26 <- Sys.time()
pred6<-predict(m3, newdata = testFPM6, type = "surv", se.fit =TRUE)
end_time_pred26 <- Sys.time()
pred_time26 <- end_time_pred26 - start_time_pred26
test6$survpred<- as.numeric(pred6$Estimate)
test6$pred<- 1-as.numeric(pred6$Estimate)

#Calibration slope and 95% CI

cslope6<-coxph(Surv(timevent, event)~lp6, data = test6)
Cslope[7]<-as.numeric(cslope6$coefficients)
CslopeL[7]<-confint(cslope6)[1]
CslopeU[7]<-confint(cslope6)[2]

#Harrell's C-statistic and 95% CI

Cstat[7]<-as.numeric(summary(cslope6)$concordance[1])
CstatL[7]<-as.numeric(summary(cslope6)$concordance[1]
                      -(1.96*summary(cslope6)$concordance[2]))
CstatU[7]<-as.numeric(summary(cslope6)$concordance[1]
                      +(1.96*summary(cslope6)$concordance[2]))

# Brier score 

BS[7]<-sbrier(obj = Surv(test6$timevent,test6$event), 
              test6$survpred, btime = 4)

# Visual assessment of calibration by risk groups

# Calibration plot

val_ests6 <- val.surv(est.surv = test6$survpred, 
                      S = Surv(test6$timevent, test6$event),
                      u=3.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.6, 1)))
plot(val_ests6,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = FALSE) 
groupkm(test6$survpred, S = Surv(test6$timevent, test6$event), 
        g=10,u=3.5, pl=T, add=T,lty=0, col=cols[3], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),
       pch=c(19,NA,NA), col=c(cols[3], "black", "black"), 
       bty="n")

pred_times_two <- c(pred_time20, pred_time21, pred_time22, pred_time23, pred_time24, pred_time25, pred_time26)
fit_times_two <- rep(fit_time2, times = 7)
FPM.two<-data.frame(LT,Cslope,CslopeL,CslopeU,Cstat,CstatL, CstatU, BS, fit_times_two, pred_times_two)
FPM.two

write.table(FPM.one, file = "fpmTA_1.csv", sep = ",", col.names = NA,
            qmethod = "double")
write.table(FPM.two, file = "fpmTA_2.csv", sep = ",", col.names = NA,
            qmethod = "double")

## Optimism estimates (only stores FPM 2)

#Bootstrap function (FPM) 

bootvadFPM<-function(data, B, seed){
  dat2 <- data
  dat2$studyno <- dat2$groupid
  
  #sub1 <- dat2 %>% filter(trt1 == 1) %>% select(studyno)
  #sub2 <- dat2 %>% filter(trt2 == 1) %>% select(studyno)
  #sub3 <- dat2 %>% filter(trt4 == 1) %>% select(studyno)
  #sub4 <- dat2 %>% filter(trt1561 == 1) %>% select(studyno)
  #sub5 <- dat2 %>% filter(trt10060 == 1) %>% select(studyno)
  
  c2<-0
  c2od<-0
  C2<-matrix(0, nrow = B, ncol = 7)
  bs2<-matrix(0,nrow = B, ncol = 7)
  C2.od<-matrix(0, nrow = B, ncol = 7)
  bs.od2<-matrix(0, nrow = B, ncol = 7)
  opt.d2<-matrix(0, nrow = B, ncol = 7)
  opt.c2<-matrix(0, nrow = B, ncol = 7)
  results2<-matrix(0, nrow = 2, ncol = 7)
  
  # original data - temporal assessment
  
  od0 <- test0[-19]
  od0pred <- od0
  od0pred$timevent <- rep(1, nrow(od0pred))
  od1 <- test1
  od1pred <- od1
  od1pred$timevent <- rep(1.5, nrow(od1pred))
  od2 <- test2
  od2pred <- od2
  od2pred$timevent <- rep(2, nrow(od2pred))
  od3 <- test3
  od3pred <- od3
  od3pred$timevent <- rep(2.5, nrow(od3pred))
  od4 <- test4
  od4pred <- od4
  od4pred$timevent <- rep(3, nrow(od4pred))
  od5 <- test5
  od5pred <- od5
  od5pred$timevent <- rep(3.5, nrow(od5pred))
  od6 <- test6
  od6pred <- od6
  od6pred$timevent <- rep(4, nrow(od6pred))
  
  od0$studyno<-od0$groupid
  od1$studyno<-od1$groupid
  od2$studyno<-od2$groupid
  od3$studyno<-od3$groupid
  od4$studyno<-od4$groupid
  od5$studyno<-od5$groupid
  od6$studyno<-od6$groupid
  
  for (j in 1:B){
    seed = seed + j
    set.seed(seed)
    smp <- sort(sample(dat2$studyno, length(dat2$studyno), replace=TRUE))
    smp.df <- data.frame(studyno=smp)
    boot2 <- as.data.frame(merge(smp.df, dat2, by = "studyno", all.x=TRUE))

    boot2$firsttreat <- factor(boot2$firsttreat, levels = c("1", "2", "4", "1561", "10060"))
    boot2$smoke <- factor(boot2$smoke, levels = c("0", "1", "2"))

    bm2<-stpm2(Surv(timevent, event) ~ age + pgen + firsttreat + 
                 disdur + previous_dmards + lung + 
                 diabetes + bmi + renal + steroids + ovmean + 
                 smoke + dascore, data = boot2)

    
    # Temporal assessment sets

    val0<-boot2
    val0pred <- val0
    val0pred$timevent <- rep(1, times=nrow(val0pred))
    val1<-na.omit(merge(smp.df, od1, by = "studyno", all.x=TRUE))
    val1pred <- val1
    head(od1)
    val1pred$timevent <- rep(1.5, times=nrow(val1pred))
    val2<-na.omit(merge(smp.df, od2, by = "studyno", all.x=TRUE))
    val2pred <- val2
    val2pred$timevent <- rep(2, times=nrow(val2pred))
    val3<-na.omit(merge(smp.df, od3, by = "studyno", all.x=TRUE))
    val3pred <- val3
    val3pred$timevent <- rep(2.5, times=nrow(val3pred))
    val4<-na.omit(merge(smp.df, od4, by = "studyno", all.x=TRUE))
    val4pred <- val4
    val4pred$timevent <- rep(3, times=nrow(val4pred))
    val5<-na.omit(merge(smp.df, od5, by = "studyno", all.x=TRUE))
    val5pred <- val5
    val5pred$timevent <- rep(3.5, times=nrow(val5pred))
    val6<-na.omit(merge(smp.df, od6, by = "studyno", all.x=TRUE))
    val6pred <- val6
    val6pred$timevent <- rep(4, times=nrow(val6pred))

    # Survival predictions at different time points (FPM 2)
    
    pred02 <-predict(bm2, newdata = val0pred, type = "surv", 
                     se.fit =TRUE)
    val0$survpred2<- as.numeric(pred02$Estimate)
    val0$pred2<- 1-as.numeric(pred02$Estimate)

    pred12 <-predict(bm2, newdata = val1pred, type = "surv", 
                     se.fit =TRUE)
    val1$survpred2<- as.numeric(pred12$Estimate)
    val1$pred2<- 1-as.numeric(pred12$Estimate)
    
    pred22 <-predict(bm2, newdata = val2pred, type = "surv", 
                     se.fit =TRUE)
    val2$survpred2<- as.numeric(pred22$Estimate)
    val2$pred2<- 1-as.numeric(pred22$Estimate)
    
    pred32 <-predict(bm2, newdata = val3pred, type = "surv", 
                     se.fit =TRUE)
    val3$survpred2<- as.numeric(pred32$Estimate)
    val3$pred2<- 1-as.numeric(pred32$Estimate)
    
    pred42 <-predict(bm2, newdata = val4pred, type = "surv", 
                     se.fit =TRUE)
    val4$survpred2<- as.numeric(pred42$Estimate)
    val4$pred2<- 1-as.numeric(pred42$Estimate)
    
    pred52 <-predict(bm2, newdata = val5pred, type = "surv", 
                     se.fit =TRUE)
    val5$survpred2<- as.numeric(pred52$Estimate)
    val5$pred2<- 1-as.numeric(pred52$Estimate)
    
    pred62 <-predict(bm2, newdata = val6pred, type = "surv", 
                     se.fit =TRUE)
    val6$survpred2<- as.numeric(pred62$Estimate)
    val6$pred2<- 1-as.numeric(pred62$Estimate)

    # Linear predictors (FPM 2)
    
    bm2coef<-as.numeric(summary(bm2)@coef[2:18,1])

    val0sub<-val0 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    lpmat0<-as.matrix(val0sub)
    
    lp02<-bm2coef %*% t(lpmat0)
    lp02<-as.vector(t(lp02))
    
    val1sub<-val1 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    lpmat1<-as.matrix(val1sub)
    

    lp12<-bm2coef %*% t(lpmat1)

    lp12<-as.vector(t(lp12))
  
    val2sub<-val2 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    lpmat2<-as.matrix(val2sub)
    

    lp22<-bm2coef %*% t(lpmat2)

    lp22<-as.vector(t(lp22))
    
    val3sub<-val3 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    lpmat3<-as.matrix(val3sub)
    

    lp32<-bm2coef %*% t(lpmat3)

    lp32<-as.vector(t(lp32))
    
    val4sub<-val4 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    lpmat4<-as.matrix(val4sub)
    

    lp42<-bm2coef %*% t(lpmat4)

    lp42<-as.vector(t(lp42))
    
    val5sub<-val5 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    lpmat5<-as.matrix(val5sub)
    

    lp52<-bm2coef %*% t(lpmat5)

    lp52<-as.vector(t(lp52))
    
    val6sub<-val6 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    lpmat6<-as.matrix(val6sub)

    lp62<-bm2coef %*% t(lpmat6)

    lp62<-as.vector(t(lp62))
    
    # Discrimination
  
    cmodel02<-coxph(Surv(timevent, event)~lp02, data = val0)

    c2[1]<-as.numeric(summary(cmodel02)$concordance[1])
    

    cmodel12<-coxph(Surv(timevent, event)~lp12, data = val1)

    c2[2]<-as.numeric(summary(cmodel12)$concordance[1])
    

    cmodel22<-coxph(Surv(timevent, event)~lp22, data = val2)

    c2[3]<-as.numeric(summary(cmodel22)$concordance[1])
    

    cmodel32<-coxph(Surv(timevent, event)~lp32, data = val3)

    c2[4]<-as.numeric(summary(cmodel32)$concordance[1])
    

    cmodel42<-coxph(Surv(timevent, event)~lp42, data = val4)

    c2[5]<-as.numeric(summary(cmodel42)$concordance[1])
    

    cmodel52<-coxph(Surv(timevent, event)~lp52, data = val5)

    c2[6]<-as.numeric(summary(cmodel52)$concordance[1])
    

    cmodel62<-coxph(Surv(timevent, event)~lp62, data = val6)

    c2[7]<-as.numeric(summary(cmodel62)$concordance[1])
    

    C2[j,]<-c2
    
    # Calibration measure (FPM 2)
              
    bs02<-sbrier(obj = Surv(val0$timevent,val0$event), 
                 val0$survpred2, btime = 1)
    bs12<-sbrier(obj = Surv(val1$timevent,val1$event), 
                 val1$survpred2, btime = 1.5)
    bs22<-sbrier(obj = Surv(val2$timevent,val2$event), 
                 val2$survpred2, btime = 2)
    bs32<-sbrier(obj = Surv(val3$timevent,val3$event), 
                 val3$survpred2, btime = 2.5)
    bs42<-sbrier(obj = Surv(val4$timevent,val4$event), 
                 val4$survpred2, btime = 3)
    bs52<-sbrier(obj = Surv(val5$timevent,val5$event), 
                 val5$survpred2, btime = 3.5)
    bs62<-sbrier(obj = Surv(val6$timevent,val6$event), 
                 val6$survpred2, btime = 4)
    
    bs2[j,]<-c(bs02, bs12, bs22, bs32, bs42, bs52, bs62)
    
    # Survival predictions - original data
    
    # Survival predictions at different time points (FPM 2)
    
    pred02od <-predict(bm2, newdata = od0pred, type = "surv", 
                     se.fit =TRUE)
    od0$survpred2<- as.numeric(pred02od$Estimate)
    od0$pred2<- 1-as.numeric(pred02od$Estimate)
    
    pred12od <-predict(bm2, newdata = od1pred, type = "surv", 
                     se.fit =TRUE)
    od1$survpred2<- as.numeric(pred12od$Estimate)
    od1$pred2<- 1-as.numeric(pred12od$Estimate)
    
    pred22od <-predict(bm2, newdata = od2pred, type = "surv", 
                     se.fit =TRUE)
    od2$survpred2<- as.numeric(pred22od$Estimate)
    od2$pred2<- 1-as.numeric(pred22od$Estimate)
    
    pred32od <-predict(bm2, newdata = od3pred, type = "surv", 
                     se.fit =TRUE)
    od3$survpred2<- as.numeric(pred32od$Estimate)
    od3$pred2<- 1-as.numeric(pred32od$Estimate)
    
    pred42od <-predict(bm2, newdata = od4pred, type = "surv", 
                     se.fit =TRUE)
    od4$survpred2<- as.numeric(pred42od$Estimate)
    od4$pred2<- 1-as.numeric(pred42od$Estimate)
    
    pred52od <-predict(bm2, newdata = od5pred, type = "surv", 
                     se.fit =TRUE)
    od5$survpred2<- as.numeric(pred52od$Estimate)
    od5$pred2<- 1-as.numeric(pred52od$Estimate)
    
    pred62od <-predict(bm2, newdata = od6pred, type = "surv", 
                     se.fit =TRUE)
    od6$survpred2<- as.numeric(pred62od$Estimate)
    od6$pred2<- 1-as.numeric(pred62od$Estimate)
    
    # Linear predictors (FPM 1 & 2)
    
    bm2coef<-as.numeric(summary(bm2)@coef[2:18,1])
    
    od0sub<-od0 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    odlpmat0<-as.matrix(od0sub)
    

    lp02od<-bm2coef %*% t(odlpmat0)
  
    lp02od<-as.vector(t(lp02od))
    
    od1sub<-od1 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    odlpmat1<-as.matrix(od1sub)
    
  
    lp12od<-bm2coef %*% t(odlpmat1)

    lp12od<-as.vector(t(lp12od))
    
    od2sub<-od2 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    odlpmat2<-as.matrix(od2sub)
    
  
    lp22od<-bm2coef %*% t(odlpmat2)

    lp22od<-as.vector(t(lp22od))
    
    od3sub<-od3 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    odlpmat3<-as.matrix(od3sub)
    
 
    lp32od<-bm2coef %*% t(odlpmat3)

    lp32od<-as.vector(t(lp32od))
    
    od4sub<-od4 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    odlpmat4<-as.matrix(od4sub)
    

    lp42od<-bm2coef %*% t(odlpmat4)

    lp42od<-as.vector(t(lp42od))
    
    od5sub<-od5 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    odlpmat5<-as.matrix(od5sub)
    

    lp52od<-bm2coef %*% t(odlpmat5)

    lp52od<-as.vector(t(lp52od))
    
    od6sub<-od6 %>%
      dplyr::select(age, pgen, trt2, trt4, trt1561, trt10060,
                    disdur, previous_dmards, lung, 
                    diabetes, bmi, renal, steroids, ovmean,
                    smoke1, smoke2, dascore)
    odlpmat6<-as.matrix(od6sub)
    

    lp62od<-bm2coef %*% t(odlpmat6)

    lp62od<-as.vector(t(lp62od))
    
    # Discrimination
    

    cmodel02od<-coxph(Surv(timevent, event)~lp02od, data = od0)

    c2od[1]<-as.numeric(summary(cmodel02od)$concordance[1])
    

    cmodel12od<-coxph(Surv(timevent, event)~lp12od, data = od1)

    c2od[2]<-as.numeric(summary(cmodel12od)$concordance[1])
    

    cmodel22od<-coxph(Surv(timevent, event)~lp22od, data = od2)

    c2od[3]<-as.numeric(summary(cmodel22od)$concordance[1])
    
  
    cmodel32od<-coxph(Surv(timevent, event)~lp32od, data = od3)

    c2od[4]<-as.numeric(summary(cmodel32od)$concordance[1])
    

    cmodel42od<-coxph(Surv(timevent, event)~lp42od, data = od4)

    c2od[5]<-as.numeric(summary(cmodel42od)$concordance[1])
    

    cmodel52od<-coxph(Surv(timevent, event)~lp52od, data = od5)

    c2od[6]<-as.numeric(summary(cmodel52od)$concordance[1])
    

    cmodel62od<-coxph(Surv(timevent, event)~lp62od, data = od6)

    c2od[7]<-as.numeric(summary(cmodel62od)$concordance[1])
    

    C2.od[j,]<-c2od
    
    # Calibration measure (FPM 2)

    bs02od<-sbrier(obj = Surv(od0$timevent,od0$event), 
                   od0$survpred2, btime = 1)
    bs12od<-sbrier(obj = Surv(od1$timevent,od1$event), 
                   od1$survpred2, btime = 1.5)
    bs22od<-sbrier(obj = Surv(od2$timevent,od2$event), 
                   od2$survpred2, btime = 2)
    bs32od<-sbrier(obj = Surv(od3$timevent,od3$event), 
                   od3$survpred2, btime = 2.5)
    bs42od<-sbrier(obj = Surv(od4$timevent,od4$event), 
                   od4$survpred2, btime = 3)
    bs52od<-sbrier(obj = Surv(od5$timevent,od5$event), 
                   od5$survpred2, btime = 3.5)
    bs62od<-sbrier(obj = Surv(od6$timevent,od6$event), 
                   od6$survpred2, btime = 4)
    
    bs.od2[j,]<-c(bs02od, bs12od, bs22od, bs32od, bs42od, bs52od, bs62od)
    
  
    opt.d2[j,]<-as.numeric(C2[j,]-C2.od[j,])
    opt.c2[j,]<-as.numeric(bs.od2[j,]-bs2[j,])
  }

  results2[1,]<-c(mean(opt.d2[,1]), mean(opt.d2[,2]), mean(opt.d2[,3]), mean(opt.d2[,4]), 
                  mean(opt.d2[,5]), mean(opt.d2[,6]), mean(opt.d2[,7]))
  results2[2,]<-c(mean(opt.c2[,1]), mean(opt.c2[,2]), mean(opt.c2[,3]), mean(opt.c2[,4]), 
                  mean(opt.c2[,5]), mean(opt.c2[,6]), mean(opt.c2[,7]))
  results2<-data.frame(LT0 = results2[,1], LT0.5 = results2[,2], LT1 = results2[,3], 
                       LT1.5 = results2[,4], LT2 = results2[,5], LT2.5 = results2[,6], 
                       LT3 = results2[,7], row.names = c("C statistic (OE)", 
                                                        "Brier Score (OE)"))
  results2<-as.data.frame(results2)
  results2
}

####################################################################

library(foreach)
library(doParallel)
seed<-199431
cores=detectCores()
cl <- makeCluster(cores[1]-1) #always useful to keep one core free
start_time<-Sys.time()
registerDoParallel(cl, outfile="")
FPM.opt <- foreach(i=1:200, .combine=rbind, .packages = c("ipred", "MASS", "rstpm2", 
                                                           "survival", "Hmisc", "dplyr")) %dopar% {
                                                             seed<-199429+i
                                                             bootvadFPM(FPM2, 1, seed=seed)
                                                           }
stopCluster(cl)
end_time<-Sys.time()
end_time - start_time

## Concatenating output from foreach loop 

FPM.opt.C=FPM.opt[seq(1,399,2),]
FPM.opt.B=FPM.opt[seq(2,400,2),]
FPM.opt2<-matrix(0, nrow=2, ncol=7)
FPM.opt2[1,]=c(mean(FPM.opt.C[,1]), mean(FPM.opt.C[,2]),
                mean(FPM.opt.C[,3]), mean(FPM.opt.C[,4]),
                mean(FPM.opt.C[,5]), mean(FPM.opt.C[,6]),
                mean(FPM.opt.C[,7]))
FPM.opt2[2,]=c(mean(FPM.opt.B[,1]), mean(FPM.opt.B[,2]),
                mean(FPM.opt.B[,3]), mean(FPM.opt.B[,4]),
                mean(FPM.opt.B[,5]), mean(FPM.opt.B[,6]),
                mean(FPM.opt.B[,7]))

FPM.OPT<-data.frame(LT0 = FPM.opt2[,1], LT0.5 = FPM.opt2[,2], LT1 = FPM.opt2[,3], 
                     LT1.5 = FPM.opt2[,4], LT2 = FPM.opt2[,5], LT2.5 = FPM.opt2[,6], 
                     LT3 = FPM.opt2[,7], row.names = c("C statistic (OE)", 
                                                        "Brier Score (OE)"))

write.table(FPM.OPT, file = "fpm1OPT.csv", sep = ",", col.names = NA,
            qmethod = "double")
