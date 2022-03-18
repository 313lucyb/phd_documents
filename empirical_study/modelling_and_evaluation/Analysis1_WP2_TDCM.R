library(foreign)
library(rstpm2)
library(mfp)
library(rms)
library(ROCR)
library(survival)
library(Hmisc)
library(MASS)
library(mfp)
library(reshape2)
library(flexsurv)
library(pROC)
library(DescTools)
library(plyr)
library(boot)
library(ROCR)
library(risksetROC)
library(colorRamps)
library(grDevices)
library(ipred)
library(dplyr)
library(RColorBrewer)

cols <-brewer.pal(11, "RdBu")

### Example with exposed subjects
## Standard - no previous infection predictor

TDCM <- read.dta("long_TDCM.dta")
TDCM$smoke<-factor(TDCM$smoke)
TDCM$trtment<-factor(TDCM$trtment)
TDCM$stop<-TDCM$end

## Model fitting

set.seed(868)
start_time <- Sys.time()
m1.1<-coxph(Surv(start, stop, serinf) ~ ovmean + dascore + trtment + pgen + age + 
            bmi + renal + lung + diabetes + smoke + steroids + 
            previous_dmards + disdur, data=TDCM)
end_time <- Sys.time()
fit_time <- end_time - start_time

# # Temporal assessment

LT<-c(0,0.5,1,1.5,2,2.5,3)
Cslope<-0
CslopeL<-0
CslopeU<-0
Cind<-0
CindL<-0
CindU<-0
BS<-0

## LT = 0

m1coef<-as.numeric(m1.1$coefficients)

test0<-read.dta("baseline_model.dta")
test0$trtment <- test0$firsttreat
test0$smoke<-factor(test0$smoke)
test0$trtment<-factor(test0$trtment)
#test0$steroids<-test0$on_steroid_at_baseline
test0$prev_sinf<-0

lpcov01<-test0 %>% select(ovmean, dascore, trt2, trt4, trt1561, trt10060, 
                          pgen, age, bmi, renal, lung, diabetes, smoke1, smoke2,
                          steroids, previous_dmards, disdur)
lpmat01<-as.matrix(lpcov01)
lp0<-m1coef %*% t(lpmat01)
lp0<-as.vector(t(lp0))

start_time0 <- Sys.time()
risk = function(model, newdata, time) {
  as.numeric(1-summary(survfit(model, newdata = newdata, 
                               se.fit = F, conf.int = F), 
                       times = time)$surv)
}
test0$surv<-1-risk(m1.1, newdata=test0, time = 1)
end_time0 <- Sys.time()
pred_time0 <- end_time0 - start_time0

## Calibration slope - 95% Confidence Interval

calslope<-coxph(Surv(timevent, event)~lp0, data = test0)
x0<-confint(calslope)
Cslope[1]<-calslope$coefficients[1]
CslopeL[1]<-x0[1]
CslopeU[1]<-x0[2]

## Brier score - Kaplan Meier IPCW

BS[1]<-sbrier(obj = Surv(test0$timevent,test0$event), 
            test0$surv, btime = 1)

## Harrell's C-index - 95% confidence interval

Cind[1]<-summary(calslope)$concordance[1]
CindL[1]<-summary(calslope)$concordance[1]-
  (1.96*summary(calslope)$concordance[2])
CindU[1]<-summary(calslope)$concordance[1]+
  (1.96*summary(calslope)$concordance[2])

## Calibration plot

val_ests <- val.surv(est.surv = test0$surv, 
                     S = Surv(test0$timevent, test0$event), 
                     u=1, fun=function(p)log(-log(p)), 
                     pred = sort(runif(100, 0.7, 1)))
par(pty="s")
plot(val_ests,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
groupkm(test0$surv, S = Surv(test0$timevent,test0$event), 
        g=10,u=1, pl=T, add=T,lty=0, col=cols[7], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),bty="n",
       col=c(cols[7], "black", "black"))

## LT= 0.5

test1<-read.dta("validation_6months.dta")
test1$trtment <- test1$firsttreat
test1$smoke<-factor(test1$smoke)
test1$trtment<-factor(test1$trtment)

lpcov11<-test1 %>% select(ovmean, dascore, trt2, trt4, trt1561, trt10060, 
                          pgen, age, bmi, renal, lung, diabetes, smoke1, smoke2,
                          steroids, previous_dmards, disdur)
lpmat11<-as.matrix(lpcov11)
lp1<-m1coef %*% t(lpmat11)
lp1<-as.vector(t(lp1))

start_time1 <- Sys.time()
test1$surv<-1-risk(m1.1, newdata=test1, time = 1.5)
end_time1 <- Sys.time()
pred_time1 <- end_time1 - start_time1

## Calibration slope - 95% Confidence Interval

calslope1<-coxph(Surv(timevent, event)~lp1, data = test1)
x1<-confint(calslope1)
Cslope[2]<-calslope1$coefficients[1]
CslopeL[2]<-x1[1]
CslopeU[2]<-x1[2]

## Brier score - Kaplan Meier IPCW
BS[2]<-sbrier(obj = Surv(test1$timevent,test1$event), 
              test1$surv, btime = 1.5)

## Harrell's C-index - 95% confidence interval

Cind[2]<-summary(calslope1)$concordance[1]
CindL[2]<-summary(calslope1)$concordance[1]-
  (1.96*summary(calslope1)$concordance[2])
CindU[2]<-summary(calslope1)$concordance[1]+
  (1.96*summary(calslope1)$concordance[2])

## Calibration plot

val_ests <- val.surv(est.surv = test1$surv, 
                     S = Surv(test1$timevent, test1$event), 
                     u=1.5, fun=function(p)log(-log(p)), 
                     pred = sort(runif(100, 0.6, 1)))
plot(val_ests,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
groupkm(test1$surv, S = Surv(test1$timevent,test1$event), 
        g=10,u=1.5, pl=T, add=T,lty=0, col=cols[7], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),bty="n",
       col=c(cols[7], "black", "black"))

## LT= 1

test2<-read.dta("validation_12months.dta")
test2$trtment <- test2$firsttreat
test2$smoke<-factor(test2$smoke)
test2$trtment<-factor(test2$trtment)

lpcov21<-test2 %>% select(ovmean, dascore, trt2, trt4, trt1561, trt10060, 
                          pgen, age, bmi, renal, lung, diabetes, smoke1, smoke2,
                          steroids, previous_dmards, disdur)

lpmat21<-as.matrix(lpcov21)
lp2<-m1coef %*% t(lpmat21)
lp2<-as.vector(t(lp2))

start_time2 <- Sys.time()
test2$surv<-1-risk(m1.1, newdata=test2, time = 2)
end_time2 <- Sys.time()
pred_time2 <- end_time2 - start_time2

## Calibration slope - 95% Confidence Interval

calslope2<-coxph(Surv(timevent, event)~lp2, data = test2)
x2<-confint(calslope2)
Cslope[3]<-calslope2$coefficients[1]
CslopeL[3]<-x2[1]
CslopeU[3]<-x2[2]

## Brier score - Kaplan Meier IPCW
BS[3]<-sbrier(obj = Surv(test2$timevent,test2$event), 
              test2$surv, btime = 2)

## Harrell's C-index - 95% confidence interval

Cind[3]<-summary(calslope2)$concordance[1]
CindL[3]<-summary(calslope2)$concordance[1]-
  (1.96*summary(calslope2)$concordance[2])
CindU[3]<-summary(calslope2)$concordance[1]+
  (1.96*summary(calslope2)$concordance[2])

## Calibration plot

val_ests <- val.surv(est.surv = test2$surv, 
                     S = Surv(test2$timevent, test2$event), 
                     u=2, fun=function(p)log(-log(p)), 
                     pred = sort(runif(100, 0.6, 1)))
plot(val_ests,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
groupkm(test2$surv, S = Surv(test2$timevent,test2$event), 
        g=10,u=2, pl=T, add=T,lty=0, col=cols[7], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),bty="n",
       col=c(cols[7], "black", "black"))

## LT= 1.5

test3<-read.dta("validation_18months.dta")
test3$trtment <- test3$firsttreat
test3$smoke<-factor(test3$smoke)
test3$trtment<-factor(test3$trtment)

lpcov31<-test3 %>% select(ovmean, dascore, trt2, trt4, trt1561, trt10060, 
                          pgen, age, bmi, renal, lung, diabetes, smoke1, smoke2,
                          steroids, previous_dmards, disdur)
lpmat31<-as.matrix(lpcov31)
lp3<-m1coef %*% t(lpmat31)
lp3<-as.vector(t(lp3))

start_time3 <- Sys.time()
test3$surv<-1-risk(m1.1, newdata=test3, time = 2.5)
end_time3 <- Sys.time()
pred_time3 <- end_time3 - start_time3

## Calibration slope - 95% Confidence Interval

calslope3<-coxph(Surv(timevent, event)~lp3, data = test3)
x3<-confint(calslope3)
Cslope[4]<-calslope3$coefficients[1]
CslopeL[4]<-x3[1]
CslopeU[4]<-x3[2]

## Brier score - Kaplan Meier IPCW
BS[4]<-sbrier(obj = Surv(test3$timevent,test3$event), 
              test3$surv, btime = 2.5)

## Harrell's C-index - 95% confidence interval

Cind[4]<-summary(calslope3)$concordance[1]
CindL[4]<-summary(calslope3)$concordance[1]-
  (1.96*summary(calslope3)$concordance[2])
CindU[4]<-summary(calslope3)$concordance[1]+
  (1.96*summary(calslope3)$concordance[2])

## Calibration plot

val_ests <- val.surv(est.surv = test3$surv, 
                     S = Surv(test3$timevent, test3$event), 
                     u=2.5, fun=function(p)log(-log(p)), 
                     pred = sort(runif(100, 0.6, 1)))
plot(val_ests,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
groupkm(test3$surv, S = Surv(test3$timevent,test3$event), 
        g=10,u=2.5, pl=T, add=T,lty=0, col=cols[7], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),bty="n",
       col=c(cols[7], "black", "black"))

## LT = 2

test4<-read.dta("validation_24months.dta")
test4$trtment <- test4$firsttreat
test4$smoke<-factor(test4$smoke)
test4$trtment<-factor(test4$trtment)

lpcov41<-test4 %>% select(ovmean, dascore, trt2, trt4, trt1561, trt10060, 
                          pgen, age, bmi, renal, lung, diabetes, smoke1, smoke2,
                          steroids, previous_dmards, disdur)
lpmat41 <- as.matrix(lpcov41)                          
lp4<-m1coef %*% t(lpmat41)
lp4<-as.vector(t(lp4))

start_time4 <- Sys.time()
test4$surv<-1-risk(m1.1, newdata=test4, time = 3)
end_time4 <- Sys.time()
pred_time4 <- end_time4 - start_time4

## Calibration slope - 95% Confidence Interval

calslope4<-coxph(Surv(timevent, event)~lp4, data = test4)
x4<-confint(calslope4)
Cslope[5]<-calslope4$coefficients[1]
CslopeL[5]<-x4[1]
CslopeU[5]<-x4[2]

## Brier score - Kaplan Meier IPCW
BS[5]<-sbrier(obj = Surv(test4$timevent,test4$event), 
              test4$surv, btime = 3)

## Harrell's C-index - 95% confidence interval

Cind[5]<-summary(calslope4)$concordance[1]
CindL[5]<-summary(calslope4)$concordance[1]-
  (1.96*summary(calslope4)$concordance[2])
CindU[5]<-summary(calslope4)$concordance[1]+
  (1.96*summary(calslope4)$concordance[2])

## Calibration plot

val_ests <- val.surv(est.surv = test4$surv, 
                     S = Surv(test4$timevent, test4$event), 
                     u=3, fun=function(p)log(-log(p)), 
                     pred = sort(runif(100, 0.6, 1)))
plot(val_ests,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
groupkm(test4$surv, S = Surv(test4$timevent,test4$event), 
        g=10,u=3, pl=T, add=T,lty=0, col=cols[7], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),bty="n",
       col=c(cols[7], "black", "black"))

## LT = 2.5

test5<-read.dta("validation_30months.dta")
test5$trtment <- test5$firsttreat
test5$smoke<-factor(test5$smoke)
test5$trtment<-factor(test5$trtment)

lpcov51<-test5 %>% select(ovmean, dascore, trt2, trt4, trt1561, trt10060, 
                          pgen, age, bmi, renal, lung, diabetes, smoke1, smoke2,
                          steroids, previous_dmards, disdur)
lpmat51<-as.matrix(lpcov51)
lp5<-m1coef %*% t(lpmat51)
lp5<-as.vector(t(lp5))

start_time5 <- Sys.time()
test5$surv<-1-risk(m1.1, newdata=test5, time = 3.5)
end_time5 <- Sys.time()
pred_time5 <- end_time5 - start_time5

## Calibration slope - 95% Confidence Interval

calslope5<-coxph(Surv(timevent, event)~lp5, data = test5)
x5<-confint(calslope5)
Cslope[6]<-calslope5$coefficients[1]
CslopeL[6]<-x5[1]
CslopeU[6]<-x5[2]

## Brier score - Kaplan Meier IPCW

BS[6]<-sbrier(obj = Surv(test5$timevent,test5$event), 
              test5$surv, btime = 3.5)

## Harrell's C-index - 95% confidence interval

Cind[6]<-summary(calslope5)$concordance[1]
CindL[6]<-summary(calslope5)$concordance[1]-
  (1.96*summary(calslope5)$concordance[2])
CindU[6]<-summary(calslope5)$concordance[1]+
  (1.96*summary(calslope5)$concordance[2])

## Calibration plot

val_ests <- val.surv(est.surv = test5$surv, 
                     S = Surv(test5$timevent, test5$event), 
                     u=3.5, fun=function(p)log(-log(p)), 
                     pred = sort(runif(100, 0.6, 1)))
plot(val_ests,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
groupkm(test5$surv, S = Surv(test5$timevent,test5$event), 
        g=10,u=3.5, pl=T, add=T,lty=0, col=cols[7], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),bty="n",
       col=c(cols[7], "black", "black"))

## LT = 3

test6<-read.dta("validation_36months.dta")
test6$trtment <- test6$firsttreat
test6$smoke<-factor(test6$smoke)
test6$trtment<-factor(test6$trtment)

lpcov61<-test6 %>% select(ovmean, dascore, trt2, trt4, trt1561, trt10060, 
                          pgen, age, bmi, renal, lung, diabetes, smoke1, smoke2,
                          steroids, previous_dmards, disdur)
lpmat61<-as.matrix(lpcov61)
lp6<-m1coef %*% t(lpmat61)
lp6<-as.vector(t(lp6))

start_time6 <- Sys.time()
test6$surv<-1-risk(m1.1, newdata=test6, time = 4)
end_time6 <- Sys.time()
pred_time6 <- end_time6 - start_time6
## Calibration slope - 95% Confidence Interval

calslope6<-coxph(Surv(timevent, event)~lp6, data = test6)
x6<-confint(calslope6)
Cslope[7]<-calslope6$coefficients[1]
CslopeL[7]<-x6[1]
CslopeU[7]<-x6[2]

## Brier score - Kaplan Meier IPCW
BS[7]<-sbrier(obj = Surv(test6$timevent,test6$event), 
              test6$surv, btime = 4)

## Harrell's C-index - 95% confidence interval

Cind[7]<-summary(calslope6)$concordance[1]
CindL[7]<-summary(calslope6)$concordance[1]-
  (1.96*summary(calslope6)$concordance[2])
CindU[7]<-summary(calslope6)$concordance[1]+
  (1.96*summary(calslope6)$concordance[2])

## Calibration plot

val_ests <- val.surv(est.surv = test6$surv, 
                     S = Surv(test6$timevent, test6$event), 
                     u=4, fun=function(p)log(-log(p)), 
                     pred = sort(runif(100, 0.6, 1)))
plot(val_ests,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
groupkm(test6$surv, S = Surv(test6$timevent,test6$event), 
        g=10,u=4, pl=T, add=T,lty=0, col=cols[7], 
        cex.subtitle=FALSE)
legend("bottomright",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),bty="n",
       col=c(cols[7], "black", "black"))

pred_times <- c(pred_time0, pred_time1, pred_time2, pred_time3, pred_time4, pred_time5, pred_time6)
fit_times <- rep(fit_time, times = 7)
  
TDCM.TA<-data.frame(LT, Cslope, CslopeL, CslopeU, Cind, CindL, CindU, BS, fit_times, pred_times)
write.table(TDCM.TA, file = "tdcmTA.csv", sep = ",", col.names = NA,
            qmethod = "double")
x<-read.table("tdcmTA.csv", header = TRUE, sep = ",", row.names = 1)
head(x)

# # Bootstrap function

bootvadTDCM<-function(data, B, seed){
  dat <- data
  dat$studyno <- dat$groupid
  indiv <- unique(dat$studyno)
  C<-0
  od.C<-0
  c<-matrix(0, nrow = B, ncol = 7)
  bs<-matrix(0,nrow = B, ncol = 7)
  od.c<-matrix(0, nrow = B, ncol = 7)
  bs.od<-matrix(0, nrow = B, ncol = 7)
  opt.d<-matrix(0, nrow = B, ncol = 7)
  opt.c<-matrix(0, nrow = B, ncol = 7)
  results<-matrix(0, nrow = 2, ncol = 7)
  k10 <- qchisq(0.1,1,lower.tail=FALSE)
  
  # original data - temporal assessment
  
  od0 <- test0
  od1 <- test1
  od2 <- test2
  od3 <- test3
  od4 <- test4
  od5 <- test5
  od6 <- test6
  print(names(od0))
  
  od0$studyno<-od0$groupid
  od1$studyno<-od1$groupid
  od2$studyno<-od2$groupid
  od3$studyno<-od3$groupid
  od4$studyno<-od4$groupid
  od5$studyno<-od5$groupid
  od6$studyno<-od6$groupid
  
  for (j in 1:B){
    seed<-seed+j
    set.seed(seed)
    smp <- sort(sample(indiv, length(indiv), replace=TRUE))
    smp.df <- data.frame(studyno=smp)

    # renaming duplicated ID numbers
    
    smp.df<-smp.df %>% group_by(studyno) %>% mutate(count = row_number())
    smp.df$count<-seq(1, nrow(smp.df), by = 1)
    smp.df$newid<-paste(smp.df$studyno, smp.df$count, sep=".")
    dat.b = merge(smp.df, dat, by = "studyno", all.x=TRUE)
    dat.b$studyno <- dat.b$newid
    dat.b<-dat.b[with(dat.b, order(newid, start)),]
    
    # fitting TDCM
    
    tdcm.b <- coxph(Surv(start, stop, serinf) ~ ovmean + trtment + pgen + age + 
                bmi + renal + lung + diabetes + smoke + steroids + 
                previous_dmards + disdur, data=dat.b)
    
    # Temporal assessment sets

    val0 <- na.omit(merge(smp.df, od0, by = "studyno", all.x=TRUE))
    val1 <- na.omit(merge(smp.df, od1, by = "studyno", all.x=TRUE))
    val2 <- na.omit(merge(smp.df, od2, by = "studyno", all.x=TRUE))
    val3 <- na.omit(merge(smp.df, od3, by = "studyno", all.x=TRUE))
    val4 <- na.omit(merge(smp.df, od4, by = "studyno", all.x=TRUE))
    val5 <- na.omit(merge(smp.df, od5, by = "studyno", all.x=TRUE))
    val6 <- na.omit(merge(smp.df, od6, by = "studyno", all.x=TRUE))

    # Survival predictions at different time points
    
    val0$surv<-1-risk(tdcm.b, newdata=val0, time = 1)
    val1$surv<-1-risk(tdcm.b, newdata=val1, time = 1.5)
    val2$surv<-1-risk(tdcm.b, newdata=val2, time = 2)
    val3$surv<-1-risk(tdcm.b, newdata=val3, time = 2.5)
    val4$surv<-1-risk(tdcm.b, newdata=val4, time = 3)
    val5$surv<-1-risk(tdcm.b, newdata=val5, time = 3.5)
    val6$surv<-1-risk(tdcm.b, newdata=val6, time = 4)
    
    # Discrimination measure
    
    C[1]<- rcorr.cens(val0$surv, 
                      Surv(val0$timevent, val0$event))["C Index"]
    C[2]<- rcorr.cens(val1$surv, 
                      Surv(val1$timevent, val1$event))["C Index"]
    C[3]<- rcorr.cens(val2$surv, 
                      Surv(val2$timevent, val2$event))["C Index"]
    C[4]<- rcorr.cens(val3$surv, 
                      Surv(val3$timevent, val3$event))["C Index"]
    C[5]<- rcorr.cens(val4$surv, 
                      Surv(val4$timevent, val4$event))["C Index"]
    C[6]<- rcorr.cens(val5$surv, 
                      Surv(val5$timevent, val5$event))["C Index"]
    C[7]<- rcorr.cens(val6$surv, 
                      Surv(val6$timevent, val6$event))["C Index"]
    
    c[j,]<-C
    
    # Calibration measure
    
    bs0<-sbrier(obj = Surv(val0$timevent,val0$event), 
                val0$surv, btime = 1)
    bs1<-sbrier(obj = Surv(val1$timevent,val1$event), 
                val1$surv, btime = 1.5)
    bs2<-sbrier(obj = Surv(val2$timevent,val2$event), 
                val2$surv, btime = 2)
    bs3<-sbrier(obj = Surv(val3$timevent,val3$event), 
                val3$surv, btime = 2.5)
    bs4<-sbrier(obj = Surv(val4$timevent,val4$event), 
                val4$surv, btime = 3)
    bs5<-sbrier(obj = Surv(val5$timevent,val5$event), 
                val5$surv, btime = 3.5)
    bs6<-sbrier(obj = Surv(val6$timevent,val6$event), 
                val6$surv, btime = 4)
    
    bs[j,]<-c(bs0, bs1, bs2, bs3, bs4, bs5, bs6)
    
    # Survival predictions - original data
    
    od0$surv<-1-risk(tdcm.b, newdata=od0, time = 1)
    od1$surv<-1-risk(tdcm.b, newdata=od1, time = 1.5)
    od2$surv<-1-risk(tdcm.b, newdata=od2, time = 2)
    od3$surv<-1-risk(tdcm.b, newdata=od3, time = 2.5)
    od4$surv<-1-risk(tdcm.b, newdata=od4, time = 3)
    od5$surv<-1-risk(tdcm.b, newdata=od5, time = 3.5)
    od6$surv<-1-risk(tdcm.b, newdata=od6, time = 4)
    
    # Discrimination measure
    
    od.C[1]<- rcorr.cens(od0$surv, 
                         Surv(od0$timevent, od0$event))["C Index"]
    od.C[2]<- rcorr.cens(od1$surv, 
                         Surv(od1$timevent, od1$event))["C Index"]
    od.C[3]<- rcorr.cens(od2$surv, 
                         Surv(od2$timevent, od2$event))["C Index"]
    od.C[4]<- rcorr.cens(od3$surv, 
                         Surv(od3$timevent, od3$event))["C Index"]
    od.C[5]<- rcorr.cens(od4$surv, 
                         Surv(od4$timevent, od4$event))["C Index"]
    od.C[6]<- rcorr.cens(od5$surv, 
                         Surv(od5$timevent, od5$event))["C Index"]
    od.C[7]<- rcorr.cens(od6$surv, 
                         Surv(od6$timevent, od6$event))["C Index"]
    
    od.c[j,]<-od.C

    # Calibration measure
    
    bs0.od<-sbrier(obj = Surv(od0$timevent,od0$event), 
                   od0$surv, btime = 1)
    bs1.od<-sbrier(obj = Surv(od1$timevent,od1$event), 
                   od1$surv, btime = 1.5)
    bs2.od<-sbrier(obj = Surv(od2$timevent,od2$event), 
                   od2$surv, btime = 2)
    bs3.od<-sbrier(obj = Surv(od3$timevent,od3$event), 
                   od3$surv, btime = 2.5)
    bs4.od<-sbrier(obj = Surv(od4$timevent,od4$event), 
                   od4$surv, btime = 3)
    bs5.od<-sbrier(obj = Surv(od5$timevent,od5$event), 
                   od5$surv, btime = 3.5)
    bs6.od<-sbrier(obj = Surv(od6$timevent,od6$event), 
                   od6$surv, btime = 4)
    
    bs.od[j,]<-c(bs0.od, bs1.od, bs2.od, bs3.od, bs4.od, bs5.od, bs6.od)

    
    opt.d[j,]<-as.numeric(c[j,]-od.c[j,])
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
cl<-makeCluster(cores[1]-1) #always useful to keep one core free
start_time<-Sys.time()
registerDoParallel(cl)
TDCM.opt <- foreach(i=1:200, .combine=rbind, .packages = c("ipred", "MASS", 
                                                   "survival", "Hmisc", "dplyr")) %dopar% {
                                                     seed<-199430+i
                                                     bootvadTDCM(TDCM,1, seed = seed)
}
stopCluster(cl)
end_time<-Sys.time()
end_time - start_time

## Time difference for 5 bootstraps: 13 mins to under 5 mins
######################################################################

#Concatenating foreach loop output 

TDCM.opt.C=TDCM.opt[seq(1,399,2),]
TDCM.opt.B=TDCM.opt[seq(2,400,2),]
TDCM.opt2<-matrix(0, nrow=2, ncol=7)
TDCM.opt2[1,]=c(mean(TDCM.opt.C[,1]), mean(TDCM.opt.C[,2]),
               mean(TDCM.opt.C[,3]), mean(TDCM.opt.C[,4]),
               mean(TDCM.opt.C[,5]), mean(TDCM.opt.C[,6]),
               mean(TDCM.opt.C[,7]))
TDCM.opt2[2,]=c(mean(TDCM.opt.B[,1]), mean(TDCM.opt.B[,2]),
                mean(TDCM.opt.B[,3]), mean(TDCM.opt.B[,4]),
                mean(TDCM.opt.B[,5]), mean(TDCM.opt.B[,6]),
                mean(TDCM.opt.B[,7]))

TDCM.OPT<-data.frame(LT0 = TDCM.opt2[,1], LT0.5 = TDCM.opt2[,2], LT1 = TDCM.opt2[,3], 
                    LT1.5 = TDCM.opt2[,4], LT2 = TDCM.opt2[,5], LT2.5 = TDCM.opt2[,6], 
                    LT3 = TDCM.opt2[,7], row.names = c("C statistic (OE)", 
                                                     "Brier Score (OE)"))

write.table(TDCM.OPT, file = "tdcmOPT.csv", sep = ",", col.names = NA,
            qmethod = "double")





















## Within-method comparison

# Option 1

par(mfrow=c(1,1))
par(pty="s")
plot(TDCM$LT, TDCM$Cslope, type = "l", xlab = "Landmark time", ylab = "Calibration Slope (95% CI)", ylim = c(0.6,1.4))
lines(TDCM$LT, TDCM$CslopeU, lty = 2)
lines(TDCM$LT, TDCM$CslopeL, lty = 2)
par(new = TRUE)
plot(TDCM$LT, TDCM$Cind, type = "l", ylab ="", xlab = "", col = "red", xaxt = "n", yaxt = "n", ylim=c(0.6,0.7))
axis(side = 4)
mtext("C-index (95% CI)", side = 4, line = 3)
lines(TDCM$LT, TDCM$CindU, lty = 2, col = "red")
lines(TDCM$LT, TDCM$CindL, lty = 2, col = "red")
legend("topleft", c("C-index", "Cal. slope"),
       col = c("red", "black"), lty = c(1, 1))

# Option 2
p<- ggplot(TDCM, aes(LT, Cind))+geom_line(size = 1.5, col="red")+geom_ribbon(data=TDCM, aes(ymin=CindL, ymax=CindU), alpha = 0.2, fill = "red")
p<-p + ylab("C index")
p <- p  + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

p2<- ggplot(TDCM, aes(LT, Cslope))+geom_line(size = 1.5, col="purple")+geom_ribbon(data=TDCM, aes(ymin=CslopeL, ymax=CslopeU), alpha = 0.2, fill = "purple")
p2<-p2 + ylab("Calibration slope") + xlab("Landmark time")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(p, vp = define_region(row = 1, col = 1))
print(p2, vp = define_region(row = 2, col = 1))

#Option 3
p1<- ggplot(TDCM, aes(LT, Cind)) + geom_line(colour = "black", size = 1.5)
p1<- p1 + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "gray70", size = 0.5), panel.grid.major.x = element_blank())
p1 <- p1 + theme(panel.background = element_blank())
p1 <- p1 + scale_y_continuous(expand = c(0,0), limits = c(0.5, 1))
p1 <- p1 + theme(axis.text.y = element_text(colour="black", size = 14), axis.text.x = element_text(size = 14))
p1 <- p1 + ggtitle("C-index\n") + labs(x="Landmark time", y= NULL, size = 14)
p1 <- p1 + theme(plot.title = element_text(hjust = -0.05, vjust=1, colour="black", size = 14))
p1 <- p1 + geom_ribbon(data=TDCM, aes(ymin=CindL, ymax=CindU), alpha = 0.2)
p1 <- p1 + theme(axis.title=element_text(size=14))

p2<- ggplot(TDCM, aes(LT, Cslope)) + geom_line(colour = "red", size = 1.5)
p2<- p2 + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.grid.major.x = element_blank())
p2 <- p2 + theme(panel.background = element_blank())
p2 <- p2 + scale_y_continuous(expand = c(0,0), limits = c(0.5, 1.5))
p2 <- p2 + theme(axis.text.y = element_text(colour="red", size = 14), axis.text.x = element_text(size = 14))
p2 <- p2 + ggtitle("Calibration slope\n") + labs(x=NULL, y= NULL)
p2 <- p2 + theme(plot.title = element_text(hjust = 0.95, vjust=1, colour="red", size = 14))
p2 <- p2 + geom_ribbon(data=TDCM, aes(ymin=CslopeL, ymax=CslopeU), alpha = 0.2, fill = "red")

g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))

# get the location of the panel of p1 
# so that the panel of p2 is positioned correctly on top of it
pp <- c(subset(g1$layout, name == "panel", se = t:r))

# superimpose p2 (the panel) on p1
g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
                     pp$l, pp$b, pp$l)

# extract the y-axis of p2
ia <- which(g2$layout$name == "axis-l")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]

# flip it horizontally
ax$widths <- rev(ax$widths)
ax$grobs <- rev(ax$grobs)

# add the flipped y-axis to the right
g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

# change text content
g$grobs[[16]]$children$GRID.text.1991$label<- c("C-index", "Calibration Slope")
# change color
g$grobs[[16]]$children$GRID.text.1991$gp$col <- c("black","red")

# change x-coordinate
g$grobs[[16]]$children$GRID.text.1991$x <- unit(c(-0.05, 0.77), "npc")
grid.draw(g)










#------------------# Extended analyses
library(RColorBrewer)
cols <-brewer.pal(11, "RdBu")

TDCM <- read.dta("long_TDCM.dta")
TDCM$smoke<-factor(TDCM$smoke)
TDCM$trtment<-factor(TDCM$trtment)
TDCM$stop<-TDCM$end

temp1<-data.frame(start=0,
                  stop=3,
                  renal=0,
                  lung=0,
                  diabetes=0,
                  pgen=0,
                  steroids=0,
                  trtment="1",
                  age=56.2,
                  ovmean=h0,
                  dascore=6.41,
                  bmi = 27.4,
                  prev_sinf =0,
                  previous_dmards=3.708045,
                  smoke="0",
                  disdur=11.81012,
                  serinf=0)
temp1$trtment<-factor(temp1$trtment, ordered = F)


d0 <- TDCM %>% filter(pyears < 0.5)
d1 <- TDCM %>% filter(pyears >= 0.5 & pyears < 1)
d2 <- TDCM %>% filter(pyears >= 1 & pyears < 1.5)
d3 <- TDCM %>% filter(pyears >= 1.5 & pyears < 2)
d4 <- TDCM %>% filter(pyears >= 2 & pyears < 2.5)
d5 <- TDCM %>% filter(pyears >= 2.5 & pyears < 3)

h0 <- mean(d0$ovmean)
h1 <- mean(d1$ovmean)
h2 <- mean(d2$ovmean)
h3 <- mean(d3$ovmean)
h4 <- mean(d4$ovmean)
h5 <- mean(d5$ovmean)

temp<-data.frame(start=c(0,0.5,1,1.5,2,2.5),
                 stop=c(0.5,1,1.5,2,2.5,3),
                 renal=c(0,0,0,0,0,0),
                 lung=c(0,0,0,0,0,0),
                 diabetes=c(0,0,0,0,0,0),
                 pgen=c(0,0,0,0,0,0),
                 steroids=c(0,0,0,0,0,0),
                 trtment=c("1","1","1","1","1","1"),
                 age=rep(56.2, times=6),
                 ovmean=c(h0,h1,h2,h3,h4,h5),
                 dascore=rep(6.41, times=6),
                 bmi = rep(27.4, times=6),
                 prev_sinf = rep(0, times=6),
                 previous_dmards=c(3.708045, 3.708045, 3.708045, 3.708045, 3.708045, 3.708045),
                 smoke=c("0","0","0","0","0","0"),
                 disdur=c(11.81012,11.81012,11.81012,11.81012,11.81012,11.81012),
                 serinf=c(0,0,0,0,0,0))
temp$trtment<-factor(temp$trtment, ordered = F)

m1<-coxph(Surv(start, stop, serinf) ~ ovmean + dascore + trtment + pgen + age + 
            bmi + renal + lung + diabetes + smoke + steroids + 
            previous_dmards + disdur + prev_sinf, data=TDCM)

sdata<-read.dta("baseline_model2.dta")
sdata$trtment <- as.factor(sdata$trtment)
sdata$smoke <- as.factor(sdata$smoke)
cox<-coxph(Surv(timevent, event) ~ ovmean + dascore + trtment + pgen + age + 
             bmi + renal + lung + diabetes + smoke + steroids + 
             previous_dmards + disdur, data=sdata)

surv1<-survfit(cox, newdata=temp1, individual = T, se.fit=F)
surv1.1<-survfit(cox, newdata=temp1, individual = T)
surv2<-survfit(Surv(timevent, event)~1, sdata, se.fit=F)
surv2.1 <- survfit(Surv(timevent, event)~1, sdata)
surv3<-survfit(m1, newdata=temp, individual = T, se.fit=F)
surv3.1 <-survfit(m1, newdata=temp, individual = T)
surv4.1 <-survfit(m1, newdata=temp1, individual = T, se.fit=F)

plot(survfit(m1, newdata=temp1, individual = T), xlab="Time since baseline (years)", ylab="Suvival probability", cex.lab=1.2, lwd=2, col = alpha(cols[11], 0.3), ylim=c(0.85,1))
lines(surv2.1, lwd=2, alpha=0.8, mark.time=F, col=alpha(cols[2], 0.3))
lines(surv1.1, lwd=2, alpha=0.8, mark.time=F, col=alpha(cols[9], 0.3))
lines(surv3.1, lwd=2, alpha=0.8, mark.time=F, col=alpha(cols[4], 0.3))
lines(surv4.1, lwd=2, alpha=0.8, mark.time=F, col=alpha(cols[11], 0.9))
lines(surv2, lwd=2, alpha=0.8, mark.time=F, col=alpha(cols[2], 0.9))
lines(surv1, lwd=2, alpha=0.8, mark.time=F, col=alpha(cols[9], 0.9))
lines(surv3, lwd=2, alpha=0.8, mark.time=F, col=alpha(cols[4], 0.9))

legend("bottomleft", legend = c("Observed", "TDCM", "Cox", "Manual"), col=c(cols[2], cols[11], cols[9], cols[4]), lty=rep(1,times=4), lwd=rep(3, times=4), cex=1.2)

mage<-mean(sdata$age)
mbmi <- mean(sdata$bmi)
mdascore <- mean(sdata$dascore)



temp1<-data.frame(start=0,
                  stop=3,
                  renal=0,
                  lung=0,
                  diabetes=0,
                  pgen=0,
                  steroids=0,
                  trtment="1",
                  age=56.2,
                  ovmean=h0,
                  dascore=6.41,
                  bmi = 27.4,
                  prev_sinf =0,
                  previous_dmards=3.708045,
                  smoke="0",
                  disdur=11.81012,
                  serinf=0)
temp1$trtment<-factor(temp1$trtment, ordered = F)
temp2<-data.frame(start=c(0,0.5,1,1.5,2,2.5),
                 stop=c(0.5,1,1.5,2,2.5,3),
                 renal=c(0,0,0,0,0,0),
                 lung=c(0,0,0,0,0,0),
                 diabetes=c(0,0,0,0,0,0),
                 pgen=c(0,0,0,0,0,0),
                 steroids=c(0,0,0,0,0,0),
                 trtment=c("1","1","1","1","1","1"),
                 age=rep(56.2, times=6),
                 ovmean=c(h0,h1,h2,h3,h4,h5),
                 dascore=rep(6.41, times=6),
                 bmi = rep(27.4, times=6),
                 prev_sinf = rep(0, times=6),
                 previous_dmards=c(3.708045, 3.708045, 3.708045, 3.708045, 3.708045, 3.708045),
                 smoke=c("0","0","0","0","0","0"),
                 disdur=c(11.81012,11.81012,11.81012,11.81012,11.81012,11.81012),
                 serinf=rep(0, times=6))
temp2$trtment<-factor(temp$trtment, ordered = F)

#baseline survival calculation (time = 1.5)

bh<-basehaz(m1.1)
bh1.5<-bh[which(bh$time<1.5),]
x<-length(bh1.5$time)
bh1.5<-bh1.5$hazard[x]
bs1.5<-exp(-bh1.5)

test0$survprob<-bs1^(exp(lp0))
m2 <- cph(Surv(start, stop, serinf) ~  ovmean + trtment + pgen + age + bmi + renal + lung + diabetes + smoke + on_steroid_at_baseline + previous_dmards + disdur, data=example,surv = T,x=T,y=T)
m2.1 <- cph(Surv(start, stop, serinf) ~  ovmean + trtment + pgen + age + lung + diabetes + smoke + on_steroid_at_baseline + previous_dmards + disdur, data=example,surv = T,x=T,y=T)

## Baseline assessment

test.data<-read.dta("baseline_model.dta")
test.data$trtment <- test.data$firsttreat
test.data$smoke<-factor(test.data$smoke, levels=c("0","1","2"), ordered=T)
test.data$trtment<-factor(test.data$trtment, ordered = F)
test.data$prob = survest(m2.1,newdata=test.data,times=1)$surv
min(test.data$prob)

## Calibration plot

attach(test.data)
val_ests <- val.surv(est.surv = test.data$prob, S = Surv(test.data$timevent, test.data$event), u=1, fun=function(p)log(-log(p)), pred = sort(runif(100, 0.7, 1)))
par(pty="s")
plot(val_ests,xlab="Expected Survival Probability",ylab="Observed Survival Probability") 
groupkm(test.data$prob, S = Surv(test.data$timevent,test.data$event), g=10,u=1, pl=T, add=T,lty=0, col="purple", cex.subtitle=FALSE)
legend(0.7,1,c("Risk groups","Reference line","95% CI"),lty=c(0,2,1),pch=c(19,NA,NA),bty="n")

## 6 months assessment

test.data<-read.dta("validation_6months.dta")
test.data$trtment <- test.data$firsttreat
test.data$smoke<-factor(test.data$smoke, levels=c("0","1","2"), ordered=T)
test.data$trtment<-factor(test.data$trtment, ordered = F)
test.data$prob = survest(m2.1,newdata=test.data,times=1.5)$surv

## Calibration plot

attach(test.data)
val_ests <- val.surv(est.surv = test.data$prob, S = Surv(test.data$timevent, test.data$event), u=1.5, fun=function(p)log(-log(p)), pred = sort(runif(100, 0.7, 1)))

par(pty="s")
plot(val_ests,xlab="Expected Survival Probability",ylab="Observed Survival Probability") 
groupkm(test.data$prob, S = Surv(test.data$timevent,test.data$event), g=10,u=1.5, pl=T, add=T,lty=0, col="purple", cex.subtitle=FALSE)
legend(0.7,1,c("Risk groups","Reference line","95% CI"),lty=c(0,2,1),pch=c(19,NA,NA),bty="n")


## for (i in (2:nrow(example))){
##   if (example$studyno[i-1] == example$studyno[(i)]){example$stop[i-1] <- example$pyears[(i)] }}
## for (i in (2:nrow(example))){
##  if (example$studyno[i-1] != example$studyno[i]){example$stop[i-1] = example$stop[i-2]+0.25}}
## example$stop[nrow(example)]=example$stop[nrow(example)-1]+0.25
## example$start <- example$pyears
