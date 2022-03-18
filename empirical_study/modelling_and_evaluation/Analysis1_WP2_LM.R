library(dynpred)
library(foreign)
library(rms)
library(colorRamps)
library(grDevices)

# Data preparation

LMdata<-read.dta("long_LM.dta")
LMdata$smoke<-factor(LMdata$smoke)
LMdata$trtment<-factor(LMdata$trtment)

# Setting landmark time points

w <- 1
LMone <- seq(0,max(LMdata$pyears)-w,by=0.1)
LMtwo <- LMdata$pyears[which(LMdata$serinf==1)]
LMtwo <- sort(LMtwo, decreasing = FALSE)
LMtwo <- round(LMtwo, digits = 1)
LMtwo <- unique(LMtwo)
LMthree <- seq(0,max(LMdata$pyears)-w,by=0.5)

# # Only using the first 3 years of predictor info - HAQ isn't measured afterwards

LMone <- LMone[which(LMone <= 3)]
LMtwo <- LMtwo[which(LMtwo <= 3)]
LMthree <- LMthree[which(LMthree <= 3)]

# Version 1 - Data Prep

SLMdata <- cutLM(data=LMdata,
                outcome=list(time="timevent", status="event"),
                LM=0,
                horizon=1,
                covs=list(fixed=c("age","pgen", "previous_dmards", "disdur", "trtment", "dascore", 
                                  "steroids", "bmi", "smoke", "renal", "lung", "diabetes"), varying="ovmean"), 
                format = "long", 
                id = "studyno", 
                rtime = "pyears", 
                right = FALSE)

for (i in 2:length(LMone)) {
  LM <- LMone[i]
  SLMdata <- rbind(SLMdata,cutLM(data=LMdata,
                               outcome=list(time="timevent",status="event"),
                               LM=LM,
                               horizon=LM+w,
                               covs=list(fixed=c("age","pgen",  "previous_dmards", "disdur", "trtment"),
                                         varying=c("ovmean", "dascore", "steroids", "bmi", "smoke", "renal", "lung", "diabetes")), 
                               format = "long", 
                               id = "studyno", 
                               rtime = "pyears", 
                               right = "FALSE"))
}

# Version 2 - Data Prep

SLMdata.2 <- cutLM(data=LMdata,
                 outcome=list(time="timevent", status="event"),
                 LM=0,
                 horizon=w,
                 covs=list(fixed=c("age","pgen", "previous_dmards", "disdur", "trtment", "dascore", 
                                   "steroids", "bmi", "smoke", "renal", "lung", "diabetes"), varying="ovmean"), 
                 format = "long", 
                 id = "studyno", 
                 rtime = "pyears", 
                 right = "FALSE")

for (i in 2:length(LMtwo)) {
  LM <- LMtwo[i]
  SLMdata.2 <- rbind(SLMdata.2,cutLM(data=LMdata,
                                 outcome=list(time="timevent",status="event"),
                                 LM=LM,
                                 horizon=LM+w,
                                 covs=list(fixed=c("age","pgen",  "previous_dmards", "disdur", "trtment"),
                                           varying=c("ovmean", "dascore", "steroids", "bmi", "smoke", "renal", "lung", "diabetes")), 
                                 format = "long", 
                                 id = "studyno", 
                                 rtime = "pyears", 
                                 right = "FALSE"))
}

# Version 3 - Data Prep

SLMdata.3 <- cutLM(data=LMdata,
                   outcome=list(time="timevent", status="event"),
                   LM=0,
                   horizon=w,
                   covs=list(fixed=c("age","pgen", "previous_dmards", "disdur", "trtment", "dascore", 
                                     "steroids", "bmi", "smoke", "renal", "lung", "diabetes", "groupid"), varying="ovmean"), 
                   format = "long", 
                   id = "studyno", 
                   rtime = "pyears", 
                   right = "FALSE")

for (i in 2:length(LMthree)) {
  LM <- LMthree[i]
  SLMdata.3 <- rbind(SLMdata.3,cutLM(data=LMdata,
                                     outcome=list(time="timevent",status="event"),
                                     LM=LM,
                                     horizon=LM+w,
                                     covs=list(fixed=c("age","pgen",  "previous_dmards", "disdur", "trtment", "groupid"),
                                               varying=c("ovmean", "dascore", "steroids", "bmi", "smoke", "renal", "lung", "diabetes")), 
                                     format = "long", 
                                     id = "studyno", 
                                     rtime = "pyears", 
                                     right = "FALSE"))
}
##################################################################################
set.seed(836536753)
# CPM development
start_time1 <- Sys.time()
slm1 <- coxph(Surv(LM, timevent,event) ~ age + pgen + previous_dmards + 
                disdur + trtment + ovmean + dascore + steroids +
                bmi + smoke + renal + lung + diabetes + strata(LM) + cluster(studyno), data=SLMdata, method="breslow")
end_time1 <- Sys.time()
fit_time1 <- end_time1 - start_time1
start_time2 <- Sys.time()
slm2 <- coxph(Surv(LM,timevent,event) ~ age + pgen + previous_dmards + 
                disdur + trtment + ovmean + dascore + steroids +
                bmi + smoke + renal + lung + diabetes + strata(LM) + cluster(studyno), data=SLMdata.2, method="breslow")
end_time2 <- Sys.time()
fit_time2 <- end_time2 - start_time2
start_time3 <- Sys.time()
slm3 <- coxph(Surv(LM,timevent,event) ~ age + pgen + previous_dmards + 
                disdur + trtment + ovmean + dascore + steroids +
                bmi + smoke + renal + lung + diabetes + strata(LM) + cluster(studyno), data=SLMdata.3, method="breslow")
end_time3 <- Sys.time()
fit_time3 <- end_time3 - start_time3
# The models

summary(slm1)
confint.default(slm1)
summary(slm2)
confint.default(slm2)
summary(slm3)
confint.default(slm3)

############################# Temporal assessment ################################

Cslope1<-0
Cslope2<-0
Cslope3<-0
Cslope1L<-0
Cslope2L<-0
Cslope3L<-0
Cslope1U<-0
Cslope2U<-0
Cslope3U<-0

Cind1<-0
Cind2<-0
Cind3<-0
Cind1L<-0
Cind2L<-0
Cind3L<-0
Cind1U<-0
Cind2U<-0
Cind3U<-0

BS1<-0
BS2<-0
BS3<-0

# LT = 0  

test0 <- read.dta("baseline_model.dta")
test0$smoke<-factor(test0$smoke)
test0$trtment<-factor(test0$firsttreat)
#test0$steroids<-test0$on_steroid_at_baseline
test0$LM<-0
test0$prev_sinf<-0
test0p<-test0
test0p$timevent<-rep(1,nrow(test0p))

# Predicted probabilities and linear predictor
start_time10 <- Sys.time()
P10<-predict(slm1, newdata = test0p, type = "surv", se.fit =TRUE)
end_time10 <- Sys.time()
pred_time10 <- end_time10 - start_time10

LP10<-predict(slm1, newdata = test0p, type = "lp", se.fit =TRUE)
start_time20 <- Sys.time()
P20<-predict(slm2, newdata = test0p, type = "surv", se.fit =TRUE)
end_time20 <- Sys.time()
pred_time20 <- end_time20 - start_time20
LP20<-predict(slm2, newdata = test0p, type = "lp", se.fit =TRUE)

start_time30 <- Sys.time()
P30<-predict(slm3, newdata = test0p, type = "surv", se.fit =TRUE)
end_time30 <- Sys.time()
pred_time30 <- end_time30 - start_time30
LP30<-predict(slm3, newdata = test0p, type = "lp", se.fit =TRUE)

test0$LP1<-as.numeric(LP10$fit)
test0$LP2<-as.numeric(LP20$fit)
test0$LP3<-as.numeric(LP30$fit)

test0$P1<-as.numeric(P10$fit)
test0$P2<-as.numeric(P20$fit)
test0$P3<-as.numeric(P30$fit)

#Calibration slope

surv10 <- coxph(Surv(timevent,event)~LP1, method="breslow", data = test0)
x10<-confint(surv10)
Cslope1[1]<-surv10$coefficients[1]
Cslope1L[1]<-x10[1]
Cslope1U[1]<-x10[2]

surv20 <- coxph(Surv(timevent,event)~LP2, method="breslow", data = test0)
x20<-confint(surv20)
Cslope2[1]<-surv20$coefficients[1]
Cslope2L[1]<-x20[1]
Cslope2U[1]<-x20[2]

surv30 <- coxph(Surv(timevent,event)~LP3, method="breslow", data = test0)
x30<-confint(surv30)
Cslope3[1]<-surv30$coefficients[1]
Cslope3L[1]<-x30[1]
Cslope3U[1]<-x30[2]

# Harrell's C-statistic

Cind1[1]<-summary(surv10)$concordance[1]
Cind1L[1]<-summary(surv10)$concordance[1]-
  (1.96*summary(surv10)$concordance[2])
Cind1U[1]<-summary(surv10)$concordance[1]+
  (1.96*summary(surv10)$concordance[2])

Cind2[1]<-summary(surv20)$concordance[1]
Cind2L[1]<-summary(surv20)$concordance[1]-
  (1.96*summary(surv20)$concordance[2])
Cind2U[1]<-summary(surv20)$concordance[1]+
  (1.96*summary(surv20)$concordance[2])

Cind3[1]<-summary(surv30)$concordance[1]
Cind3L[1]<-summary(surv30)$concordance[1]-
  (1.96*summary(surv30)$concordance[2])
Cind3U[1]<-summary(surv30)$concordance[1]+
  (1.96*summary(surv30)$concordance[2])

# Brier score

BS1[1]<-sbrier(obj = Surv(test0$timevent,test0$event), 
              test0$P1, btime = 1)
BS2[1]<-sbrier(obj = Surv(test0$timevent,test0$event), 
               test0$P2, btime = 1)
BS3[1]<-sbrier(obj = Surv(test0$timevent,test0$event), 
               test0$P3, btime = 1)

# Calibration plot

par(pty="s")
par(mfrow=c(1,3))
val_ests1 <- val.surv(est.surv = P10$fit, S = Surv(test0$timevent, 
                                                   test0$event), 
                      u=1, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests1,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P10$fit), S = Surv(test0$timevent, test0$event), 
           g=10,u=1, pl=T, add=T,lty=0, col=cols[9], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),col=c(cols[9], "black", "black"), bty="n")

val_ests2 <- val.surv(est.surv = P20$fit, S = Surv(test0$timevent, 
                                                   test0$event), 
                      u=1, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests2,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P20$fit), S = Surv(test0$timevent, test0$event), 
           g=10,u=1, pl=T, add=T,lty=0, col=cols[10], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),col=c(cols[10], "black", "black"), bty="n")

val_ests3 <- val.surv(est.surv = P30$fit, S = Surv(test0$timevent, 
                                                   test0$event), 
                      u=1, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests3,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P30$fit), S = Surv(test0$timevent, test0$event), 
           g=10,u=1, pl=T, add=T,lty=0, col=cols[11], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),col=c(cols[11], "black", "black"), bty="n")


# LT = 0.5  

test1 <- read.dta("validation_6months.dta")
test1$smoke<-factor(test1$smoke)
test1$trtment<-factor(test1$firsttreat)

test1$LM<-0.5

test1p<-test1
test1p$timevent<-rep(1.5,nrow(test1p))

# Predicted probabilities and linear predictor
start_time11 <- Sys.time()
P11<-predict(slm1, newdata = test1p, type = "surv", se.fit =TRUE)
end_time11 <- Sys.time()
pred_time11 <- end_time11 - start_time11

LP11<-predict(slm1, newdata = test1p, type = "lp", se.fit =TRUE)
start_time21 <- Sys.time()
P21<-predict(slm2, newdata = test1p, type = "surv", se.fit =TRUE)
end_time21 <- Sys.time()
pred_time21 <- end_time21 - start_time21

LP21<-predict(slm2, newdata = test1p, type = "lp", se.fit =TRUE)
start_time31 <- Sys.time()
P31<-predict(slm3, newdata = test1p, type = "surv", se.fit =TRUE)
end_time31 <- Sys.time()
pred_time31 <- end_time31 - start_time31
LP31<-predict(slm3, newdata = test1p, type = "lp", se.fit =TRUE)

test1$LP1<-as.numeric(LP11$fit)
test1$LP2<-as.numeric(LP21$fit)
test1$LP3<-as.numeric(LP31$fit)

test1$P1<-as.numeric(P11$fit)
test1$P2<-as.numeric(P21$fit)
test1$P3<-as.numeric(P31$fit)

#Calibration slope

surv11 <- coxph(Surv(timevent,event)~LP1, method="breslow", data = test1)
x11<-confint(surv11)
Cslope1[2]<-surv11$coefficients[1]
Cslope1L[2]<-x11[1]
Cslope1U[2]<-x11[2]

surv21 <- coxph(Surv(timevent,event)~LP2, method="breslow", data = test1)
x21<-confint(surv21)
Cslope2[2]<-surv21$coefficients[1]
Cslope2L[2]<-x21[1]
Cslope2U[2]<-x21[2]

surv31 <- coxph(Surv(timevent,event)~LP3, method="breslow", data = test1)
x31<-confint(surv31)
Cslope3[2]<-surv31$coefficients[1]
Cslope3L[2]<-x31[1]
Cslope3U[2]<-x31[2]

# Harrell's C-statistic

Cind1[2]<-summary(surv11)$concordance[1]
Cind1L[2]<-summary(surv11)$concordance[1]-
  (1.96*summary(surv11)$concordance[2])
Cind1U[2]<-summary(surv11)$concordance[1]+
  (1.96*summary(surv11)$concordance[2])

Cind2[2]<-summary(surv21)$concordance[1]
Cind2L[2]<-summary(surv21)$concordance[1]-
  (1.96*summary(surv21)$concordance[2])
Cind2U[2]<-summary(surv21)$concordance[1]+
  (1.96*summary(surv21)$concordance[2])

Cind3[2]<-summary(surv31)$concordance[1]
Cind3L[2]<-summary(surv31)$concordance[1]-
  (1.96*summary(surv31)$concordance[2])
Cind3U[2]<-summary(surv31)$concordance[1]+
  (1.96*summary(surv31)$concordance[2])

# Brier score

BS1[2]<-sbrier(obj = Surv(test1$timevent,test1$event), 
               test1$P1, btime = 1.5)
BS2[2]<-sbrier(obj = Surv(test1$timevent,test1$event), 
               test1$P2, btime = 1.5)
BS3[2]<-sbrier(obj = Surv(test1$timevent,test1$event), 
               test1$P3, btime = 1.5)

# Calibration plot

par(pty="s")
par(mfrow=c(1,3))
val_ests1 <- val.surv(est.surv = P11$fit, S = Surv(test1$timevent, 
                                                   test1$event), 
                      u=1.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests1,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P11$fit), S = Surv(test1$timevent, test1$event), 
           g=10,u=1.5, pl=T, add=T,lty=0, col=cols[9], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),col=c(cols[9], "black", "black"), bty="n")

val_ests2 <- val.surv(est.surv = P21$fit, S = Surv(test1$timevent, 
                                                   test1$event), 
                      u=1.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests2,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P21$fit), S = Surv(test1$timevent, test1$event), 
           g=10,u=1.5, pl=T, add=T,lty=0, col=cols[10], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[10], "black", "black"), bty="n")

val_ests3 <- val.surv(est.surv = P31$fit, S = Surv(test1$timevent, 
                                                   test1$event), 
                      u=1.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests3,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P31$fit), S = Surv(test1$timevent, test1$event), 
           g=10,u=1.5, pl=T, add=T,lty=0, col=cols[11], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA),col=c(cols[11], "black", "black"), bty="n")


# LT = 1  

test2 <- read.dta("validation_12months.dta")
test2$smoke<-factor(test2$smoke)
test2$trtment<-factor(test2$firsttreat)
test2$LM<-1
test2p<-test2
test2p$timevent<-rep(2, nrow(test2p))

# Predicted probabilities and linear predictor
start_time12<-Sys.time()
P12<-predict(slm1, newdata = test2p, type = "surv", se.fit =TRUE)
end_time12 <- Sys.time()
pred_time12 <- end_time12 - start_time12
LP12<-predict(slm1, newdata = test2p, type = "lp", se.fit =TRUE)

start_time22 <- Sys.time()
P22<-predict(slm2, newdata = test2p, type = "surv", se.fit =TRUE)
end_time22 <- Sys.time()
pred_time22 <- end_time22 - start_time22
LP22<-predict(slm2, newdata = test2p, type = "lp", se.fit =TRUE)
start_time32 <- Sys.time()
P32<-predict(slm3, newdata = test2p, type = "surv", se.fit =TRUE)
end_time32 <- Sys.time()
pred_time32 <- end_time32 - start_time32
LP32<-predict(slm3, newdata = test2p, type = "lp", se.fit =TRUE)

test2$LP1<-as.numeric(LP12$fit)
test2$LP2<-as.numeric(LP22$fit)
test2$LP3<-as.numeric(LP32$fit)

test2$P1<-as.numeric(P12$fit)
test2$P2<-as.numeric(P22$fit)
test2$P3<-as.numeric(P32$fit)

#Calibration slope

surv12 <- coxph(Surv(timevent,event)~LP1, method="breslow", data = test2)
x12<-confint(surv12)
Cslope1[3]<-surv12$coefficients[1]
Cslope1L[3]<-x12[1]
Cslope1U[3]<-x12[2]

surv22 <- coxph(Surv(timevent,event)~LP2, method="breslow", data = test2)
x22<-confint(surv22)
Cslope2[3]<-surv22$coefficients[1]
Cslope2L[3]<-x22[1]
Cslope2U[3]<-x22[2]

surv32 <- coxph(Surv(timevent,event)~LP3, method="breslow", data = test2)
x32<-confint(surv32)
Cslope3[3]<-surv32$coefficients[1]
Cslope3L[3]<-x32[1]
Cslope3U[3]<-x32[2]

# Harrell's C-statistic

Cind1[3]<-summary(surv12)$concordance[1]
Cind1L[3]<-summary(surv12)$concordance[1]-
  (1.96*summary(surv12)$concordance[2])
Cind1U[3]<-summary(surv12)$concordance[1]+
  (1.96*summary(surv12)$concordance[2])

Cind2[3]<-summary(surv22)$concordance[1]
Cind2L[3]<-summary(surv22)$concordance[1]-
  (1.96*summary(surv22)$concordance[2])
Cind2U[3]<-summary(surv22)$concordance[1]+
  (1.96*summary(surv22)$concordance[2])

Cind3[3]<-summary(surv32)$concordance[1]
Cind3L[3]<-summary(surv32)$concordance[1]-
  (1.96*summary(surv32)$concordance[2])
Cind3U[3]<-summary(surv32)$concordance[1]+
  (1.96*summary(surv32)$concordance[2])

# Brier score

BS1[3]<-sbrier(obj = Surv(test2$timevent,test2$event), 
               test2$P1, btime = 2)
BS2[3]<-sbrier(obj = Surv(test2$timevent,test2$event), 
               test2$P2, btime = 2)
BS3[3]<-sbrier(obj = Surv(test2$timevent,test2$event), 
               test2$P3, btime = 2)
# Calibration plot

par(pty="s")
par(mfrow=c(1,3))
val_ests1 <- val.surv(est.surv = P12$fit, S = Surv(test2$timevent, 
                                                   test2$event), 
                      u=2, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests1,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P12$fit), S = Surv(test2$timevent, test2$event), 
           g=10,u=2, pl=T, add=T,lty=0, col=cols[9], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA) , col=c(cols[9], "black", "black"), bty="n")

val_ests2 <- val.surv(est.surv = P22$fit, S = Surv(test2$timevent, 
                                                   test2$event), 
                      u=2, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests2,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P22$fit), S = Surv(test2$timevent, test2$event), 
           g=10,u=2, pl=T, add=T,lty=0, col=cols[10], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[10], "black", "black"), bty="n")

val_ests3 <- val.surv(est.surv = P32$fit, S = Surv(test2$timevent, 
                                                   test2$event), 
                      u=2, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests3,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P32$fit), S = Surv(test2$timevent, test2$event), 
           g=10,u=2, pl=T, add=T,lty=0, col=cols[11], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[11], "black", "black"),  bty="n")


# LT = 1.5  

test3 <- read.dta("validation_18months.dta")
test3$smoke<-factor(test3$smoke)
test3$trtment<-factor(test3$firsttreat)
test3$LM<-1.5
test3p<-test3
test3p$timevent<-rep(2.5,nrow(test3p))

# Predicted probabilities and linear predictor

start_time13 <- Sys.time()
P13<-predict(slm1, newdata = test3p, type = "surv", se.fit =TRUE)
end_time13 <- Sys.time()
pred_time13 <- end_time13 - start_time13
LP13<-predict(slm1, newdata = test3p, type = "lp", se.fit =TRUE)
start_time23 <- Sys.time()
P23<-predict(slm2, newdata = test3p, type = "surv", se.fit =TRUE)
end_time23 <- Sys.time()
pred_time23 <- end_time23 - start_time23
LP23<-predict(slm2, newdata = test3p, type = "lp", se.fit =TRUE)
start_time33 <- Sys.time()
P33<-predict(slm3, newdata = test3p, type = "surv", se.fit =TRUE)
end_time33 <- Sys.time()
pred_time33 <- end_time33 - start_time33
LP33<-predict(slm3, newdata = test3p, type = "lp", se.fit =TRUE)

test3$LP1<-as.numeric(LP13$fit)
test3$LP2<-as.numeric(LP23$fit)
test3$LP3<-as.numeric(LP33$fit)

test3$P1<-as.numeric(P13$fit)
test3$P2<-as.numeric(P23$fit)
test3$P3<-as.numeric(P33$fit)

#Calibration slope

surv13 <- coxph(Surv(timevent,event)~LP1, method="breslow", data = test3)
x13<-confint(surv13)
Cslope1[4]<-surv13$coefficients[1]
Cslope1L[4]<-x13[1]
Cslope1U[4]<-x13[2]

surv23 <- coxph(Surv(timevent,event)~LP2, method="breslow", data = test3)
x23<-confint(surv23)
Cslope2[4]<-surv23$coefficients[1]
Cslope2L[4]<-x23[1]
Cslope2U[4]<-x23[2]

surv33 <- coxph(Surv(timevent,event)~LP3, method="breslow", data = test3)
x33<-confint(surv33)
Cslope3[4]<-surv33$coefficients[1]
Cslope3L[4]<-x33[1]
Cslope3U[4]<-x33[2]

# Harrell's C-statistic

Cind1[4]<-summary(surv13)$concordance[1]
Cind1L[4]<-summary(surv13)$concordance[1]-
  (1.96*summary(surv13)$concordance[2])
Cind1U[4]<-summary(surv13)$concordance[1]+
  (1.96*summary(surv13)$concordance[2])

Cind2[4]<-summary(surv23)$concordance[1]
Cind2L[4]<-summary(surv23)$concordance[1]-
  (1.96*summary(surv23)$concordance[2])
Cind2U[4]<-summary(surv23)$concordance[1]+
  (1.96*summary(surv23)$concordance[2])

Cind3[4]<-summary(surv33)$concordance[1]
Cind3L[4]<-summary(surv33)$concordance[1]-
  (1.96*summary(surv33)$concordance[2])
Cind3U[4]<-summary(surv33)$concordance[1]+
  (1.96*summary(surv33)$concordance[2])

# Brier score

BS1[4]<-sbrier(obj = Surv(test3$timevent,test3$event), 
               test3$P1, btime = 2.5)
BS2[4]<-sbrier(obj = Surv(test3$timevent,test3$event), 
               test3$P2, btime = 2.5)
BS3[4]<-sbrier(obj = Surv(test3$timevent,test3$event), 
               test3$P3, btime = 2.5)
# Calibration plot

par(pty="s")
par(mfrow=c(1,3))
val_ests1 <- val.surv(est.surv = P13$fit, S = Surv(test3$timevent, 
                                                   test3$event), 
                      u=2.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests1,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P13$fit), S = Surv(test3$timevent, test3$event), 
           g=10,u=2.5, pl=T, add=T,lty=0, col=cols[9], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[9], "black", "black"), bty="n")

val_ests2 <- val.surv(est.surv = P23$fit, S = Surv(test3$timevent, 
                                                   test3$event), 
                      u=2.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests2,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P23$fit), S = Surv(test3$timevent, test3$event), 
           g=10,u=2.5, pl=T, add=T,lty=0, col=cols[10], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[10], "black", "black"), bty="n")

val_ests3 <- val.surv(est.surv = P33$fit, S = Surv(test3$timevent, 
                                                   test3$event), 
                      u=2.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests3,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P33$fit), S = Surv(test3$timevent, test3$event), 
           g=10,u=2.5, pl=T, add=T,lty=0, col=cols[11], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[11], "black", "black"), bty="n")


# LT = 2  

test4 <- read.dta("validation_24months.dta")
test4$smoke<-factor(test4$smoke)
test4$trtment<-factor(test4$firsttreat)
test4$LM<-2
test4p<-test4
test4p$timevent<-rep(3,nrow(test4p))

# Predicted probabilities and linear predictor
start_time14 <- Sys.time()
P14<-predict(slm1, newdata = test4p, type = "surv", se.fit =TRUE)
end_time14 <- Sys.time()
pred_time14 <- end_time14 - start_time14
LP14<-predict(slm1, newdata = test4p, type = "lp", se.fit =TRUE)

start_time24 <- Sys.time()
P24<-predict(slm2, newdata = test4p, type = "surv", se.fit =TRUE)
end_time24 <- Sys.time()
pred_time24 <- end_time24 - start_time24
LP24<-predict(slm2, newdata = test4p, type = "lp", se.fit =TRUE)

start_time34 <- Sys.time()
P34<-predict(slm3, newdata = test4p, type = "surv", se.fit =TRUE)
end_time34 <- Sys.time()
pred_time34 <- end_time34 - start_time34
LP34<-predict(slm3, newdata = test4p, type = "lp", se.fit =TRUE)

test4$LP1<-as.numeric(LP14$fit)
test4$LP2<-as.numeric(LP24$fit)
test4$LP3<-as.numeric(LP34$fit)

test4$P1<-as.numeric(P14$fit)
test4$P2<-as.numeric(P24$fit)
test4$P3<-as.numeric(P34$fit)

#Calibration slope

surv14 <- coxph(Surv(timevent,event)~LP1, method="breslow", data = test4)
x14<-confint(surv14)
Cslope1[5]<-surv14$coefficients[1]
Cslope1L[5]<-x14[1]
Cslope1U[5]<-x14[2]

surv24 <- coxph(Surv(timevent,event)~LP2, method="breslow", data = test4)
x24<-confint(surv24)
Cslope2[5]<-surv24$coefficients[1]
Cslope2L[5]<-x24[1]
Cslope2U[5]<-x24[2]

surv34 <- coxph(Surv(timevent,event)~LP3, method="breslow", data = test4)
x34<-confint(surv34)
Cslope3[5]<-surv34$coefficients[1]
Cslope3L[5]<-x34[1]
Cslope3U[5]<-x34[2]

# Harrell's C-statistic

Cind1[5]<-summary(surv14)$concordance[1]
Cind1L[5]<-summary(surv14)$concordance[1]-
  (1.96*summary(surv14)$concordance[2])
Cind1U[5]<-summary(surv14)$concordance[1]+
  (1.96*summary(surv14)$concordance[2])

Cind2[5]<-summary(surv24)$concordance[1]
Cind2L[5]<-summary(surv24)$concordance[1]-
  (1.96*summary(surv24)$concordance[2])
Cind2U[5]<-summary(surv24)$concordance[1]+
  (1.96*summary(surv24)$concordance[2])

Cind3[5]<-summary(surv34)$concordance[1]
Cind3L[5]<-summary(surv34)$concordance[1]-
  (1.96*summary(surv34)$concordance[2])
Cind3U[5]<-summary(surv34)$concordance[1]+
  (1.96*summary(surv34)$concordance[2])

# Brier score

BS1[5]<-sbrier(obj = Surv(test4$timevent,test4$event), 
               test4$P1, btime = 3)
BS2[5]<-sbrier(obj = Surv(test4$timevent,test4$event), 
               test4$P2, btime = 3)
BS3[5]<-sbrier(obj = Surv(test4$timevent,test4$event), 
               test4$P3, btime = 3)
# Calibration plot

par(pty="s")
par(mfrow=c(1,3))

val_ests1 <- val.surv(est.surv = P14$fit, S = Surv(test4$timevent, 
                                                   test4$event), 
                      u=3, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests1,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P14$fit), S = Surv(test4$timevent, test4$event), 
           g=10,u=3, pl=T, add=T,lty=0, col=cols[9], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[9], "black", "black"), bty="n")

val_ests2 <- val.surv(est.surv = P24$fit, S = Surv(test4$timevent, 
                                                   test4$event), 
                      u=3, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests2,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P24$fit), S = Surv(test4$timevent, test4$event), 
           g=10,u=3, pl=T, add=T,lty=0, col=cols[10], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[10], "black", "black"), bty="n")

val_ests3 <- val.surv(est.surv = P34$fit, S = Surv(test4$timevent, 
                                                   test4$event), 
                      u=3, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests3,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P34$fit), S = Surv(test4$timevent, test4$event), 
           g=10,u=3, pl=T, add=T,lty=0, col=cols[11], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[11], "black", "black"), bty="n")



# LT = 2.5  

test5 <- read.dta("validation_30months.dta")
test5$smoke<-factor(test5$smoke)
test5$trtment<-factor(test5$firsttreat)
test5$LM<-2.5
test5p<-test5
test5p$timevent<-rep(3.5,nrow(test5p))

# Predicted probabilities and linear predictor
start_time15 <- Sys.time()
P15<-predict(slm1, newdata = test5p, type = "surv", se.fit =TRUE)
end_time15 <- Sys.time()
pred_time15 <- end_time15 - start_time15
LP15<-predict(slm1, newdata = test5p, type = "lp", se.fit =TRUE)
start_time25 <- Sys.time()
P25<-predict(slm2, newdata = test5p, type = "surv", se.fit =TRUE)
end_time25 <- Sys.time()
pred_time25 <- end_time25 - start_time25
LP25<-predict(slm2, newdata = test5p, type = "lp", se.fit =TRUE)
start_time35 <- Sys.time()
P35<-predict(slm3, newdata = test5p, type = "surv", se.fit =TRUE)
end_time35 <- Sys.time()
pred_time35 <- end_time35 - start_time35
LP35<-predict(slm3, newdata = test5p, type = "lp", se.fit =TRUE)

test5$LP1<-as.numeric(LP15$fit)
test5$LP2<-as.numeric(LP25$fit)
test5$LP3<-as.numeric(LP35$fit)

test5$P1<-as.numeric(P15$fit)
test5$P2<-as.numeric(P25$fit)
test5$P3<-as.numeric(P35$fit)

#Calibration slope

surv15 <- coxph(Surv(timevent,event)~LP1, method="breslow", data = test5)
x15<-confint(surv15)
Cslope1[6]<-surv15$coefficients[1]
Cslope1L[6]<-x15[1]
Cslope1U[6]<-x15[2]

surv25 <- coxph(Surv(timevent,event)~LP2, method="breslow", data = test5)
x25<-confint(surv25)
Cslope2[6]<-surv25$coefficients[1]
Cslope2L[6]<-x25[1]
Cslope2U[6]<-x25[2]

surv35 <- coxph(Surv(timevent,event)~LP3, method="breslow", data = test5)
x35<-confint(surv35)
Cslope3[6]<-surv35$coefficients[1]
Cslope3L[6]<-x35[1]
Cslope3U[6]<-x35[2]

# Harrell's C-statistic

Cind1[6]<-summary(surv15)$concordance[1]
Cind1L[6]<-summary(surv15)$concordance[1]-
  (1.96*summary(surv15)$concordance[2])
Cind1U[6]<-summary(surv15)$concordance[1]+
  (1.96*summary(surv15)$concordance[2])

Cind2[6]<-summary(surv25)$concordance[1]
Cind2L[6]<-summary(surv25)$concordance[1]-
  (1.96*summary(surv25)$concordance[2])
Cind2U[6]<-summary(surv25)$concordance[1]+
  (1.96*summary(surv25)$concordance[2])

Cind3[6]<-summary(surv35)$concordance[1]
Cind3L[6]<-summary(surv35)$concordance[1]-
  (1.96*summary(surv35)$concordance[2])
Cind3U[6]<-summary(surv35)$concordance[1]+
  (1.96*summary(surv35)$concordance[2])

# Brier score

BS1[6]<-sbrier(obj = Surv(test5$timevent,test5$event), 
               test5$P1, btime = 3.5)
BS2[6]<-sbrier(obj = Surv(test5$timevent,test5$event), 
               test5$P2, btime = 3.5)
BS3[6]<-sbrier(obj = Surv(test5$timevent,test5$event), 
               test5$P3, btime = 3.5)
# Calibration plot

par(pty="s")
par(mfrow=c(1,3))
val_ests1 <- val.surv(est.surv = P15$fit, S = Surv(test5$timevent, 
                                                   test5$event), 
                      u=3.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests1,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P15$fit), S = Surv(test5$timevent, test5$event), 
           g=10,u=3.5, pl=T, add=T,lty=0, col=cols[9], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[9], "black", "black"), bty="n")

val_ests2 <- val.surv(est.surv = P25$fit, S = Surv(test5$timevent, 
                                                   test5$event), 
                      u=3.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests2,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P25$fit), S = Surv(test5$timevent, test5$event), 
           g=10,u=3.5, pl=T, add=T,lty=0, col=cols[10], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[10], "black", "black"), bty="n")

val_ests3 <- val.surv(est.surv = P35$fit, S = Surv(test5$timevent, 
                                                   test5$event), 
                      u=3.5, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests3,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P35$fit), S = Surv(test5$timevent, test5$event), 
           g=10,u=3.5, pl=T, add=T,lty=0, col=cols[11], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[11], "black", "black"), bty="n")



# LT = 3 

test6 <- read.dta("validation_36months.dta")
test6$smoke<-factor(test6$smoke)
test6$trtment<-factor(test6$firsttreat)
test6$LM<-3
test6p<-test6
test6p$timevent<-rep(4,nrow(test6p))

# Predicted probabilities and linear predictor
start_time16 <- Sys.time()
P16<-predict(slm1, newdata = test6p, type = "surv", se.fit =TRUE)
end_time16 <- Sys.time()
pred_time16 <- end_time16 - start_time16
LP16<-predict(slm1, newdata = test6p, type = "lp", se.fit =TRUE)

start_time26 <- Sys.time()
P26<-predict(slm2, newdata = test6p, type = "surv", se.fit =TRUE)
end_time26 <- Sys.time()
pred_time26 <- end_time26 - start_time26
LP26<-predict(slm2, newdata = test6p, type = "lp", se.fit =TRUE)

start_time36 <- Sys.time()
P36<-predict(slm3, newdata = test6p, type = "surv", se.fit =TRUE)
end_time36 <- Sys.time()
pred_time36 <- end_time36 - start_time36
LP36<-predict(slm3, newdata = test6p, type = "lp", se.fit =TRUE)

test6$LP1<-as.numeric(LP16$fit)
test6$LP2<-as.numeric(LP26$fit)
test6$LP3<-as.numeric(LP36$fit)

test6$P1<-as.numeric(P16$fit)
test6$P2<-as.numeric(P26$fit)
test6$P3<-as.numeric(P36$fit)

#Calibration slope

surv16 <- coxph(Surv(timevent,event)~LP1, method="breslow", data = test6)
x16<-confint(surv16)
Cslope1[7]<-surv16$coefficients[1]
Cslope1L[7]<-x16[1]
Cslope1U[7]<-x16[2]

surv26 <- coxph(Surv(timevent,event)~LP2, method="breslow", data = test6)
x26<-confint(surv26)
Cslope2[7]<-surv26$coefficients[1]
Cslope2L[7]<-x26[1]
Cslope2U[7]<-x26[2]

surv36 <- coxph(Surv(timevent,event)~LP3, method="breslow", data = test6)
x36<-confint(surv36)
Cslope3[7]<-surv36$coefficients[1]
Cslope3L[7]<-x36[1]
Cslope3U[7]<-x36[2]

# Harrell's C-statistic

Cind1[7]<-summary(surv16)$concordance[1]
Cind1L[7]<-summary(surv16)$concordance[1]-
  (1.96*summary(surv16)$concordance[2])
Cind1U[7]<-summary(surv16)$concordance[1]+
  (1.96*summary(surv16)$concordance[2])

Cind2[7]<-summary(surv26)$concordance[1]
Cind2L[7]<-summary(surv26)$concordance[1]-
  (1.96*summary(surv26)$concordance[2])
Cind2U[7]<-summary(surv26)$concordance[1]+
  (1.96*summary(surv26)$concordance[2])

Cind3[7]<-summary(surv36)$concordance[1]
Cind3L[7]<-summary(surv36)$concordance[1]-
  (1.96*summary(surv36)$concordance[2])
Cind3U[7]<-summary(surv36)$concordance[1]+
  (1.96*summary(surv36)$concordance[2])

# Brier score

BS1[7]<-sbrier(obj = Surv(test6$timevent,test6$event), 
               test6$P1, btime = 4)
BS2[7]<-sbrier(obj = Surv(test6$timevent,test6$event), 
               test6$P2, btime = 4)
BS3[7]<-sbrier(obj = Surv(test6$timevent,test6$event), 
               test6$P3, btime = 4)
# Calibration plot

par(pty="s")
par(mfrow=c(1,3))
val_ests1 <- val.surv(est.surv = P16$fit, S = Surv(test6$timevent, 
                                                   test6$event), 
                      u=4, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests1,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P16$fit), S = Surv(test6$timevent, test6$event), 
           g=10,u=4, pl=T, add=T,lty=0, col=cols[9], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[9], "black", "black"), bty="n")

val_ests2 <- val.surv(est.surv = P26$fit, S = Surv(test6$timevent, 
                                                   test6$event), 
                      u=4, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests2,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P26$fit), S = Surv(test6$timevent, test6$event), 
           g=10,u=4, pl=T, add=T,lty=0, col=cols[10], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[10], "black", "black"), bty="n")

val_ests3 <- val.surv(est.surv = P36$fit, S = Surv(test6$timevent, 
                                                   test6$event), 
                      u=4, fun=function(p)log(-log(p)), 
                      pred = sort(runif(100, 0.85, 1)))
plot(val_ests3,xlab="Expected Survival Probability",
     ylab="Observed Survival Probability", riskdist = F) 
x<-groupkm((P36$fit), S = Surv(test6$timevent, test6$event), 
           g=10,u=4, pl=T, add=T,lty=0, col=cols[11], cex.subtitle=FALSE)
legend("topleft",c("Risk groups","Reference line","95% CI"),
       lty=c(0,2,1),pch=c(19,NA,NA), col=c(cols[11], "black", "black"), bty="n")


#### Temporal assessment results
fit_times1 <- rep(fit_time1, times=7)
fit_times2 <- rep(fit_time2, times=7)
fit_times3 <- rep(fit_time3, times=7)

pred_times1 <- c(pred_time10, pred_time11, pred_time12, pred_time13, pred_time14, pred_time15, pred_time16)
pred_times2 <- c(pred_time20, pred_time21, pred_time22, pred_time23, pred_time24, pred_time25, pred_time26)
pred_times3 <- c(pred_time30, pred_time31, pred_time32, pred_time33, pred_time34, pred_time35, pred_time36)

LM1.TA<-cbind(Cslope1, Cslope1L, Cslope1U, Cind1, Cind1L, Cind1U, BS1, fit_times1, pred_times1)
LM2.TA<-cbind(Cslope2, Cslope2L, Cslope2U, Cind2, Cind2L, Cind2U, BS2, fit_times2, pred_times2)
LM3.TA<-cbind(Cslope3, Cslope3L, Cslope3U, Cind3, Cind3L, Cind3U, BS3, fit_times3, pred_times3)

LM.TA<-rbind(LM1.TA, LM2.TA, LM3.TA)
LM.TA<-as.data.frame(LM.TA)
write.table(LM.TA, file = "lmTA.csv", sep = ",", col.names = NA,
            qmethod = "double")
x<-read.table("lmTA.csv", header = TRUE, sep = ",", row.names = 1)
head(x)


########################## Bootstrap function #################################

bootvadLM<-function(SLMdata, B, seed){
  dat <- SLMdata
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
  
  # original data - temporal assessment
  
  od0 <- test0
  od0p <- test0p
  od1 <- test1
  od1p <- test1p
  od2 <- test2
  od2p <- test2p
  od3 <- test3
  od3p <- test3p
  od4 <- test4
  od4p <- test4p
  od5 <- test5
  od5p <- test5p
  od6 <- test6
  od6p <- test6p
  
  od0$studyno<-od0$groupid
  od1$studyno<-od1$groupid
  od2$studyno<-od2$groupid
  od3$studyno<-od3$groupid
  od4$studyno<-od4$groupid
  od5$studyno<-od5$groupid
  od6$studyno<-od6$groupid
  
  od0p$studyno<-od0p$groupid
  od1p$studyno<-od1p$groupid
  od2p$studyno<-od2p$groupid
  od3p$studyno<-od3p$groupid
  od4p$studyno<-od4p$groupid
  od5p$studyno<-od5p$groupid
  od6p$studyno<-od6p$groupid
  
  for (j in 1:B){
    seed<-seed+j
    set.seed(seed)
    smp <- sort(sample(indiv, length(indiv), replace=TRUE))
    smp.df <- data.frame(studyno=smp)
    print(length(smp.df$studyno))

    # renaming duplicated ID numbers
    
    smp.df<-smp.df %>% group_by(studyno) %>% mutate(count = row_number())
    smp.df$count<-seq(1, nrow(smp.df), by = 1)
    smp.df$newid<-paste(smp.df$studyno, smp.df$count, sep=".")
    dat.b = merge(smp.df, dat, by = "studyno", all.x=TRUE)
    dat.b$studyno <- dat.b$newid
    dat.b<-dat.b[with(dat.b, order(studyno, LM)),]
    print(length(dat.b$studyno))
    
    # fitting LM
    
    slm <- coxph(Surv(LM,timevent,event) ~ age + pgen + previous_dmards + 
                    disdur + as.factor(trtment) + ovmean + dascore + steroids +
                    bmi + as.factor(smoke) + renal + lung + diabetes + strata(LM) + 
                    cluster(studyno), 
                  data=dat.b, method="breslow")
    
    # Temporal assessment sets
    
    val0 <- na.omit(merge(smp.df, od0, by = "studyno", all.x=TRUE))
    val0p <- na.omit(merge(smp.df, od0p, by = "studyno", all.x=TRUE))
    val1 <- na.omit(merge(smp.df, od1, by = "studyno", all.x=TRUE))
    val1p <- na.omit(merge(smp.df, od1p, by = "studyno", all.x=TRUE))
    val2 <- na.omit(merge(smp.df, od2, by = "studyno", all.x=TRUE))
    val2p <- na.omit(merge(smp.df, od2p, by = "studyno", all.x=TRUE))
    val3 <- na.omit(merge(smp.df, od3, by = "studyno", all.x=TRUE))
    val3p <- na.omit(merge(smp.df, od3p, by = "studyno", all.x=TRUE))
    val4 <- na.omit(merge(smp.df, od4, by = "studyno", all.x=TRUE))
    val4p <- na.omit(merge(smp.df, od4p, by = "studyno", all.x=TRUE))
    val5 <- na.omit(merge(smp.df, od5, by = "studyno", all.x=TRUE))
    val5p <- na.omit(merge(smp.df, od5p, by = "studyno", all.x=TRUE))
    val6 <- na.omit(merge(smp.df, od6, by = "studyno", all.x=TRUE))
    val6p <- na.omit(merge(smp.df, od6p, by = "studyno", all.x=TRUE))
    
    # Survival predictions at different time points
    
    P0<-predict(slm, newdata = val0p, type = "surv", se.fit =TRUE)
    LP0<-predict(slm, newdata = val0p, type = "lp", se.fit =TRUE)
    P1<-predict(slm, newdata = val1p, type = "surv", se.fit =TRUE)
    LP1<-predict(slm, newdata = val1p, type = "lp", se.fit =TRUE)
    P2<-predict(slm, newdata = val2p, type = "surv", se.fit =TRUE)
    LP2<-predict(slm, newdata = val2p, type = "lp", se.fit =TRUE)
    P3<-predict(slm, newdata = val3p, type = "surv", se.fit =TRUE)
    LP3<-predict(slm, newdata = val3p, type = "lp", se.fit =TRUE)
    P4<-predict(slm, newdata = val4p, type = "surv", se.fit =TRUE)
    LP4<-predict(slm, newdata = val4p, type = "lp", se.fit =TRUE)
    P5<-predict(slm, newdata = val5p, type = "surv", se.fit =TRUE)
    LP5<-predict(slm, newdata = val5p, type = "lp", se.fit =TRUE)
    P6<-predict(slm, newdata = val6p, type = "surv", se.fit =TRUE)
    LP6<-predict(slm, newdata = val6p, type = "lp", se.fit =TRUE)
    
    val0$LP0<-as.numeric(LP0$fit)
    val1$LP1<-as.numeric(LP1$fit)
    val2$LP2<-as.numeric(LP2$fit)
    val3$LP3<-as.numeric(LP3$fit)
    val4$LP4<-as.numeric(LP4$fit)
    val5$LP5<-as.numeric(LP5$fit)
    val6$LP6<-as.numeric(LP6$fit)
    
    val0$P0<-as.numeric(P0$fit)
    val1$P1<-as.numeric(P1$fit)
    val2$P2<-as.numeric(P2$fit)
    val3$P3<-as.numeric(P3$fit)
    val4$P4<-as.numeric(P4$fit)
    val5$P5<-as.numeric(P5$fit)
    val6$P6<-as.numeric(P6$fit)
    
    # Discrimination measure
    
    C[1]<- rcorr.cens(val0$P0, 
                      Surv(val0$timevent, val0$event))["C Index"]
    C[2]<- rcorr.cens(val1$P1, 
                      Surv(val1$timevent, val1$event))["C Index"]
    C[3]<- rcorr.cens(val2$P2, 
                      Surv(val2$timevent, val2$event))["C Index"]
    C[4]<- rcorr.cens(val3$P3, 
                      Surv(val3$timevent, val3$event))["C Index"]
    C[5]<- rcorr.cens(val4$P4, 
                      Surv(val4$timevent, val4$event))["C Index"]
    C[6]<- rcorr.cens(val5$P5, 
                      Surv(val5$timevent, val5$event))["C Index"]
    C[7]<- rcorr.cens(val6$P6, 
                      Surv(val6$timevent, val6$event))["C Index"]
    
    c[j,]<-C
    
    # Calibration measure
    
    bs0<-sbrier(obj = Surv(val0$timevent,val0$event), 
                val0$P0, btime = 1)
    bs1<-sbrier(obj = Surv(val1$timevent,val1$event), 
                val1$P1, btime = 1.5)
    bs2<-sbrier(obj = Surv(val2$timevent,val2$event), 
                val2$P2, btime = 2)
    bs3<-sbrier(obj = Surv(val3$timevent,val3$event), 
                val3$P3, btime = 2.5)
    bs4<-sbrier(obj = Surv(val4$timevent,val4$event), 
                val4$P4, btime = 3)
    bs5<-sbrier(obj = Surv(val5$timevent,val5$event), 
                val5$P5, btime = 3.5)
    bs6<-sbrier(obj = Surv(val6$timevent,val6$event), 
                val6$P6, btime = 4)
    
    bs[j,]<-c(bs0, bs1, bs2, bs3, bs4, bs5, bs6)
    
    # Survival predictions - original data
    
    oP0<-predict(slm, newdata = od0p, type = "surv", se.fit =TRUE)
    oLP0<-predict(slm, newdata = od0p, type = "lp", se.fit =TRUE)
    oP1<-predict(slm, newdata = od1p, type = "surv", se.fit =TRUE)
    oLP1<-predict(slm, newdata = od1p, type = "lp", se.fit =TRUE)
    oP2<-predict(slm, newdata = od2p, type = "surv", se.fit =TRUE)
    oLP2<-predict(slm, newdata = od2p, type = "lp", se.fit =TRUE)
    oP3<-predict(slm, newdata = od3p, type = "surv", se.fit =TRUE)
    oLP3<-predict(slm, newdata = od3p, type = "lp", se.fit =TRUE)
    oP4<-predict(slm, newdata = od4p, type = "surv", se.fit =TRUE)
    oLP4<-predict(slm, newdata = od4p, type = "lp", se.fit =TRUE)
    oP5<-predict(slm, newdata = od5p, type = "surv", se.fit =TRUE)
    oLP5<-predict(slm, newdata = od5p, type = "lp", se.fit =TRUE)
    oP6<-predict(slm, newdata = od6p, type = "surv", se.fit =TRUE)
    oLP6<-predict(slm, newdata = od6p, type = "lp", se.fit =TRUE)
    
    od0$oLP0<-as.numeric(oLP0$fit)
    od1$oLP1<-as.numeric(oLP1$fit)
    od2$oLP2<-as.numeric(oLP2$fit)
    od3$oLP3<-as.numeric(oLP3$fit)
    od4$oLP4<-as.numeric(oLP4$fit)
    od5$oLP5<-as.numeric(oLP5$fit)
    od6$oLP6<-as.numeric(oLP6$fit)
    
    od0$oP0<-as.numeric(oP0$fit)
    od1$oP1<-as.numeric(oP1$fit)
    od2$oP2<-as.numeric(oP2$fit)
    od3$oP3<-as.numeric(oP3$fit)
    od4$oP4<-as.numeric(oP4$fit)
    od5$oP5<-as.numeric(oP5$fit)
    od6$oP6<-as.numeric(oP6$fit)
    
    # Discrimination measure
    
    od.C[1]<- rcorr.cens(od0$oP0, 
                         Surv(od0$timevent, od0$event))["C Index"]
    od.C[2]<- rcorr.cens(od1$oP1, 
                         Surv(od1$timevent, od1$event))["C Index"]
    od.C[3]<- rcorr.cens(od2$oP2, 
                         Surv(od2$timevent, od2$event))["C Index"]
    od.C[4]<- rcorr.cens(od3$oP3, 
                         Surv(od3$timevent, od3$event))["C Index"]
    od.C[5]<- rcorr.cens(od4$oP4, 
                         Surv(od4$timevent, od4$event))["C Index"]
    od.C[6]<- rcorr.cens(od5$oP5, 
                         Surv(od5$timevent, od5$event))["C Index"]
    od.C[7]<- rcorr.cens(od6$oP6, 
                         Surv(od6$timevent, od6$event))["C Index"]
    
    od.c[j,]<-od.C

    # Calibration measure
    
    bs0.od<-sbrier(obj = Surv(od0$timevent,od0$event), 
                   od0$oP0, btime = 1)
    bs1.od<-sbrier(obj = Surv(od1$timevent,od1$event), 
                   od1$oP1, btime = 1.5)
    bs2.od<-sbrier(obj = Surv(od2$timevent,od2$event), 
                   od2$oP2, btime = 2)
    bs3.od<-sbrier(obj = Surv(od3$timevent,od3$event), 
                   od3$oP3, btime = 2.5)
    bs4.od<-sbrier(obj = Surv(od4$timevent,od4$event), 
                   od4$oP4, btime = 3)
    bs5.od<-sbrier(obj = Surv(od5$timevent,od5$event), 
                   od5$oP5, btime = 3.5)
    bs6.od<-sbrier(obj = Surv(od6$timevent,od6$event), 
                   od6$oP6, btime = 4)
    
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

#bootvadLM(SLMdata.3, 1, 767)

########################## Parallelisation ####################################

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1) #always useful to keep one core free
start_time<-Sys.time()
registerDoParallel(cl)
LM.opt <- foreach(i=1:200, .combine=rbind, .packages = c("ipred", "dplyr", 
                                                         "dynpred", "survival", "Hmisc")) %dopar% {
                                                           seed<-199429+i
                                                           bootvadLM(SLMdata.3, 1, seed)
                                                         }
stopCluster(cl)
end_time<-Sys.time()
end_time - start_time

LM.opt.C=LM.opt[seq(1,399,2),]
LM.opt.B=LM.opt[seq(2,400,2),]
LM.opt2<-matrix(0, nrow=2, ncol=7)
LM.opt2[1,]=c(mean(LM.opt.C[,1]), mean(LM.opt.C[,2]),
               mean(LM.opt.C[,3]), mean(LM.opt.C[,4]),
               mean(LM.opt.C[,5]), mean(LM.opt.C[,6]),
               mean(LM.opt.C[,7]))
LM.opt2[2,]=c(mean(LM.opt.B[,1]), mean(LM.opt.B[,2]),
               mean(LM.opt.B[,3]), mean(LM.opt.B[,4]),
               mean(LM.opt.B[,5]), mean(LM.opt.B[,6]),
               mean(LM.opt.B[,7]))

LM.OPT<-data.frame(LT0 = LM.opt2[,1], LT0.5 = LM.opt2[,2], LT1 = LM.opt2[,3], 
                    LT1.5 = LM.opt2[,4], LT2 = LM.opt2[,5], LT2.5 = LM.opt2[,6], 
                    LT3 = LM.opt2[,7], row.names = c("C statistic (OE)", 
                                                      "Brier Score (OE)"))

write.table(LM.OPT, file = "lmOPT.csv", sep = ",", col.names = NA,
            qmethod = "double")
