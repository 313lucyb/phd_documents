library(foreign)
library(mfp)
library(caret)
library(rms)
library(ROCR)
library(survival)
library(Hmisc)
library(MASS)
library(pROC)
library(DescTools)
library(plyr)
library(boot)
library(party)
library(risksetROC)
library(colorRamps)
library(grDevices)
library(dplyr)
library(RColorBrewer)
cols <-brewer.pal(11, "RdBu")

# CPM Development

LR<-read.dta("baseline_binary.dta")
LR$smoke<-factor(LR$smoke)
LR$firsttreat<-factor(LR$firsttreat)

attach(LR)

# Relationships between continuous predictors

# Scatterplot matrix
par(pty="m")
pairs(cbind(ovmean,age, bmi, dascore, previous_dmards, disdur),
      col="indianred",labels=c("HAQ","Age","BMI","DAS28", "No.\n prev. DMARDs", 
                         "Disease\n duration"))

# Correlation coefficients
round(rcorr(cbind(ovmean,age, bmi, dascore, previous_dmards, 
                  disdur),type="pearson")$r,6)

# Checking nonlinearities

# Univariable analysis - no nonlinear relationships

set.seed(12)

LR_fp1<-mfp(event~fp(age), family = binomial, data = LR)
LR_fp2<-mfp(event~fp(dascore), family = binomial, data = LR)
LR_fp3<-mfp(event~fp(ovmean), family = binomial, data = LR)
LR_fp5<-mfp(event~fp(bmi), family = binomial, data = LR)
LR_fp6<-mfp(event~fp(disdur), family = binomial, data = LR)
LR_fp7<-mfp(event~fp(previous_dmards), family = binomial, data = LR)

# multivariable analysis

LR_fp<-mfp(event~fp(age)+pgen+fp(dascore)+fp(ovmean)+
             steroids+firsttreat+fp(bmi)+smoke+lung+
             renal+diabetes+fp(disdur)+fp(previous_dmards), 
           family = binomial, data = LR)


# Final model

set.seed(123)
start_time <- Sys.time()
m1<-glm(event~age + pgen + firsttreat + 
          disdur + previous_dmards + lung + 
          diabetes + bmi + renal + steroids + ovmean + 
          smoke + dascore, data=LR, family=binomial)
end_time <- Sys.time()
fit_time <- end_time - start_time


###########################################################
##                  Temporal Assessment                  ##
###########################################################

## PA measure vectors

LT<-c(0,0.5,1,1.5,2, 2.5, 3)
Clarge<-0
ClargeL<-0
ClargeU<-0
Cslope<-0
CslopeL<-0
CslopeU<-0
Cstat<-0
CstatL<-0
CstatU<-0
BS<-0


##### LT = 0 #####


# Obtain the predicted probabilities for each subject
start_time0<-Sys.time()
pred_prob0 <- predict(m1,type="response")
end_time0<-Sys.time()
pred_time0 <- end_time0 - start_time0
# Obtain the linear predictor/PI for each subject
pred_LP0 <- predict(m1,type="link")

# # Discrimination

# Obtain the c statistic / AUC
c0 <- roc(LR$event~pred_prob0,ci=TRUE)
Cstat[1]<-as.numeric(c0$auc[1])
CstatL[1]<-as.numeric(c0$ci[1])
CstatU[1]<-as.numeric(c0$ci[3])

# # Calibration

# Calibration and other point estimates

# Brier score
BS[1]<-as.numeric(BrierScore(m1))

# Calibration In The Large 
set.seed(1234)
m1_0 <- glm(LR$event~offset(pred_LP0),family="binomial")
Clarge[1]<-as.numeric(m1_0$coef)
x0<-confint(m1_0)
ClargeL[1]<-as.numeric(x0[1])
ClargeU[1]<-as.numeric(x0[2])

# C-slope
m1_00 <- glm(LR$event~pred_LP0,family="binomial",x=TRUE,y=TRUE)
Cslope[1]<-m1_00$coefficients[2]
x00<-confint(m1_00)
CslopeL[1]<-as.numeric(x00[2,1])
CslopeU[1]<-as.numeric(x00[2,2])

# Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(pred_prob0,breaks=quantile(pred_prob0, 
                                         prob = c(0,0.1,0.2,0.3,
                                                  0.4,0.5,0.6,0.7,
                                                  0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(LR,groups,pred_prob0)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(pred_prob0))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

par(pty="s")
plot(obs~exp[,2],xlim=c(0,0.2),ylim=c(0,0.2),col=cols[1], 
     pch =16, ylab="Observed",xlab="Expected", cex.lab = 1.2)
lines(c(0,0.2),c(0,0.2),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1])
}
h <- hist(pred_prob0, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(event))~pred_prob0,span=1))
lines_data <- data.frame(pred_prob0,obs_all)
lines_data2 <- lines_data[order(pred_prob0),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend(0.0,0.2,c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")



##### LT = 0.5 #####



# Data for predictions/validation

sixmth<-read.dta("validation_6binary.dta")
sixmth$smoke<-factor(sixmth$smoke)
sixmth$firsttreat<-factor(sixmth$firsttreat)



# Obtain the predicted probabilities for each patient
start_time1 <- Sys.time()
pred_prob1 <- predict(m1, newdata = sixmth, type="response")
end_time1<-Sys.time()
pred_time1 <- end_time1 - start_time1
# Obtain the linear predictor/PI for each patient
pred_LP1 <- predict(m1, newdata = sixmth, type="link")


# # Discrimination

# Obtain the c statistic / AUC
c1 <- roc(sixmth$event~pred_prob1,ci=TRUE)
Cstat[2]<-as.numeric(c1$auc[1])
CstatL[2]<-as.numeric(c1$ci[1])
CstatU[2]<-as.numeric(c1$ci[3])

# # Calibration

# Calibration and other point estimates
BS[2]<-BrierScore(resp= sixmth$event, pred = pred_prob1, 
                  scaled = FALSE)

# Calibration In The Large 
set.seed(12345)
m1_1 <- glm(sixmth$event~offset(pred_LP1),family="binomial")
Clarge[2]<-as.numeric(m1_1$coef)
x1<-confint(m1_1)
ClargeL[2]<-as.numeric(x1[1])
ClargeU[2]<-as.numeric(x1[2])

# C-slope
m1_11 <- glm(sixmth$event~pred_LP1,family="binomial",x=TRUE,y=TRUE)
Cslope[2]<-m1_11$coefficients[2]
x11<-confint(m1_11)
CslopeL[2]<-as.numeric(x11[2,1])
CslopeU[2]<-as.numeric(x11[2,2])

# Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(pred_prob1,breaks=quantile(pred_prob1, 
                                         prob = c(0,0.1,0.2,0.3,
                                                  0.4,0.5,0.6,0.7,
                                                  0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(sixmth,groups,pred_prob1)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(pred_prob1))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

par(pty="s")
plot(obs~exp[,2],xlim=c(0,0.2),ylim=c(0,0.2),col=cols[1], 
     pch =16, ylab="Observed",xlab="Expected", cex.lab = 1.2)
lines(c(0,0.2),c(0,0.2),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1])
}
h <- hist(pred_prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(event))~pred_prob1,span=1))
lines_data <- data.frame(pred_prob1,obs_all)
lines_data2 <- lines_data[order(pred_prob1),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend(0.0,0.2,c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


##### LT = 1 #####



# Data for predictions/validation

twmth<-read.dta("validation_12binary.dta")
twmth$smoke<-factor(twmth$smoke)
twmth$firsttreat<-factor(twmth$firsttreat)


# Obtain the predicted probabilities for each patient
start_time2 <- Sys.time()
pred_prob2 <- predict(m1, newdata = twmth, type="response")
end_time2<-Sys.time()
pred_time2 <- end_time2 - start_time2

# Obtain the linear predictor/PI for each patient
pred_LP2 <- predict(m1, newdata = twmth, type="link")


# # Discrimination

# Obtain the c statistic / AUC
c2 <- roc(twmth$event~pred_prob2,ci=TRUE)
Cstat[3]<-as.numeric(c2$auc[1])
CstatL[3]<-as.numeric(c2$ci[1])
CstatU[3]<-as.numeric(c2$ci[3])

# # Calibration

# Calibration and other point estimates
BS[3]<-BrierScore(resp= twmth$event, pred = pred_prob2, 
                  scaled = FALSE)

# Calibration In The Large 
set.seed(123456)
m1_2 <- glm(twmth$event~offset(pred_LP2),family="binomial")
Clarge[3]<-as.numeric(m1_2$coef)
x2<-confint(m1_2)
ClargeL[3]<-as.numeric(x2[1])
ClargeU[3]<-as.numeric(x2[2])

# C-slope
m1_22 <- glm(twmth$event~pred_LP2,family="binomial",x=TRUE,y=TRUE)
Cslope[3]<-m1_22$coefficients[2]
x22<-confint(m1_22)
CslopeL[3]<-as.numeric(x22[2,1])
CslopeU[3]<-as.numeric(x22[2,2])

# Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(pred_prob2,breaks=quantile(pred_prob2, 
                                         prob = c(0,0.1,0.2,0.3,
                                                  0.4,0.5,0.6,0.7,
                                                  0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(twmth,groups,pred_prob2)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(pred_prob2))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

par(pty="s")
plot(obs~exp[,2],xlim=c(0,0.2),ylim=c(0,0.2),col=cols[1], 
     pch =16, ylab="Observed",xlab="Expected", cex.lab = 1.2)
lines(c(0,0.2),c(0,0.2),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1])
}
h <- hist(pred_prob2, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(event))~pred_prob2,span=1))
lines_data <- data.frame(pred_prob2,obs_all)
lines_data2 <- lines_data[order(pred_prob2),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend(0.0,0.2,c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")



##### LT = 1.5 #####



# Data for predictions/validation

emth<-read.dta("validation_18binary.dta")
emth$smoke<-factor(emth$smoke)
emth$firsttreat<-factor(emth$firsttreat)



# Obtain the predicted probabilities for each patient
start_time3 <- Sys.time()
pred_prob3 <- predict(m1, newdata = emth, type="response")
end_time3 <- Sys.time()
pred_time3 <- end_time3 - start_time3

# Obtain the linear predictor/PI for each patient
pred_LP3 <- predict(m1, newdata = emth, type="link")


# # Discrimination

# Obtain the c statistic / AUC
c3 <- roc(emth$event~pred_prob3,ci=TRUE)
Cstat[4]<-as.numeric(c3$auc[1])
CstatL[4]<-as.numeric(c3$ci[1])
CstatU[4]<-as.numeric(c3$ci[3])

# # Calibration

# Calibration and other point estimates
BS[4]<-BrierScore(resp= emth$event, pred = pred_prob3, 
                  scaled = FALSE)

# Calibration In The Large 
set.seed(1234567)
m1_3 <- glm(emth$event~offset(pred_LP3),family="binomial")
Clarge[4]<-as.numeric(m1_3$coef)
x3<-confint(m1_3)
ClargeL[4]<-as.numeric(x3[1])
ClargeU[4]<-as.numeric(x3[2])

# C-slope
m1_33 <- glm(emth$event~pred_LP3,family="binomial",x=TRUE,y=TRUE)
Cslope[4]<-m1_33$coefficients[2]
x33<-confint(m1_33)
CslopeL[4]<-as.numeric(x33[2,1])
CslopeU[4]<-as.numeric(x33[2,2])

# Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(pred_prob3,breaks=quantile(pred_prob3, 
                                         prob = c(0,0.1,0.2,0.3,
                                                  0.4,0.5,0.6,0.7,
                                                  0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(emth,groups,pred_prob3)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(pred_prob3))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

par(pty="s")
plot(obs~exp[,2],xlim=c(0,0.2),ylim=c(0,0.2),col=cols[1], 
     pch =16, ylab="Observed",xlab="Expected", cex.lab = 1.2)
lines(c(0,0.2),c(0,0.2),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1])
}
h <- hist(pred_prob3, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(event))~pred_prob3,span=1))
lines_data <- data.frame(pred_prob3,obs_all)
lines_data2 <- lines_data[order(pred_prob3),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend(0.0,0.2,c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")



##### LT = 2 #####


# Data for predictions/validation

tfmth<-read.dta("validation_24binary.dta")
tfmth$smoke<-factor(tfmth$smoke)
tfmth$firsttreat<-factor(tfmth$firsttreat)



# Obtain the predicted probabilities for each patient
start_time4 <- Sys.time()
pred_prob4 <- predict(m1, newdata = tfmth, type="response")
end_time4 <- Sys.time()
pred_time4 <- end_time4 - start_time4

# Obtain the linear predictor/PI for each patient
pred_LP4 <- predict(m1, newdata = tfmth, type="link")


# # Discrimination

# Obtain the c statistic / AUC
c4 <- roc(tfmth$event~pred_prob4,ci=TRUE)
Cstat[5]<-as.numeric(c4$auc[1])
CstatL[5]<-as.numeric(c4$ci[1])
CstatU[5]<-as.numeric(c4$ci[3])

# # Calibration

# Calibration and other point estimates
BS[5]<-BrierScore(resp= tfmth$event, pred = pred_prob4, 
                  scaled = FALSE)

# Calibration In The Large 
set.seed(12345678)
m1_4 <- glm(tfmth$event~offset(pred_LP4),family="binomial")
Clarge[5]<-as.numeric(m1_4$coef)
x4<-confint(m1_4)
ClargeL[5]<-as.numeric(x4[1])
ClargeU[5]<-as.numeric(x4[2])

# C-slope
m1_44 <- glm(tfmth$event~pred_LP4,family="binomial",x=TRUE,y=TRUE)
Cslope[5]<-m1_44$coefficients[2]
x44<-confint(m1_44)
CslopeL[5]<-as.numeric(x44[2,1])
CslopeU[5]<-as.numeric(x44[2,2])

# Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(pred_prob4,breaks=quantile(pred_prob4, 
                                         prob = c(0,0.1,0.2,0.3,
                                                  0.4,0.5,0.6,0.7,
                                                  0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(tfmth,groups,pred_prob4)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(pred_prob4))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

par(pty="s")
plot(obs~exp[,2],xlim=c(0,0.2),ylim=c(0,0.2),col=cols[1], 
     pch =16, ylab="Observed",xlab="Expected", cex.lab=1.2)
lines(c(0,0.2),c(0,0.2),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1])
}
h <- hist(pred_prob4, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(event))~pred_prob4,span=1))
lines_data <- data.frame(pred_prob4,obs_all)
lines_data2 <- lines_data[order(pred_prob4),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend(0.0,0.2,c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


##### LT = 2.5 #####


# Data for predictions/validation

thmth<-read.dta("validation_30binary.dta")
thmth$smoke<-factor(thmth$smoke)
thmth$firsttreat<-factor(thmth$firsttreat)


# Obtain the predicted probabilities for each patient
start_time5 <- Sys.time()
pred_prob5 <- predict(m1, newdata = thmth, type="response")
end_time5<-Sys.time()
pred_time5 <- end_time5 - start_time5
# Obtain the linear predictor/PI for each patient
pred_LP5 <- predict(m1, newdata = thmth, type="link")


# # Discrimination

# Obtain the c statistic / AUC
c5 <- roc(thmth$event~pred_prob5,ci=TRUE)
Cstat[6]<-as.numeric(c5$auc[1])
CstatL[6]<-as.numeric(c5$ci[1])
CstatU[6]<-as.numeric(c5$ci[3])

# # Calibration

# Calibration and other point estimates
BS[6]<-BrierScore(resp= thmth$event, pred = pred_prob5, 
                  scaled = FALSE)

# Calibration In The Large 
set.seed(123456789)
m1_5 <- glm(thmth$event~offset(pred_LP5),family="binomial")
Clarge[6]<-as.numeric(m1_5$coef)
x5<-confint(m1_5)
ClargeL[6]<-as.numeric(x5[1])
ClargeU[6]<-as.numeric(x5[2])

# C-slope
m1_55 <- glm(thmth$event~pred_LP5,family="binomial",x=TRUE,y=TRUE)
Cslope[6]<-m1_55$coefficients[2]
x55<-confint(m1_55)
CslopeL[6]<-as.numeric(x55[2,1])
CslopeU[6]<-as.numeric(x55[2,2])

# Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(pred_prob5,breaks=quantile(pred_prob5, 
                                         prob = c(0,0.1,0.2,0.3,
                                                  0.4,0.5,0.6,0.7,
                                                  0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(thmth,groups,pred_prob5)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(pred_prob5))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

par(pty="s")
plot(obs~exp[,2],xlim=c(0,0.2),ylim=c(0,0.2),col=cols[1], 
     pch =16, ylab="Observed",xlab="Expected", cex.lab=1.2)
lines(c(0,0.2),c(0,0.2),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1])
}
h <- hist(pred_prob5, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(event))~pred_prob5,span=1))
lines_data <- data.frame(pred_prob5,obs_all)
lines_data2 <- lines_data[order(pred_prob5),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend(0.0,0.2,c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


##### LT = 3 #####


# Data for predictions/validation

tsmth<-read.dta("validation_36binary.dta")
tsmth$smoke<-factor(tsmth$smoke)
tsmth$firsttreat<-factor(tsmth$firsttreat)



# Obtain the predicted probabilities for each patient
start_time6 <- Sys.time()
pred_prob6 <- predict(m1, newdata = tsmth, type="response")
end_time6<-Sys.time()
pred_time6 <- end_time6 - start_time6
# Obtain the linear predictor/PI for each patient
pred_LP6 <- predict(m1, newdata = tsmth, type="link")


# # Discrimination

# Obtain the c statistic / AUC
c6 <- roc(tsmth$event~pred_prob6,ci=TRUE)
Cstat[7]<-as.numeric(c6$auc[1])
CstatL[7]<-as.numeric(c6$ci[1])
CstatU[7]<-as.numeric(c6$ci[3])

# # Calibration

# Calibration and other point estimates
BS[7]<-BrierScore(resp= tsmth$event, pred = pred_prob6, 
                  scaled = FALSE)

# Calibration In The Large 
set.seed(1234567896)
m1_6 <- glm(tsmth$event~offset(pred_LP6),family="binomial")
Clarge[7]<-as.numeric(m1_6$coef)
x6<-confint(m1_6)
ClargeL[7]<-as.numeric(x6[1])
ClargeU[7]<-as.numeric(x6[2])

# C-slope
m1_66 <- glm(tsmth$event~pred_LP6,family="binomial",x=TRUE,y=TRUE)
Cslope[7]<-m1_66$coefficients[2]
x66<-confint(m1_66)
CslopeL[7]<-as.numeric(x66[2,1])
CslopeU[7]<-as.numeric(x66[2,2])

# Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(pred_prob6,breaks=quantile(pred_prob6, 
                                         prob = c(0,0.1,0.2,0.3,
                                                  0.4,0.5,0.6,0.7,
                                                  0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(tsmth,groups,pred_prob6)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(pred_prob6))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

par(pty="s")
plot(obs~exp[,2],xlim=c(0,0.2),ylim=c(0,0.2),col=cols[1], 
     pch =16, ylab="Observed",xlab="Expected", cex.lab = 1.2)
lines(c(0,0.2),c(0,0.2),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1])
}
h <- hist(pred_prob6, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(event))~pred_prob6,span=1))
lines_data <- data.frame(pred_prob6,obs_all)
lines_data2 <- lines_data[order(pred_prob6),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend(0.0,0.2,c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# Performance data frame

pred_time <- c(pred_time0, pred_time1, pred_time2, pred_time3, pred_time4, pred_time5, pred_time6)
fit_vec <- rep(fit_time, times = 7)
LRM.TA<-cbind(LT, Cslope, CslopeL, CslopeU, Clarge, ClargeL, ClargeU, Cstat, CstatL, CstatU, BS, fit_vec, pred_time)
LRM.TA<-as.data.frame(LRM.TA)
write.table(LRM.TA, file = "lrmTA.csv", sep = ",", col.names = NA,
            qmethod = "double")
x<-read.table("lrmTA.csv", header = TRUE, sep = ",", row.names = 1)
head(x)

#Bootstrap function

LRM.opt<-bootvadLRM(LR, event, 200, 199430)
LRM.opt<-as.data.frame(LRM.opt)
write.table(LRM.opt, file = "lrmOPT.csv", sep = ",", col.names = NA,
            qmethod = "double")
x<-read.table("lrmOPT.csv", header = TRUE, sep = ",", row.names = 1)
head(x)

###########################################################
##                  Temporal Assessment (survival)                  ##
###########################################################

## PA measure vectors

LT<-c(0,0.5,1,1.5,2, 2.5, 3)
Cstat<-0
CstatL<-0
CstatU<-0
BS<-0


##### LT = 0 #####


# Obtain the predicted probabilities for each subject
test0<-read.dta("baseline_model.dta")
test0$smoke <- as.factor(test0$smoke)
test0$firsttreat <- as.factor(test0$firsttreat)
test0$pred_prob0 <- predict(m1,newdata=test0,type="response")
test0$pred_LP0 <- predict(m1,newdata=test0, type="link")
test0$surv_pred0 <- 1-test0$pred_prob0

# # Discrimination

# C-index
c0 <- coxph(Surv(timevent,event)~pred_LP0, data=test0)
Cstat[1]<-as.numeric(summary(c0)$concordance[1])
CstatL[1]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[1]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

# Brier
BS[1]<-sbrier(obj = Surv(test0$timevent,test0$event), 
              test0$surv_pred0, btime = 1)

##### LT = 0.5 #####


# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_6months.dta")
test0$smoke<-factor(test0$smoke, levels=c("0","1","2"), 
                    ordered=F)
test0$firsttreat<-factor(test0$firsttreat, ordered = F)

pred_prob0 <- predict(m1,newdata=test0,type="response")
test0$pred_LP0 <- predict(m1,newdata=test0, type="link")
surv_pred0 <- 1-pred_prob0

# # Discrimination

# C-index
c0 <- coxph(Surv(timevent,event)~pred_LP0, data=test0)
Cstat[2]<-as.numeric(summary(c0)$concordance[1])
CstatL[2]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[2]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

# Brier
BS[2]<-sbrier(obj = Surv(test0$timevent,test0$event), 
              surv_pred0, btime = 1.5)


##### LT = 1 #####

# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_12months.dta")
test0$smoke<-factor(test0$smoke, levels=c("0","1","2"), 
                    ordered=F)
test0$firsttreat<-factor(test0$firsttreat, ordered = F)

pred_prob0 <- predict(m1,newdata=test0,type="response")
test0$pred_LP0 <- predict(m1,newdata=test0, type="link")
surv_pred0 <- 1-pred_prob0

# # Discrimination

# C-index
c0 <- coxph(Surv(timevent,event)~pred_LP0, data=test0)
Cstat[3]<-as.numeric(summary(c0)$concordance[1])
CstatL[3]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[3]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

# Brier
BS[3]<-sbrier(obj = Surv(test0$timevent,test0$event), 
              surv_pred0, btime = 2)



##### LT = 1.5 #####

# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_18months.dta")
test0$smoke<-factor(test0$smoke, levels=c("0","1","2"), 
                    ordered=F)
test0$firsttreat<-factor(test0$firsttreat, ordered = F)

pred_prob0 <- predict(m1,newdata=test0,type="response")
test0$pred_LP0 <- predict(m1,newdata=test0, type="link")
surv_pred0 <- 1-pred_prob0

# # Discrimination

# C-index
c0 <- coxph(Surv(timevent,event)~pred_LP0, data=test0)
Cstat[4]<-as.numeric(summary(c0)$concordance[1])
CstatL[4]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[4]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

# Brier
BS[4]<-sbrier(obj = Surv(test0$timevent,test0$event), 
              surv_pred0, btime = 2.5)

##### LT = 2 #####


# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_24months.dta")
test0$smoke<-factor(test0$smoke, levels=c("0","1","2"), 
                    ordered=F)
test0$firsttreat<-factor(test0$firsttreat, ordered = F)

pred_prob0 <- predict(m1,newdata=test0,type="response")
test0$pred_LP0 <- predict(m1,newdata=test0, type="link")
surv_pred0 <- 1-pred_prob0

# # Discrimination

# C-index
c0 <- coxph(Surv(timevent,event)~pred_LP0, data=test0)
Cstat[5]<-as.numeric(summary(c0)$concordance[1])
CstatL[5]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[5]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

# Brier
BS[5]<-sbrier(obj = Surv(test0$timevent,test0$event), 
              surv_pred0, btime = 3)


##### LT = 2.5 #####


# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_30months.dta")
test0$smoke<-factor(test0$smoke, levels=c("0","1","2"), 
                    ordered=F)
test0$firsttreat<-factor(test0$firsttreat, ordered = F)

pred_prob0 <- predict(m1,newdata=test0,type="response")
test0$pred_LP0 <- predict(m1,newdata=test0, type="link")
surv_pred0 <- 1-pred_prob0

# # Discrimination

# C-index
c0 <- coxph(Surv(timevent,event)~pred_LP0, data=test0)
Cstat[6]<-as.numeric(summary(c0)$concordance[1])
CstatL[6]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[6]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

# Brier
BS[6]<-sbrier(obj = Surv(test0$timevent,test0$event), 
              surv_pred0, btime = 3.5)


##### LT = 3 #####


# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_36months.dta")
test0$smoke<-factor(test0$smoke, levels=c("0","1","2"), 
                    ordered=F)
test0$firsttreat<-factor(test0$firsttreat, ordered = F)

pred_prob0 <- predict(m1,newdata=test0,type="response")
test0$pred_LP0 <- predict(m1,newdata=test0, type="link")
surv_pred0 <- 1-pred_prob0

# # Discrimination

# C-index
c0 <- coxph(Surv(timevent,event)~pred_LP0, data=test0)
Cstat[7]<-as.numeric(summary(c0)$concordance[1])
CstatL[7]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[7]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

# Brier
BS[7]<-sbrier(obj = Surv(test0$timevent,test0$event), 
              surv_pred0, btime = 4)
# Performance data frame


LRM.TA_surv<-cbind(LT, Cstat, CstatL, CstatU, BS)
LRM.TA_surv<-as.data.frame(LRM.TA_surv)
write.table(LRM.TA_surv, file = "lrmTA_surv.csv", sep = ",", col.names = NA,
            qmethod = "double")

FPM<-read.dta("baseline_model.dta")
LRM.opt_surv<-bootvadLRM_surv(LR, FPM, 200, 199430)
LRM.opt_surv<-as.data.frame(LRM.opt_surv)
write.table(LRM.opt_surv, file = "lrmOPT_surv.csv", sep = ",", col.names = NA,
            qmethod = "double")




##########################################




bootvadLRM<-function(data, event, B, seed){
  dat <- data
  dat$studyno <- dat$groupid
  
  indiv <- unique(dat$studyno)
  c<-matrix(0, nrow = B, ncol = 7)
  bs<-matrix(0,nrow = B, ncol = 7)
  od.c<-matrix(0, nrow = B, ncol = 7)
  bs.od<-matrix(0, nrow = B, ncol = 7)
  opt.d<-matrix(0, nrow = B, ncol = 7)
  opt.c<-matrix(0, nrow = B, ncol = 7)
  results<-matrix(0, nrow = 2, ncol = 7)
  
  # original data - temporal assessment
  
  od0 <- LR
  od1 <- sixmth
  od2 <- twmth
  od3 <- emth
  od4 <- tfmth
  od5 <- thmth
  od6 <- tsmth
  
  od0$studyno<-od0$groupid
  od1$studyno<-od1$groupid
  od2$studyno<-od2$groupid
  od3$studyno<-od3$groupid
  od4$studyno<-od4$groupid
  od5$studyno<-od5$groupid
  od6$studyno<-od6$groupid
  
  #k10 <- qchisq(0.1,1,lower.tail=FALSE)
  
  for (j in 1:B){
    new_seed = seed +j
    set.seed(new_seed)
    smp <- sort(sample(indiv, length(indiv), replace=TRUE))
    smp.df <- data.frame(studyno=smp)
    print("step")
    boot <- merge(smp.df, dat, by = "studyno", all.x=TRUE)
    
    #fitting logistic regression model 
    print("step")
    stepLRM<-glm(event~age + pgen + firsttreat + 
                    disdur + previous_dmards + lung + 
                    diabetes + bmi + renal + steroids + ovmean + 
                    smoke + dascore, data=boot, family=binomial)
    
    # Temporal assessment sets
    print("step")
    val0<-na.omit(merge(smp.df, od0, by="studyno", all.x=TRUE))
    val1<-na.omit(merge(smp.df, od1, by = "studyno", all.x=TRUE))
    val2<-na.omit(merge(smp.df, od2, by = "studyno", all.x=TRUE))
    val3<-na.omit(merge(smp.df, od3, by = "studyno", all.x=TRUE))
    val4<-na.omit(merge(smp.df, od4, by = "studyno", all.x=TRUE))
    val5<-na.omit(merge(smp.df, od5, by = "studyno", all.x=TRUE))
    val6<-na.omit(merge(smp.df, od6, by = "studyno", all.x=TRUE))
    
    # Survival predictions at different time points
    print("step")
    pred0 <- predict(stepLRM, newdata = val0, type="response")
    pred1 <- predict(stepLRM, newdata = val1, type="response")
    pred2 <- predict(stepLRM, newdata = val2, type="response")
    pred3 <- predict(stepLRM, newdata = val3, type="response")
    pred4 <- predict(stepLRM, newdata = val4, type="response")
    pred5 <- predict(stepLRM, newdata = val5, type="response")
    pred6 <- predict(stepLRM, newdata = val6, type="response")
    print("step")
    lp0 <- predict(stepLRM, newdata = val0, type="link")
    lp1 <- predict(stepLRM, newdata = val1, type="link")
    lp2 <- predict(stepLRM, newdata = val2, type="link")
    lp3 <- predict(stepLRM, newdata = val3, type="link")
    lp4 <- predict(stepLRM, newdata = val4, type="link")
    lp5 <- predict(stepLRM, newdata = val5, type="link")
    lp6 <- predict(stepLRM, newdata = val6, type="link")
    print("step")
    # Discrimination measure
    c0 <- roc(val0$event~pred0)
    c0 <-as.numeric(c0$auc[1])
    c1 <- roc(val1$event~pred1)
    c1 <-as.numeric(c1$auc[1])
    c2 <- roc(val2$event~pred2)
    c2 <-as.numeric(c2$auc[1])
    c3 <- roc(val3$event~pred3)
    c3 <-as.numeric(c3$auc[1])
    c4 <- roc(val4$event~pred4)
    c4 <-as.numeric(c4$auc[1])
    c5 <- roc(val5$event~pred5)
    c5 <-as.numeric(c5$auc[1])
    c6 <- roc(val6$event~pred6)
    c6 <-as.numeric(c6$auc[1])
    print("step")
    c[j,]<-c(c0, c1, c2, c3, c4, c5, c6)

    # Calibration measure
    
    BS0<-BrierScore(resp= val0$event, pred = pred0, 
                      scaled = FALSE)
    BS1<-BrierScore(resp= val1$event, pred = pred1, 
                    scaled = FALSE)
    BS2<-BrierScore(resp= val2$event, pred = pred2, 
                    scaled = FALSE)
    BS3<-BrierScore(resp= val3$event, pred = pred3, 
                    scaled = FALSE)
    BS4<-BrierScore(resp= val4$event, pred = pred4, 
                    scaled = FALSE)
    BS5<-BrierScore(resp= val5$event, pred = pred5, 
                    scaled = FALSE)
    BS6<-BrierScore(resp= val6$event, pred = pred6, 
                    scaled = FALSE)
    
    bs[j,]<-c(BS0, BS1, BS2, BS3, BS4, BS5, BS6)

    # Survival predictions - original data
    
    pred0od <- predict(stepLRM, newdata = od0, type="response")
    pred1od <- predict(stepLRM, newdata = od1, type="response")
    pred2od <- predict(stepLRM, newdata = od2, type="response")
    pred3od <- predict(stepLRM, newdata = od3, type="response")
    pred4od <- predict(stepLRM, newdata = od4, type="response")
    pred5od <- predict(stepLRM, newdata = od5, type="response")
    pred6od <- predict(stepLRM, newdata = od6, type="response")
    
    lp0od <- predict(stepLRM, newdata = od0, type="link")
    lp1od <- predict(stepLRM, newdata = od1, type="link")
    lp2od <- predict(stepLRM, newdata = od2, type="link")
    lp3od <- predict(stepLRM, newdata = od3, type="link")
    lp4od <- predict(stepLRM, newdata = od4, type="link")
    lp5od <- predict(stepLRM, newdata = od5, type="link")
    lp6od <- predict(stepLRM, newdata = od6, type="link")
    
    # Discrimination measure
    
    c0od <- roc(od0$event~pred0od)
    c0od <-as.numeric(c0od$auc[1])
    c1od <- roc(od1$event~pred1od)
    c1od <-as.numeric(c1od$auc[1])
    c2od <- roc(od2$event~pred2od)
    c2od <-as.numeric(c2od$auc[1])
    c3od <- roc(od3$event~pred3od)
    c3od <-as.numeric(c3od$auc[1])
    c4od <- roc(od4$event~pred4od)
    c4od <-as.numeric(c4od$auc[1])
    c5od <- roc(od5$event~pred5od)
    c5od <-as.numeric(c5od$auc[1])
    c6od <- roc(od6$event~pred6od)
    c6od <-as.numeric(c6od$auc[1])
    
    od.c[j,]<-c(c0od, c1od, c2od, c3od, c4od, c5od, c6od)

    # Calibration measure
    
    BS0od<-BrierScore(resp= od0$event, pred = pred0od, 
                    scaled = FALSE)
    BS1od<-BrierScore(resp= od1$event, pred = pred1od, 
                    scaled = FALSE)
    BS2od<-BrierScore(resp= od2$event, pred = pred2od, 
                    scaled = FALSE)
    BS3od<-BrierScore(resp= od3$event, pred = pred3od, 
                    scaled = FALSE)
    BS4od<-BrierScore(resp= od4$event, pred = pred4od, 
                    scaled = FALSE)
    BS5od<-BrierScore(resp= od5$event, pred = pred5od, 
                    scaled = FALSE)
    BS6od<-BrierScore(resp= od6$event, pred = pred6od, 
                    scaled = FALSE)
    
    bs.od[j,]<-c(BS0od, BS1od, BS2od, BS3od, BS4od, BS5od, BS6od)

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


bootvadLRM_surv<-function(data, surv_data, B, seed){
  dat <- data
  dat$studyno <- dat$groupid
  surv_dat <- surv_data
  surv_dat$studyno <- surv_dat$groupid
  
  indiv <- unique(surv_dat$studyno)
  
  c<-matrix(0, nrow = B, ncol = 7)
  bs<-matrix(0,nrow = B, ncol = 7)
  od.c<-matrix(0, nrow = B, ncol = 7)
  bs.od<-matrix(0, nrow = B, ncol = 7)
  opt.d<-matrix(0, nrow = B, ncol = 7)
  opt.c<-matrix(0, nrow = B, ncol = 7)
  results<-matrix(0, nrow = 2, ncol = 7)
  
  # original data - temporal assessment
  
  od0 <- read.dta("baseline_model.dta")
  od0$firsttreat<-as.factor(od0$firsttreat)
  od0$smoke<-as.factor(od0$smoke)
  od1 <- read.dta("validation_6months.dta")
  od1$firsttreat<-as.factor(od1$firsttreat)
  od1$smoke<-as.factor(od1$smoke)
  od2 <- read.dta("validation_12months.dta")
  od2$firsttreat<-as.factor(od2$firsttreat)
  od2$smoke<-as.factor(od2$smoke)
  od3 <- read.dta("validation_18months.dta")
  od3$firsttreat<-as.factor(od3$firsttreat)
  od3$smoke<-as.factor(od3$smoke)
  od4 <- read.dta("validation_24months.dta")
  od4$firsttreat<-as.factor(od4$firsttreat)
  od4$smoke<-as.factor(od4$smoke)
  od5 <- read.dta("validation_30months.dta")
  od5$firsttreat<-as.factor(od5$firsttreat)
  od5$smoke<-as.factor(od5$smoke)
  od6 <- read.dta("validation_36months.dta")
  od6$firsttreat<-as.factor(od6$firsttreat)
  od6$smoke<-as.factor(od6$smoke)
  
  od0$studyno<-od0$groupid
  od1$studyno<-od1$groupid
  od2$studyno<-od2$groupid
  od3$studyno<-od3$groupid
  od4$studyno<-od4$groupid
  od5$studyno<-od5$groupid
  od6$studyno<-od6$groupid
  
  #k10 <- qchisq(0.1,1,lower.tail=FALSE)
  
  for (j in 1:B){
    new_seed = seed +j
    set.seed(new_seed)
    smp <- sort(sample(indiv, length(indiv), replace=TRUE))
    smp.df <- data.frame(studyno=smp)

    boot <- merge(smp.df, dat, by = "studyno", all.x=TRUE)
    boot <- boot %>% filter((timevent == 1 && event == 0) | event == 1)
    
    #fitting logistic regression model 

    stepLRM<-glm(event~age + pgen + firsttreat + 
                   disdur + previous_dmards + lung + 
                   diabetes + bmi + renal + steroids + ovmean + 
                   smoke + dascore, data=boot, family=binomial)
    
    # Temporal assessment sets

    val0<-na.omit(merge(smp.df, od0, by="studyno", all.x=TRUE))
    val1<-na.omit(merge(smp.df, od1, by = "studyno", all.x=TRUE))
    val2<-na.omit(merge(smp.df, od2, by = "studyno", all.x=TRUE))
    val3<-na.omit(merge(smp.df, od3, by = "studyno", all.x=TRUE))
    val4<-na.omit(merge(smp.df, od4, by = "studyno", all.x=TRUE))
    val5<-na.omit(merge(smp.df, od5, by = "studyno", all.x=TRUE))
    val6<-na.omit(merge(smp.df, od6, by = "studyno", all.x=TRUE))
    
    # Survival predictions at different time points

    pred0 <- predict(stepLRM, newdata = val0, type="response")
    pred1 <- predict(stepLRM, newdata = val1, type="response")
    pred2 <- predict(stepLRM, newdata = val2, type="response")
    pred3 <- predict(stepLRM, newdata = val3, type="response")
    pred4 <- predict(stepLRM, newdata = val4, type="response")
    pred5 <- predict(stepLRM, newdata = val5, type="response")
    pred6 <- predict(stepLRM, newdata = val6, type="response")
    surv_pred0 <- 1-pred0
    surv_pred1 <- 1-pred1
    surv_pred2 <- 1-pred2
    surv_pred3 <- 1-pred3
    surv_pred4 <- 1-pred4
    surv_pred5 <- 1-pred5
    surv_pred6 <- 1-pred6

    lp0 <- predict(stepLRM, newdata = val0, type="link")
    lp1 <- predict(stepLRM, newdata = val1, type="link")
    lp2 <- predict(stepLRM, newdata = val2, type="link")
    lp3 <- predict(stepLRM, newdata = val3, type="link")
    lp4 <- predict(stepLRM, newdata = val4, type="link")
    lp5 <- predict(stepLRM, newdata = val5, type="link")
    lp6 <- predict(stepLRM, newdata = val6, type="link")

    # Discrimination measure
    c0 <- coxph(Surv(val0$timevent, val0$event)~lp0)
    c0 <-as.numeric(summary(c0)$concordance[1])
    c1 <- coxph(Surv(val1$timevent, val1$event)~lp1)
    c1 <-as.numeric(summary(c1)$concordance[1])
    c2 <- coxph(Surv(val2$timevent, val2$event)~lp2)
    c2 <-as.numeric(summary(c2)$concordance[1])
    c3 <- coxph(Surv(val3$timevent, val3$event)~lp3)
    c3 <-as.numeric(summary(c3)$concordance[1])
    c4 <- coxph(Surv(val4$timevent, val4$event)~lp4)
    c4 <-as.numeric(summary(c4)$concordance[1])
    c5 <- coxph(Surv(val5$timevent, val5$event)~lp5)
    c5 <-as.numeric(summary(c5)$concordance[1])
    c6 <- coxph(Surv(val6$timevent, val6$event)~lp6)
    c6 <-as.numeric(summary(c6)$concordance[1])

    c[j,]<-c(c0, c1, c2, c3, c4, c5, c6)
    
    # Calibration measure
    
    BS0<-sbrier(obj=Surv(val0$timevent, val0$event), surv_pred0, 
                    btime=1)
    BS1<-sbrier(obj=Surv(val1$timevent, val1$event), surv_pred1, 
                btime=1.5)
    BS2<-sbrier(obj=Surv(val2$timevent, val2$event), surv_pred2, 
                btime=2)
    BS3<-sbrier(obj=Surv(val3$timevent, val3$event), surv_pred3, 
                btime=2.5)
    BS4<-sbrier(obj=Surv(val4$timevent, val4$event), surv_pred4, 
                btime=3)
    BS5<-sbrier(obj=Surv(val5$timevent, val5$event), surv_pred5, 
                btime=3.5)
    BS6<-sbrier(obj=Surv(val6$timevent, val6$event), surv_pred6, 
                btime=4)
    
    bs[j,]<-c(BS0, BS1, BS2, BS3, BS4, BS5, BS6)
    
    # Survival predictions - original data
    
    pred0od <- predict(stepLRM, newdata = od0, type="response")
    pred1od <- predict(stepLRM, newdata = od1, type="response")
    pred2od <- predict(stepLRM, newdata = od2, type="response")
    pred3od <- predict(stepLRM, newdata = od3, type="response")
    pred4od <- predict(stepLRM, newdata = od4, type="response")
    pred5od <- predict(stepLRM, newdata = od5, type="response")
    pred6od <- predict(stepLRM, newdata = od6, type="response")
    
    surv_pred0od <- 1-pred0od
    surv_pred1od <- 1-pred1od
    surv_pred2od <- 1-pred2od
    surv_pred3od <- 1-pred3od
    surv_pred4od <- 1-pred4od
    surv_pred5od <- 1-pred5od
    surv_pred6od <- 1-pred6od
    
    lp0od <- predict(stepLRM, newdata = od0, type="link")
    lp1od <- predict(stepLRM, newdata = od1, type="link")
    lp2od <- predict(stepLRM, newdata = od2, type="link")
    lp3od <- predict(stepLRM, newdata = od3, type="link")
    lp4od <- predict(stepLRM, newdata = od4, type="link")
    lp5od <- predict(stepLRM, newdata = od5, type="link")
    lp6od <- predict(stepLRM, newdata = od6, type="link")
    
    # Discrimination measure
    
    c0od <- coxph(Surv(od0$timevent, od0$event)~lp0od)
    c0od <-as.numeric(summary(c0od)$concordance[1])
    c1od <- coxph(Surv(od1$timevent, od1$event)~lp1od)
    c1od <-as.numeric(summary(c1od)$concordance[1])
    c2od <- coxph(Surv(od2$timevent, od2$event)~lp2od)
    c2od <-as.numeric(summary(c2od)$concordance[1])
    c3od <- coxph(Surv(od3$timevent, od3$event)~lp3od)
    c3od <-as.numeric(summary(c3od)$concordance[1])
    c4od <- coxph(Surv(od4$timevent, od4$event)~lp4od)
    c4od <-as.numeric(summary(c4od)$concordance[1])
    c5od <- coxph(Surv(od5$timevent, od5$event)~lp5od)
    c5od <-as.numeric(summary(c5od)$concordance[1])
    c6od <- coxph(Surv(od6$timevent, od6$event)~lp6od)
    c6od <-as.numeric(summary(c6od)$concordance[1])
    
    od.c[j,]<-c(c0od, c1od, c2od, c3od, c4od, c5od, c6od)
    
    # Calibration measure
    
    BS0od<-sbrier(obj=Surv(od0$timevent, od0$event), surv_pred0od, 
                btime=1)
    BS1od<-sbrier(obj=Surv(od1$timevent, od1$event), surv_pred1od, 
                btime=1.5)
    BS2od<-sbrier(obj=Surv(od2$timevent, od2$event), surv_pred2od, 
                btime=2)
    BS3od<-sbrier(obj=Surv(od3$timevent, od3$event), surv_pred3od, 
                btime=2.5)
    BS4od<-sbrier(obj=Surv(od4$timevent, od4$event), surv_pred4od, 
                btime=3)
    BS5od<-sbrier(obj=Surv(od5$timevent, od5$event), surv_pred5od, 
                btime=3.5)
    BS6od<-sbrier(obj=Surv(od6$timevent, od6$event), surv_pred6od, 
                btime=4)
    
    bs.od[j,]<-c(BS0od, BS1od, BS2od, BS3od, BS4od, BS5od, BS6od)
    
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
FPM <-read.dta("baseline_model.dta")
LRM.opt_surv<-bootvadLRM_surv(LR, FPM, 200, 199430)
