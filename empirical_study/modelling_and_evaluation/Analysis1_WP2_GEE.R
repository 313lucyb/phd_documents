 
library(geepack)
library(foreign)
library(rms)
library(ROCR)
library(pROC)
library(joineR)
library(pROC)
library(plyr)
library(mfp)
library(colorRamps)
library(grDevices)
library(DescTools)
library(dummies)
library(RColorBrewer)

# GLM cross-sectional dataset for nonlinear relationship
set.seed(37546)
cols <-brewer.pal(11, "RdBu")
cdata<-read.dta("cloglog_model.dta")
nonlinear<- mfp(event ~ fp(timevent, df = 4, select = 0.05) + 
                  offset(log(timevent)), family = binomial(link ="cloglog"), 
                data = cdata)
print(nonlinear)
plot(cdata$timevent, cdata$event)
cdata$NL <- I(cdata$timevent^3)+I(cdata$timevent^3*log(cdata$timevent))
cdata$pred <- predict(nonlinear, type = "response")
plot(cdata$timevent, cdata$pred)
plot(cdata$NL, cdata$pred, col = cols[2], pch=20, 
     xlab = "(Event time)^3 + [(Event time)^3]*log(Event time)", 
     ylab = "Probability of at least 1 event per year", cex.lab = 1.2)


# GEE longitudinal dataset

Gdata<-read.dta("long_GEE.dta")
Gdata$timevent<-Gdata$end
Gdata$smoke<-factor(Gdata$smoke)
Gdata$trtment<-factor(Gdata$firsttreat)
Gdata$roundtime<-as.factor(Gdata$roundtime)
Gdata$NL<-I(Gdata$timevent^3)+I(Gdata$timevent^3*log(Gdata$timevent))
set.seed(7765)

# complementary log-log model without time since baseline

start_time1<-Sys.time()
gee1<-geeglm(serinf~ age + pgen + previous_dmards + 
                     disdur + trtment + ovmean + dascore + steroids +
                     bmi + renal + lung + diabetes + smoke +  offset(log(interval)), 
                         family=binomial(link = "cloglog"), data=Gdata, 
                         id=groupid, waves=pyears, corstr = "exchangeable", 
                         std.err="san.se")
end_time1<-Sys.time()
fit_time1<- end_time1 - start_time1

# complementary log-log model with time since baseline (categorical)

start_time2<-Sys.time()
gee2<-geeglm(serinf~ roundtime + age + pgen + previous_dmards + 
               disdur + trtment + ovmean + dascore + steroids +
               bmi + renal + lung + diabetes + smoke +  offset(log(interval)), 
             family=binomial(link = "cloglog"), data=Gdata, 
             id=groupid, waves=pyears, corstr = "exchangeable", 
             std.err="san.se")
end_time2<-Sys.time()
fit_time2<- end_time2 - start_time2

# complementary log-log with linear continuous variable

start_time3<-Sys.time()
gee3<-geeglm(serinf~ timevent + age + pgen + previous_dmards + 
               disdur + trtment + ovmean + dascore + steroids +
               bmi + renal + lung + diabetes + smoke + 
                offset(log(interval)), 
             family=binomial(link = "cloglog"), data=Gdata, 
             id=groupid, waves=pyears, corstr = "exchangeable", 
             std.err="san.se")
end_time3<-Sys.time()
fit_time3<- end_time3 - start_time3

# complementary log-log with nonlinear continuous variable

start_time4<-Sys.time()
gee4<-geeglm(serinf~ I(timevent^3)+I(timevent^3*log(timevent)) + 
               age + pgen + previous_dmards + 
               disdur + trtment + ovmean + dascore + steroids +
               bmi + renal + lung + diabetes + smoke +  offset(log(interval)), 
             family=binomial(link = "cloglog"), data=Gdata, 
             id=groupid, waves=pyears, corstr = "exchangeable", 
             std.err="san.se")
end_time4<-Sys.time()
fit_time4<- end_time4 - start_time4
fit_times <- c(fit_time1, fit_time2, fit_time3, fit_time4)

## Performance assessment

gee1coef<-as.numeric(gee1$coefficients)
gee2coef<-as.numeric(gee2$coefficients)
gee3coef<-as.numeric(gee3$coefficients)
gee4coef<-as.numeric(gee4$coefficients)

## ## LT = 0 

test0<-read.dta("baseline_binary.dta")
test0$trtment <- test0$firsttreat
test0$trtment<-factor(test0$trtment)
test0$smoke<-factor(test0$smoke)
test0$prev_sinf<-0
test0$roundtime<-rep(1, nrow(test0))
test0$roundtime<-as.factor(test0$roundtime)

test0$interval<-rep(1, nrow(test0))
test0$timevent<-rep(1, nrow(test0))
test0$timevent1<-I(test0$timevent^3)
test0$timevent2<-I(test0$timevent^3*log(test0$timevent))

test0$one<-rep(1, nrow(test0))
test0$roundtime2<-rep(0, nrow(test0))
test0$roundtime3<-rep(0, nrow(test0))
test0$roundtime4<-rep(0, nrow(test0))

#Gee1 - No time since baseline covariate

lpcov10<-test0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat10<-as.matrix(lpcov10)
lp10<-gee1coef %*% t(lpmat10)
lp10<-as.vector(t(lp10))
  
#Gee2 - Categorical time since baseline covariate

lpcov20<- test0 %>% select(one, roundtime2, roundtime3, roundtime4, age, pgen, previous_dmards, 
                           disdur, trt2, trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
  
lpmat20<-as.matrix(lpcov20)
lp20<-gee2coef %*% t(lpmat20)
lp20<-as.vector(t(lp20))

#Gee3 - Linear continuous time since baseline

lpcov30<-test0 %>% select(one, timevent, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat30<-as.matrix(lpcov30)
lp30<-gee3coef %*% t(lpmat30)
lp30<-as.vector(t(lp30))

#Gee4 - Nonlinear continuous time since baseline

lpcov40<-test0 %>% select(one, timevent1, timevent2, age, pgen, previous_dmards, disdur, trt2, 
                          trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat40<-as.matrix(lpcov40)
lp40<-gee4coef %*% t(lpmat40)
lp40<-as.vector(t(lp40))

## Probability of at least one serious infection per year

p_start <- Sys.time()
test0$prob1<-1-exp(-1*exp(lp10))
p_end <- Sys.time()
p2_start <- Sys.time()
test0$prob2<-1-exp(-1*exp(lp20))
p2_end <- Sys.time()
p3_start<-Sys.time()
test0$prob3<-1-exp(-1*exp(lp30))
p3_end<- Sys.time()
p4_start <- Sys.time()
test0$prob4<-1-exp(-1*exp(lp40))
p4_end <- Sys.time()

p_time0 <- p_end - p_start
p2_time0 <- p2_end - p2_start
p3_time0 <- p3_end - p3_start
p4_time0 <- p4_end - p4_start

pred_times0 <- c(p_time0, p2_time0, p3_time0, p4_time0)

BS0<-0
BS0[1]<-BrierScore(resp= test0$event, pred = test0$prob1, 
                  scaled = FALSE)
BS0[2]<-BrierScore(resp= test0$event, pred = test0$prob2, 
                   scaled = FALSE)
BS0[3]<-BrierScore(resp= test0$event, pred = test0$prob3, 
                   scaled = FALSE)
BS0[4]<-BrierScore(resp= test0$event, pred = test0$prob4, 
                   scaled = FALSE)

# C statistic / AUC
c0<-0
c0.L<-0
c0.U<-0

C1 <- roc(test0$event~test0$prob1,ci=TRUE)
c0[1]<-as.numeric(C1$auc)
c0.L[1]<-as.numeric(C1$ci[1])
c0.U[1]<-as.numeric(C1$ci[3])
C2 <- roc(test0$event~test0$prob2,ci=TRUE)
c0[2]<-as.numeric(C2$auc)
c0.L[2]<-as.numeric(C2$ci[1])
c0.U[2]<-as.numeric(C2$ci[3])
C3 <- roc(test0$event~test0$prob3,ci=TRUE)
c0[3]<-as.numeric(C3$auc)
c0.L[3]<-as.numeric(C3$ci[1])
c0.U[3]<-as.numeric(C3$ci[3])
C4 <- roc(test0$event~test0$prob4,ci=TRUE)
c0[4]<-as.numeric(C4$auc)
c0.L[4]<-as.numeric(C4$ci[1])
c0.U[4]<-as.numeric(C4$ci[3])

#  Calibration

# Calibration In The Large
set.seed(112)

Clge0<-0
Clge0.L<-0
Clge0.U<-0

m1 <- glm(test0$event~offset(lp10), family="binomial")
Clge0[1] <- as.numeric(m1$coef)
x10<-confint(m1)
Clge0.L[1]<-as.numeric(x10[1])
Clge0.U[1]<-as.numeric(x10[2])

m2 <- glm(test0$event~offset(lp20), family="binomial")
Clge0[2] <- as.numeric(m2$coef)
x20<-confint(m2)
Clge0.L[2]<-as.numeric(x20[1])
Clge0.U[2]<-as.numeric(x20[2])

m3 <- glm(test0$event~offset(lp30), family="binomial")
Clge0[3] <- as.numeric(m3$coef)
x30<-confint(m3)
Clge0.L[3]<-as.numeric(x30[1])
Clge0.U[3]<-as.numeric(x30[2])

m4 <- glm(test0$event~offset(lp40), family="binomial")
Clge0[4] <- as.numeric(m4$coef)
x40<-confint(m4)
Clge0.L[4]<-as.numeric(x40[1])
Clge0.U[4]<-as.numeric(x40[2])


# Calibration slope
set.seed(182)

Cslope0<-0
Cslope0.L<-0
Cslope0.U<-0

m11 <- glm(test0$event~lp10,family="binomial", x=TRUE,y=TRUE)
Cslope0[1]<-m11$coefficients[2]
s10<-confint(m11)
Cslope0.L[1]<-as.numeric(s10[2,1])
Cslope0.U[1]<-as.numeric(s10[2,2])

m12 <- glm(test0$event~lp20,family="binomial", x=TRUE,y=TRUE)
Cslope0[2]<-m12$coefficients[2]
s20<-confint(m12)
Cslope0.L[2]<-as.numeric(s20[2,1])
Cslope0.U[2]<-as.numeric(s20[2,2])

m13 <- glm(test0$event~lp30,family="binomial", x=TRUE,y=TRUE)
Cslope0[3]<-m13$coefficients[2]
s30<-confint(m13)
Cslope0.L[3]<-as.numeric(s30[2,1])
Cslope0.U[3]<-as.numeric(s30[2,2])

m14 <- glm(test0$event~lp40,family="binomial", x=TRUE,y=TRUE)
Cslope0[4]<-m14$coefficients[2]
s40<-confint(m14)
Cslope0.L[4]<-as.numeric(s40[2,1])
Cslope0.U[4]<-as.numeric(s40[2,2])

par(mfrow=c(2,2))
par(pty="s")

# 1 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test0$prob1,breaks=quantile(test0$prob1, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test0,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob1))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1], lwd = 2)
}
h <- hist(test0$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test0$event))~test0$prob1,span=1))
lines_data <- data.frame(test0$prob1,obs_all)
lines_data2 <- lines_data[order(test0$prob1),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 2 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test0$prob2,breaks=quantile(test0$prob2, 
                                          prob = 
                                            c(0,0.1,0.2,0.3,0.4,0.5,
                                              0.6,0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test0,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob2))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[2],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[2], lwd = 2)
}
h <- hist(test0$prob2, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test0$event))~test0$prob2,span=1))
lines_data <- data.frame(test0$prob2,obs_all)
lines_data2 <- lines_data[order(test0$prob2),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),
       lty=c(0,2,1,1),pch=c(16,NA,NA,NA),bty="n")

# 3 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test0$prob3,breaks=quantile(test0$prob3, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,
                                                   0.4,0.5,0.6,
                                                   0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of 
#patients in each risk group 
gpdata <- cbind(test0,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob3))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[3],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[3], lwd = 2)
}
h <- hist(test0$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test0$event))~test0$prob3,span=1))
lines_data <- data.frame(test0$prob3,obs_all)
lines_data2 <- lines_data[order(test0$prob3),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 4 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test0$prob4,breaks=quantile(test0$prob4, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test0,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob4))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[4],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[4], lwd = 2)
}
h <- hist(test0$prob4, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test0$event))~test0$prob4,span=1))
lines_data <- data.frame(test0$prob4,obs_all)
lines_data2 <- lines_data[order(test0$prob4),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[4],"black",cols[4],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


LT0<-c(0,0,0,0)
Mod<-c(1,2,3,4)
LT.0<-cbind(LT0, Mod, BS0, c0, c0.L, c0.U, Clge0, 
            Clge0.L, Clge0.U,Cslope0, Cslope0.L,
            Cslope0.U, fit_times, pred_times0)



## ## LT = 0.5 

test1<-read.dta("validation_6binary.dta")
test1$trtment<-factor(test1$trtment)
test1$smoke<-factor(test1$smoke)
test1$roundtime<-rep(2, nrow(test1))
test1$roundtime<-as.factor(test1$roundtime)
test1$interval<-rep(1, nrow(test1))
test1$timevent<-rep(1.5, nrow(test1))
test1$timevent1<-I(test1$timevent^3)
test1$timevent2<-I(test1$timevent^3*log(test1$timevent))

test1$one<-rep(1, nrow(test1))
test1$roundtime2<-rep(0, nrow(test1))
test1$roundtime3<-rep(0, nrow(test1))
test1$roundtime4<-rep(0, nrow(test1))

#Gee1 - No time since baseline covariate

lpcov11<-test1 %>% select(one, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat11<-as.matrix(lpcov11)
lp11<-gee1coef %*% t(lpmat11)
lp11<-as.vector(t(lp11))

#Gee2 - Categorical time since baseline covariate

lpcov21<-test1 %>% select(one, roundtime2, roundtime3, roundtime4, age, pgen, previous_dmards, 
                          disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat21<-as.matrix(lpcov21)
lp21<-gee2coef %*% t(lpmat21)
lp21<-as.vector(t(lp21))

#Gee3 - Linear continuous time since baseline

lpcov31<-test1 %>% select(one, timevent, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat31<-as.matrix(lpcov31)
lp31<-gee3coef %*% t(lpmat31)
lp31<-as.vector(t(lp31))


#Gee4 - Nonlinear continuous time since baseline

lpcov41<-test1 %>% select(one, timevent1, timevent2, age, pgen, previous_dmards, disdur, trt2, 
                          trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat41<-as.matrix(lpcov41)
lp41<-gee4coef %*% t(lpmat41)
lp41<-as.vector(t(lp41))


## Probability of at least one serious infection per year
p_start <- Sys.time()
test1$prob1<-1-exp(-1*exp(lp11))
p_end <- Sys.time()
p_start2 <- Sys.time()
test1$prob2<-1-exp(-1*exp(lp21))
p_end2 <- Sys.time()
p_start3 <- Sys.time()
test1$prob3<-1-exp(-1*exp(lp31))
p_end3 <- Sys.time()
p_start4 <- Sys.time()
test1$prob4<-1-exp(-1*exp(lp41))
p_end4 <- Sys.time()

p_time1 <- p_end - p_start
p2_time1 <- p2_end - p2_start
p3_time1 <- p3_end - p3_start
p4_time1 <- p4_end - p4_start

pred_times1 <- c(p_time1, p2_time1, p3_time1, p4_time1)

BS1<-0
BS1[1]<-BrierScore(resp= test1$event, pred = test1$prob1, 
                   scaled = FALSE)
BS1[2]<-BrierScore(resp= test1$event, pred = test1$prob2, 
                   scaled = FALSE)
BS1[3]<-BrierScore(resp= test1$event, pred = test1$prob3, 
                   scaled = FALSE)
BS1[4]<-BrierScore(resp= test1$event, pred = test1$prob4, 
                   scaled = FALSE)



# C statistic / AUC
c1<-0
c1.L<-0
c1.U<-0

C1 <- roc(test1$event~test1$prob1,ci=TRUE)
c1[1]<-as.numeric(C1$auc)
c1.L[1]<-as.numeric(C1$ci[1])
c1.U[1]<-as.numeric(C1$ci[3])
C2 <- roc(test1$event~test1$prob2,ci=TRUE)
c1[2]<-as.numeric(C2$auc)
c1.L[2]<-as.numeric(C2$ci[1])
c1.U[2]<-as.numeric(C2$ci[3])
C3 <- roc(test1$event~test1$prob3,ci=TRUE)
c1[3]<-as.numeric(C3$auc)
c1.L[3]<-as.numeric(C3$ci[1])
c1.U[3]<-as.numeric(C3$ci[3])
C4 <- roc(test1$event~test1$prob4,ci=TRUE)
c1[4]<-as.numeric(C4$auc)
c1.L[4]<-as.numeric(C4$ci[1])
c1.U[4]<-as.numeric(C4$ci[3])

#  Calibration

# Calibration In The Large
set.seed(112)

Clge1<-0
Clge1.L<-0
Clge1.U<-0

m1 <- glm(test1$event~offset(lp11), family="binomial")
Clge1[1] <- as.numeric(m1$coef)
x11<-confint(m1)
Clge1.L[1]<-as.numeric(x11[1])
Clge1.U[1]<-as.numeric(x11[2])

m2 <- glm(test1$event~offset(lp21), family="binomial")
Clge1[2] <- as.numeric(m2$coef)
x21<-confint(m2)
Clge1.L[2]<-as.numeric(x21[1])
Clge1.U[2]<-as.numeric(x21[2])

m3 <- glm(test1$event~offset(lp31), family="binomial")
Clge1[3] <- as.numeric(m3$coef)
x31<-confint(m3)
Clge1.L[3]<-as.numeric(x31[1])
Clge1.U[3]<-as.numeric(x31[2])

m4 <- glm(test1$event~offset(lp41), family="binomial")
Clge1[4] <- as.numeric(m4$coef)
x41<-confint(m4)
Clge1.L[4]<-as.numeric(x41[1])
Clge1.U[4]<-as.numeric(x41[2])


# Calibration slope
set.seed(182)

Cslope1<-0
Cslope1.L<-0
Cslope1.U<-0

m11 <- glm(test1$event~lp11,family="binomial", x=TRUE,y=TRUE)
Cslope1[1]<-m11$coefficients[2]
s11<-confint(m11)
Cslope1.L[1]<-as.numeric(s11[2,1])
Cslope1.U[1]<-as.numeric(s11[2,2])

m12 <- glm(test1$event~lp21,family="binomial", x=TRUE,y=TRUE)
Cslope1[2]<-m12$coefficients[2]
s21<-confint(m12)
Cslope1.L[2]<-as.numeric(s21[2,1])
Cslope1.U[2]<-as.numeric(s21[2,2])

m13 <- glm(test1$event~lp31,family="binomial", x=TRUE,y=TRUE)
Cslope1[3]<-m13$coefficients[2]
s31<-confint(m13)
Cslope1.L[3]<-as.numeric(s31[2,1])
Cslope1.U[3]<-as.numeric(s31[2,2])

m14 <- glm(test1$event~lp41,family="binomial", x=TRUE,y=TRUE)
Cslope1[4]<-m14$coefficients[2]
s41<-confint(m14)
Cslope1.L[4]<-as.numeric(s41[2,1])
Cslope1.U[4]<-as.numeric(s41[2,2])

par(mfrow=c(2,2))
par(pty="s")

# 1 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test1$prob1,breaks=quantile(test1$prob1, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test1,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob1))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1], lwd = 2)
}
h <- hist(test1$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test1$event))~test1$prob1,span=1))
lines_data <- data.frame(test1$prob1,obs_all)
lines_data2 <- lines_data[order(test1$prob1),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 2 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test1$prob2,breaks=quantile(test1$prob2, 
                                          prob = 
                                            c(0,0.1,0.2,0.3,0.4,0.5,
                                              0.6,0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test1,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob2))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[2],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[2], lwd = 2)
}
h <- hist(test1$prob2, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test1$event))~test1$prob2,span=1))
lines_data <- data.frame(test1$prob2,obs_all)
lines_data2 <- lines_data[order(test1$prob2),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),
       lty=c(0,2,1,1),pch=c(16,NA,NA,NA),bty="n")

# 3 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test1$prob3,breaks=quantile(test1$prob3, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,
                                                   0.4,0.5,0.6,
                                                   0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of 
#patients in each risk group 
gpdata <- cbind(test1,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob3))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[3],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[3], lwd = 2)
}
h <- hist(test1$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test1$event))~test1$prob3,span=1))
lines_data <- data.frame(test1$prob3,obs_all)
lines_data2 <- lines_data[order(test1$prob3),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 4 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test1$prob4,breaks=quantile(test1$prob4, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test1,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob4))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[4],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[4], lwd = 2)
}
h <- hist(test1$prob4, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test1$event))~test1$prob4,span=1))
lines_data <- data.frame(test1$prob4,obs_all)
lines_data2 <- lines_data[order(test1$prob4),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[4],"black",cols[4],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

LT1<-c(0.5,0.5,0.5,0.5)
Mod<-c(1,2,3,4)
LT.1<-cbind(LT1, Mod, BS1, c1, c1.L, c1.U, Clge1, 
            Clge1.L, Clge1.U,Cslope1, Cslope1.L,
            Cslope1.U, fit_times, pred_times1)



## ## LT = 1 

test2<-read.dta("validation_12binary.dta")
test2$trtment<-factor(test2$trtment)
test2$smoke<-factor(test2$smoke)
test2$interval<-rep(1, nrow(test2))
test2$timevent<-rep(2, nrow(test2))
test2$timevent1<-I(test2$timevent^3)
test2$timevent2<-I(test2$timevent^3*log(test2$timevent))

test2$one<-rep(1, nrow(test2))
test2$roundtime2<-rep(1, nrow(test2))
test2$roundtime3<-rep(0, nrow(test2))
test2$roundtime4<-rep(0, nrow(test2))

#Gee1 - No time since baseline covariate

lpcov12<-test2 %>% select(one, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat12<-as.matrix(lpcov12)
lp12<-gee1coef %*% t(lpmat12)
lp12<-as.vector(t(lp12))

#Gee2 - Categorical time since baseline covariate

lpcov22<-test2 %>% select(one, roundtime2, roundtime3, roundtime4, age, pgen, previous_dmards, 
                          disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat22<-as.matrix(lpcov22)
lp22<-gee2coef %*% t(lpmat22)
lp22<-as.vector(t(lp22))

#Gee3 - Linear continuous time since baseline

lpcov32<-test2 %>% select(one, timevent, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat32<-as.matrix(lpcov32)
lp32<-gee3coef %*% t(lpmat32)
lp32<-as.vector(t(lp32))


#Gee4 - Nonlinear continuous time since baseline

lpcov42<-test2 %>% select(one, timevent1, timevent2, age, pgen, previous_dmards, disdur, trt2, 
                          trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat42<-as.matrix(lpcov42)
lp42<-gee4coef %*% t(lpmat42)
lp42<-as.vector(t(lp42))

## Probability of at least one serious infection per year
p_start <- Sys.time()
test2$prob1<-1-exp(-1*exp(lp12))
p_end <- Sys.time()
p_start2 <- Sys.time()
test2$prob2<-1-exp(-1*exp(lp22))
p_end2 <- Sys.time()
p_start3 <- Sys.time()
test2$prob3<-1-exp(-1*exp(lp32))
p_end3 <- Sys.time()
p_start4 <- Sys.time()
test2$prob4<-1-exp(-1*exp(lp42))
p_end4 <- Sys.time()

p_time2 <- p_end - p_start
p2_time2 <- p2_end - p2_start
p3_time2 <- p3_end - p3_start
p4_time2 <- p4_end - p4_start

pred_times2 <- c(p_time2, p2_time2, p3_time2, p4_time2)

BS2<-0
BS2[1]<-BrierScore(resp= test2$event, pred = test2$prob1, 
                   scaled = FALSE)
BS2[2]<-BrierScore(resp= test2$event, pred = test2$prob2, 
                   scaled = FALSE)
BS2[3]<-BrierScore(resp= test2$event, pred = test2$prob3, 
                   scaled = FALSE)
BS2[4]<-BrierScore(resp= test2$event, pred = test2$prob4, 
                   scaled = FALSE)

# C statistic / AUC
c2<-0
c2.L<-0
c2.U<-0

C1 <- roc(test2$event~test2$prob1,ci=TRUE)
c2[1]<-as.numeric(C1$auc)
c2.L[1]<-as.numeric(C1$ci[1])
c2.U[1]<-as.numeric(C1$ci[3])
C2 <- roc(test2$event~test2$prob2,ci=TRUE)
c2[2]<-as.numeric(C2$auc)
c2.L[2]<-as.numeric(C2$ci[1])
c2.U[2]<-as.numeric(C2$ci[3])
C3 <- roc(test2$event~test2$prob3,ci=TRUE)
c2[3]<-as.numeric(C3$auc)
c2.L[3]<-as.numeric(C3$ci[1])
c2.U[3]<-as.numeric(C3$ci[3])
C4 <- roc(test2$event~test2$prob4,ci=TRUE)
c2[4]<-as.numeric(C4$auc)
c2.L[4]<-as.numeric(C4$ci[1])
c2.U[4]<-as.numeric(C4$ci[3])

#  Calibration

# Calibration In The Large
set.seed(112)

Clge2<-0
Clge2.L<-0
Clge2.U<-0

m1 <- glm(test2$event~offset(lp12), family="binomial")
Clge2[1] <- as.numeric(m1$coef)
x12<-confint(m1)
Clge2.L[1]<-as.numeric(x12[1])
Clge2.U[1]<-as.numeric(x12[2])

m2 <- glm(test2$event~offset(lp22), family="binomial")
Clge2[2] <- as.numeric(m2$coef)
x22<-confint(m2)
Clge2.L[2]<-as.numeric(x22[1])
Clge2.U[2]<-as.numeric(x22[2])

m3 <- glm(test2$event~offset(lp32), family="binomial")
Clge2[3] <- as.numeric(m3$coef)
x32<-confint(m3)
Clge2.L[3]<-as.numeric(x32[1])
Clge2.U[3]<-as.numeric(x32[2])

m4 <- glm(test2$event~offset(lp42), family="binomial")
Clge2[4] <- as.numeric(m4$coef)
x42<-confint(m4)
Clge2.L[4]<-as.numeric(x42[1])
Clge2.U[4]<-as.numeric(x42[2])


# Calibration slope
set.seed(182)

Cslope2<-0
Cslope2.L<-0
Cslope2.U<-0

m21 <- glm(test2$event~lp12,family="binomial", x=TRUE,y=TRUE)
Cslope2[1]<-m21$coefficients[2]
s12<-confint(m21)
Cslope2.L[1]<-as.numeric(s12[2,1])
Cslope2.U[1]<-as.numeric(s12[2,2])

m22 <- glm(test2$event~lp22,family="binomial", x=TRUE,y=TRUE)
Cslope2[2]<-m22$coefficients[2]
s22<-confint(m22)
Cslope2.L[2]<-as.numeric(s22[2,1])
Cslope2.U[2]<-as.numeric(s22[2,2])

m23 <- glm(test2$event~lp32,family="binomial", x=TRUE,y=TRUE)
Cslope2[3]<-m23$coefficients[2]
s32<-confint(m23)
Cslope2.L[3]<-as.numeric(s32[2,1])
Cslope2.U[3]<-as.numeric(s32[2,2])

m24 <- glm(test2$event~lp42,family="binomial", x=TRUE,y=TRUE)
Cslope2[4]<-m24$coefficients[2]
s42<-confint(m24)
Cslope2.L[4]<-as.numeric(s42[2,1])
Cslope2.U[4]<-as.numeric(s41[2,2])

par(mfrow=c(2,2))
par(pty="s")

# 1 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test2$prob1,breaks=quantile(test2$prob1, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test2,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob1))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1], lwd = 2)
}
h <- hist(test2$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test2$event))~test2$prob1,span=1))
lines_data <- data.frame(test2$prob1,obs_all)
lines_data2 <- lines_data[order(test2$prob1),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 2 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test2$prob2,breaks=quantile(test2$prob2, 
                                          prob = 
                                            c(0,0.1,0.2,0.3,0.4,0.5,
                                              0.6,0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test2,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob2))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[2],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[2], lwd = 2)
}
h <- hist(test2$prob2, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test2$event))~test2$prob2,span=1))
lines_data <- data.frame(test2$prob2,obs_all)
lines_data2 <- lines_data[order(test2$prob2),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),
       lty=c(0,2,1,1),pch=c(16,NA,NA,NA),bty="n")

# 3 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test2$prob3,breaks=quantile(test2$prob3, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,
                                                   0.4,0.5,0.6,
                                                   0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of 
#patients in each risk group 
gpdata <- cbind(test2,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob3))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[3],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[3], lwd = 2)
}
h <- hist(test2$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test2$event))~test2$prob3,span=1))
lines_data <- data.frame(test2$prob3,obs_all)
lines_data2 <- lines_data[order(test2$prob3),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 4 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test2$prob4,breaks=quantile(test2$prob4, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test2,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob4))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[4],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[4], lwd = 2)
}
h <- hist(test2$prob4, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test2$event))~test2$prob4,span=1))
lines_data <- data.frame(test2$prob4,obs_all)
lines_data2 <- lines_data[order(test2$prob4),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[4],"black",cols[4],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


LT2<-c(1,1,1,1)
Mod<-c(1,2,3,4)
LT.2<-cbind(LT2, Mod, BS2, c2, c2.L, c2.U, Clge2, 
            Clge2.L, Clge2.U,Cslope2, Cslope2.L,
            Cslope2.U, fit_times, pred_times2)



## ## LT = 1.5 

test3<-read.dta("validation_18binary.dta")
test3$trtment<-factor(test3$trtment)
test3$smoke<-factor(test3$smoke)
test3$interval<-rep(1, nrow(test3))
test3$timevent<-rep(2.5, nrow(test3))
test3$timevent1<-I(test3$timevent^3)
test3$timevent2<-I(test3$timevent^3*log(test3$timevent))

test3$one<-rep(1, nrow(test3))
test3$roundtime2<-rep(1, nrow(test3))
test3$roundtime3<-rep(0, nrow(test3))
test3$roundtime4<-rep(0, nrow(test3))

#Gee1 - No time since baseline covariate

lpcov13<-test3 %>% select(one, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat13<-as.matrix(lpcov13)
lp13<-gee1coef %*% t(lpmat13)
lp13<-as.vector(t(lp13))

#Gee2 - Categorical time since baseline covariate

lpcov23<-test3 %>% select(one, roundtime2, roundtime3, roundtime4, age, pgen, previous_dmards, 
                          disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat23<-as.matrix(lpcov23)
lp23<-gee2coef %*% t(lpmat23)
lp23<-as.vector(t(lp23))

#Gee3 - Linear continuous time since baseline

lpcov33<-test3 %>% select(one, timevent, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat33<-as.matrix(lpcov33)
lp33<-gee3coef %*% t(lpmat33)
lp33<-as.vector(t(lp33))


#Gee4 - Nonlinear continuous time since baseline

lpcov43<- test3 %>% select(one, timevent1, timevent2, age, pgen, previous_dmards, disdur, trt2, 
                          trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat43<-as.matrix(lpcov43)
lp43<-gee4coef %*% t(lpmat43)
lp43<-as.vector(t(lp43))

## Probability of at least one serious infection per year
p_start <- Sys.time()
test3$prob1<-1-exp(-1*exp(lp13))
p_end <- Sys.time()
p_start2 <- Sys.time()
test3$prob2<-1-exp(-1*exp(lp23))
p_end2 <- Sys.time()
p_start3 <- Sys.time()
test3$prob3<-1-exp(-1*exp(lp33))
p_end3 <- Sys.time()
p_start4 <- Sys.time()
test3$prob4<-1-exp(-1*exp(lp43))
p_end4 <- Sys.time()

p_time3 <- p_end - p_start
p2_time3 <- p2_end - p2_start
p3_time3 <- p3_end - p3_start
p4_time3 <- p4_end - p4_start

pred_times3 <- c(p_time3, p2_time3, p3_time3, p4_time3)

BS3<-0
BS3[1]<-BrierScore(resp= test3$event, pred = test3$prob1, 
                   scaled = FALSE)
BS3[2]<-BrierScore(resp= test3$event, pred = test3$prob2, 
                   scaled = FALSE)
BS3[3]<-BrierScore(resp= test3$event, pred = test3$prob3, 
                   scaled = FALSE)
BS3[4]<-BrierScore(resp= test3$event, pred = test3$prob4, 
                   scaled = FALSE)

# C statistic / AUC
c3<-0
c3.L<-0
c3.U<-0

C1 <- roc(test3$event~test3$prob1,ci=TRUE)
c3[1]<-as.numeric(C1$auc)
c3.L[1]<-as.numeric(C1$ci[1])
c3.U[1]<-as.numeric(C1$ci[3])
C2 <- roc(test3$event~test3$prob2,ci=TRUE)
c3[2]<-as.numeric(C2$auc)
c3.L[2]<-as.numeric(C2$ci[1])
c3.U[2]<-as.numeric(C2$ci[3])
C3 <- roc(test3$event~test3$prob3,ci=TRUE)
c3[3]<-as.numeric(C3$auc)
c3.L[3]<-as.numeric(C3$ci[1])
c3.U[3]<-as.numeric(C3$ci[3])
C4 <- roc(test3$event~test3$prob4,ci=TRUE)
c3[4]<-as.numeric(C4$auc)
c3.L[4]<-as.numeric(C4$ci[1])
c3.U[4]<-as.numeric(C4$ci[3])

#  Calibration

# Calibration In The Large
set.seed(112)

Clge3<-0
Clge3.L<-0
Clge3.U<-0

m1 <- glm(test3$event~offset(lp13), family="binomial")
Clge3[1] <- as.numeric(m1$coef)
x13<-confint(m1)
Clge3.L[1]<-as.numeric(x13[1])
Clge3.U[1]<-as.numeric(x13[2])

m2 <- glm(test3$event~offset(lp23), family="binomial")
Clge3[2] <- as.numeric(m2$coef)
x23<-confint(m2)
Clge3.L[2]<-as.numeric(x23[1])
Clge3.U[2]<-as.numeric(x23[2])

m3 <- glm(test3$event~offset(lp33), family="binomial")
Clge3[3] <- as.numeric(m3$coef)
x33<-confint(m3)
Clge3.L[3]<-as.numeric(x33[1])
Clge3.U[3]<-as.numeric(x33[2])

m4 <- glm(test3$event~offset(lp43), family="binomial")
Clge3[4] <- as.numeric(m4$coef)
x43<-confint(m4)
Clge3.L[4]<-as.numeric(x43[1])
Clge3.U[4]<-as.numeric(x43[2])


# Calibration slope
set.seed(182)

Cslope3<-0
Cslope3.L<-0
Cslope3.U<-0

m31 <- glm(test3$event~lp13,family="binomial", x=TRUE,y=TRUE)
Cslope3[1]<-m31$coefficients[2]
s13<-confint(m31)
Cslope3.L[1]<-as.numeric(s13[2,1])
Cslope3.U[1]<-as.numeric(s13[2,2])

m32 <- glm(test3$event~lp23,family="binomial", x=TRUE,y=TRUE)
Cslope3[2]<-m32$coefficients[2]
s23<-confint(m32)
Cslope3.L[2]<-as.numeric(s23[2,1])
Cslope3.U[2]<-as.numeric(s23[2,2])

m33 <- glm(test3$event~lp33,family="binomial", x=TRUE,y=TRUE)
Cslope3[3]<-m33$coefficients[2]
s33<-confint(m33)
Cslope3.L[3]<-as.numeric(s33[2,1])
Cslope3.U[3]<-as.numeric(s33[2,2])

m34 <- glm(test3$event~lp43,family="binomial", x=TRUE,y=TRUE)
Cslope3[4]<-m34$coefficients[2]
s43<-confint(m34)
Cslope3.L[4]<-as.numeric(s43[2,1])
Cslope3.U[4]<-as.numeric(s43[2,2])

par(mfrow=c(2,2))
par(pty="s")

# 1 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test3$prob1,breaks=quantile(test3$prob1, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test3,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob1))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch =16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1], lwd = 2)
}
h <- hist(test3$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test3$event))~test3$prob1,span=1))
lines_data <- data.frame(test3$prob1,obs_all)
lines_data2 <- lines_data[order(test3$prob1),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 2 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test3$prob2,breaks=quantile(test3$prob2, 
                                          prob = 
                                            c(0,0.1,0.2,0.3,0.4,0.5,
                                              0.6,0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test3,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob2))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[2],
     ylab="Observed",xlab="Expected", pch =16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[2], lwd = 2)
}
h <- hist(test3$prob2, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test3$event))~test3$prob2,span=1))
lines_data <- data.frame(test3$prob2,obs_all)
lines_data2 <- lines_data[order(test3$prob2),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),
       lty=c(0,2,1,1),pch=c(16,NA,NA,NA),bty="n")

# 3 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test3$prob3,breaks=quantile(test3$prob3, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,
                                                   0.4,0.5,0.6,
                                                   0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of 
#patients in each risk group 
gpdata <- cbind(test3,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob3))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[3],
     ylab="Observed",xlab="Expected", pch =16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[3], lwd = 2)
}
h <- hist(test3$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test3$event))~test3$prob3,span=1))
lines_data <- data.frame(test3$prob3,obs_all)
lines_data2 <- lines_data[order(test3$prob3),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 4 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test3$prob4,breaks=quantile(test3$prob4, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test3,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob4))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[4],ylab="Observed",xlab="Expected", pch =16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[4], lwd = 2)
}
h <- hist(test3$prob4, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test3$event))~test3$prob4,span=1))
lines_data <- data.frame(test3$prob4,obs_all)
lines_data2 <- lines_data[order(test3$prob4),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[4],"black",cols[4],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

LT3<-c(1.5,1.5,1.5,1.5)
Mod<-c(1,2,3,4)
LT.3<-cbind(LT3, Mod, BS3, c3, c3.L, c3.U, Clge3, 
            Clge3.L, Clge3.U,Cslope3, Cslope3.L,
            Cslope3.U, fit_times, pred_times3)



## ## LT = 2

test4<-read.dta("validation_24binary.dta")
test4$trtment<-factor(test4$trtment)
test4$smoke<-factor(test4$smoke)
test4$interval<-rep(1, nrow(test4))
test4$timevent<-rep(3, nrow(test4))
test4$timevent1<-I(test4$timevent^3)
test4$timevent2<-I(test4$timevent^3*log(test4$timevent))

test4$one<-rep(1, nrow(test4))
test4$roundtime2<-rep(0, nrow(test4))
test4$roundtime3<-rep(1, nrow(test4))
test4$roundtime4<-rep(0, nrow(test4))

#Gee1 - No time since baseline covariate

lpcov14<-test4 %>% select(one, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat14<-as.matrix(lpcov14)
lp14<-gee1coef %*% t(lpmat14)
lp14<-as.vector(t(lp14))

#Gee2 - Categorical time since baseline covariate

lpcov24<-test4 %>% select(one, roundtime2, roundtime3, roundtime4, age, pgen, previous_dmards, 
                          disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat24<-as.matrix(lpcov24)
lp24<-gee2coef %*% t(lpmat24)
lp24<-as.vector(t(lp24))

#Gee3 - Linear continuous time since baseline

lpcov34<-test4 %>% select(one, timevent, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat34<-as.matrix(lpcov34)
lp34<-gee3coef %*% t(lpmat34)
lp34<-as.vector(t(lp34))


#Gee4 - Nonlinear continuous time since baseline

lpcov44<-test4 %>% select(one, timevent1, timevent2, age, pgen, previous_dmards, disdur, trt2, 
                          trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat44<-as.matrix(lpcov44)
lp44<-gee4coef %*% t(lpmat44)
lp44<-as.vector(t(lp44))


## Probability of at least one serious infection per year
p_start <- Sys.time()
test4$prob1<-1-exp(-1*exp(lp14))
p_end <- Sys.time()
p_start2 <- Sys.time()
test4$prob2<-1-exp(-1*exp(lp24))
p_end2 <- Sys.time()
p_start3 <- Sys.time()
test4$prob3<-1-exp(-1*exp(lp34))
p_end3 <- Sys.time()
p_start4 <- Sys.time()
test4$prob4<-1-exp(-1*exp(lp44))
p_end4 <- Sys.time()

p_time4 <- p_end - p_start
p2_time4 <- p2_end - p2_start
p3_time4 <- p3_end - p3_start
p4_time4 <- p4_end - p4_start

pred_times4 <- c(p_time4, p2_time4, p3_time4, p4_time4)

BS4<-0
BS4[1]<-BrierScore(resp= test4$event, pred = test4$prob1, 
                   scaled = FALSE)
BS4[2]<-BrierScore(resp= test4$event, pred = test4$prob2, 
                   scaled = FALSE)
BS4[3]<-BrierScore(resp= test4$event, pred = test4$prob3, 
                   scaled = FALSE)
BS4[4]<-BrierScore(resp= test4$event, pred = test4$prob4, 
                   scaled = FALSE)


# C statistic / AUC
c4<-0
c4.L<-0
c4.U<-0

C1 <- roc(test4$event~test4$prob1,ci=TRUE)
c4[1]<-as.numeric(C1$auc)
c4.L[1]<-as.numeric(C1$ci[1])
c4.U[1]<-as.numeric(C1$ci[3])
C2 <- roc(test4$event~test4$prob2,ci=TRUE)
c4[2]<-as.numeric(C2$auc)
c4.L[2]<-as.numeric(C2$ci[1])
c4.U[2]<-as.numeric(C2$ci[3])
C3 <- roc(test4$event~test4$prob3,ci=TRUE)
c4[3]<-as.numeric(C3$auc)
c4.L[3]<-as.numeric(C3$ci[1])
c4.U[3]<-as.numeric(C3$ci[3])
C4 <- roc(test4$event~test4$prob4,ci=TRUE)
c4[4]<-as.numeric(C4$auc)
c4.L[4]<-as.numeric(C4$ci[1])
c4.U[4]<-as.numeric(C4$ci[3])

#  Calibration

# Calibration In The Large
set.seed(112)

Clge4<-0
Clge4.L<-0
Clge4.U<-0

m1 <- glm(test4$event~offset(lp14), family="binomial")
Clge4[1] <- as.numeric(m1$coef)
x14<-confint(m1)
Clge4.L[1]<-as.numeric(x14[1])
Clge4.U[1]<-as.numeric(x14[2])

m2 <- glm(test4$event~offset(lp24), family="binomial")
Clge4[2] <- as.numeric(m2$coef)
x24<-confint(m2)
Clge4.L[2]<-as.numeric(x24[1])
Clge4.U[2]<-as.numeric(x24[2])

m3 <- glm(test4$event~offset(lp34), family="binomial")
Clge4[3] <- as.numeric(m3$coef)
x34<-confint(m3)
Clge4.L[3]<-as.numeric(x34[1])
Clge4.U[3]<-as.numeric(x34[2])

m4 <- glm(test4$event~offset(lp44), family="binomial")
Clge4[4] <- as.numeric(m4$coef)
x44<-confint(m4)
Clge4.L[4]<-as.numeric(x44[1])
Clge4.U[4]<-as.numeric(x44[2])


# Calibration slope
set.seed(182)

Cslope4<-0
Cslope4.L<-0
Cslope4.U<-0

m41 <- glm(test4$event~lp14,family="binomial", x=TRUE,y=TRUE)
Cslope4[1]<-m41$coefficients[2]
s14<-confint(m41)
Cslope4.L[1]<-as.numeric(s14[2,1])
Cslope4.U[1]<-as.numeric(s14[2,2])

m42 <- glm(test4$event~lp24,family="binomial", x=TRUE,y=TRUE)
Cslope4[2]<-m42$coefficients[2]
s24<-confint(m42)
Cslope4.L[2]<-as.numeric(s24[2,1])
Cslope4.U[2]<-as.numeric(s24[2,2])

m43 <- glm(test4$event~lp34,family="binomial", x=TRUE,y=TRUE)
Cslope4[3]<-m43$coefficients[2]
s34<-confint(m43)
Cslope4.L[3]<-as.numeric(s34[2,1])
Cslope4.U[3]<-as.numeric(s34[2,2])

m44 <- glm(test4$event~lp44,family="binomial", x=TRUE,y=TRUE)
Cslope4[4]<-m44$coefficients[2]
s44<-confint(m44)
Cslope4.L[4]<-as.numeric(s44[2,1])
Cslope4.U[4]<-as.numeric(s44[2,2])

par(mfrow=c(2,2))
par(pty="s")

# 1 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test4$prob1,breaks=quantile(test4$prob1, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test4,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob1))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1], lwd = 2)
}
h <- hist(test4$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test4$event))~test4$prob1,span=1))
lines_data <- data.frame(test4$prob1,obs_all)
lines_data2 <- lines_data[order(test4$prob1),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 2 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test4$prob2,breaks=quantile(test4$prob2, 
                                          prob = 
                                            c(0,0.1,0.2,0.3,0.4,0.5,
                                              0.6,0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test4,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob2))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[2],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[2], lwd = 2)
}
h <- hist(test4$prob2, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test4$event))~test4$prob2,span=1))
lines_data <- data.frame(test4$prob2,obs_all)
lines_data2 <- lines_data[order(test4$prob2),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),
       lty=c(0,2,1,1),pch=c(16,NA,NA,NA),bty="n")

# 3 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test4$prob3,breaks=quantile(test4$prob3, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,
                                                   0.4,0.5,0.6,
                                                   0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of 
#patients in each risk group 
gpdata <- cbind(test4,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob3))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[3],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[3], lwd = 2)
}
h <- hist(test4$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test4$event))~test4$prob3,span=1))
lines_data <- data.frame(test4$prob3,obs_all)
lines_data2 <- lines_data[order(test4$prob3),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 4 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test4$prob4,breaks=quantile(test4$prob4, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test4,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob4))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[4],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[4], lwd = 2)
}
h <- hist(test4$prob4, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test4$event))~test4$prob4,span=1))
lines_data <- data.frame(test4$prob4,obs_all)
lines_data2 <- lines_data[order(test4$prob4),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[4],"black",cols[4],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

LT4<-c(2,2,2,2)
Mod<-c(1,2,3,4)
LT.4<-cbind(LT4, Mod, BS4, c4, c4.L, c4.U, Clge4, 
            Clge4.L, Clge4.U,Cslope4, Cslope4.L,
            Cslope4.U, fit_times, pred_times4)



## ## LT = 2.5

test5<-read.dta("validation_30binary.dta")
test5$trtment<-factor(test5$trtment)
test5$smoke<-factor(test5$smoke)
test5$interval<-rep(1, nrow(test5))
test5$timevent<-rep(3.5, nrow(test5))
test5$timevent1<-I(test5$timevent^3)
test5$timevent2<-I(test5$timevent^3*log(test5$timevent))

test5$one<-rep(1, nrow(test5))
test5$roundtime2<-rep(0, nrow(test5))
test5$roundtime3<-rep(1, nrow(test5))
test5$roundtime4<-rep(0, nrow(test5))

#Gee1 - No time since baseline covariate

lpcov15<-test5 %>% select(one, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat15<-as.matrix(lpcov15)
lp15<-gee1coef %*% t(lpmat15)
lp15<-as.vector(t(lp15))

#Gee2 - Categorical time since baseline covariate

lpcov25<-test5 %>% select(one, roundtime2, roundtime3, roundtime4, age, pgen, previous_dmards, 
                          disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat25<-as.matrix(lpcov25)
lp25<-gee2coef %*% t(lpmat25)
lp25<-as.vector(t(lp25))

#Gee3 - Linear continuous time since baseline

lpcov35<-test5 %>% select(one, timevent, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat35<-as.matrix(lpcov35)
lp35<-gee3coef %*% t(lpmat35)
lp35<-as.vector(t(lp35))


#Gee4 - Nonlinear continuous time since baseline

lpcov45<-test5 %>% select(one, timevent1, timevent2, age, pgen, previous_dmards, disdur, trt2, 
                          trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat45<-as.matrix(lpcov45)
lp45<-gee4coef %*% t(lpmat45)
lp45<-as.vector(t(lp45))

## Probability of at least one serious infection per year
p_start <- Sys.time()
test5$prob1<-1-exp(-1*exp(lp15))
p_end <- Sys.time()
p_start2 <- Sys.time()
test5$prob2<-1-exp(-1*exp(lp25))
p_end2 <- Sys.time()
p_start3 <- Sys.time()
test5$prob3<-1-exp(-1*exp(lp35))
p_end3 <- Sys.time()
p_start4 <- Sys.time()
test5$prob4<-1-exp(-1*exp(lp45))
p_end4 <- Sys.time()

p_time5 <- p_end - p_start
p2_time5 <- p2_end - p2_start
p3_time5 <- p3_end - p3_start
p4_time5 <- p4_end - p4_start

pred_times5 <- c(p_time5, p2_time5, p3_time5, p4_time5)


BS5<-0
BS5[1]<-BrierScore(resp= test5$event, pred = test5$prob1, 
                   scaled = FALSE)
BS5[2]<-BrierScore(resp= test5$event, pred = test5$prob2, 
                   scaled = FALSE)
BS5[3]<-BrierScore(resp= test5$event, pred = test5$prob3, 
                   scaled = FALSE)
BS5[4]<-BrierScore(resp= test5$event, pred = test5$prob4, 
                   scaled = FALSE)

# C statistic / AUC
c5<-0
c5.L<-0
c5.U<-0

C1 <- roc(test5$event~test5$prob1,ci=TRUE)
c5[1]<-as.numeric(C1$auc)
c5.L[1]<-as.numeric(C1$ci[1])
c5.U[1]<-as.numeric(C1$ci[3])
C2 <- roc(test5$event~test5$prob2,ci=TRUE)
c5[2]<-as.numeric(C2$auc)
c5.L[2]<-as.numeric(C2$ci[1])
c5.U[2]<-as.numeric(C2$ci[3])
C3 <- roc(test5$event~test5$prob3,ci=TRUE)
c5[3]<-as.numeric(C3$auc)
c5.L[3]<-as.numeric(C3$ci[1])
c5.U[3]<-as.numeric(C3$ci[3])
C4 <- roc(test5$event~test5$prob4,ci=TRUE)
c5[4]<-as.numeric(C4$auc)
c5.L[4]<-as.numeric(C4$ci[1])
c5.U[4]<-as.numeric(C4$ci[3])


#  Calibration

# Calibration In The Large
set.seed(112)

Clge5<-0
Clge5.L<-0
Clge5.U<-0

m1 <- glm(test5$event~offset(lp15), family="binomial")
Clge5[1] <- as.numeric(m1$coef)
x15<-confint(m1)
Clge5.L[1]<-as.numeric(x15[1])
Clge5.U[1]<-as.numeric(x15[2])

m2 <- glm(test5$event~offset(lp25), family="binomial")
Clge5[2] <- as.numeric(m2$coef)
x25<-confint(m2)
Clge5.L[2]<-as.numeric(x25[1])
Clge5.U[2]<-as.numeric(x25[2])

m3 <- glm(test5$event~offset(lp35), family="binomial")
Clge5[3] <- as.numeric(m3$coef)
x35<-confint(m3)
Clge5.L[3]<-as.numeric(x35[1])
Clge5.U[3]<-as.numeric(x35[2])

m4 <- glm(test5$event~offset(lp45), family="binomial")
Clge5[4] <- as.numeric(m4$coef)
x45<-confint(m4)
Clge5.L[4]<-as.numeric(x45[1])
Clge5.U[4]<-as.numeric(x45[2])


# Calibration slope
set.seed(182)

Cslope5<-0
Cslope5.L<-0
Cslope5.U<-0

m51 <- glm(test5$event~lp15,family="binomial", x=TRUE,y=TRUE)
Cslope5[1]<-m51$coefficients[2]
s15<-confint(m51)
Cslope5.L[1]<-as.numeric(s15[2,1])
Cslope5.U[1]<-as.numeric(s15[2,2])

m52 <- glm(test5$event~lp25,family="binomial", x=TRUE,y=TRUE)
Cslope5[2]<-m52$coefficients[2]
s25<-confint(m52)
Cslope5.L[2]<-as.numeric(s25[2,1])
Cslope5.U[2]<-as.numeric(s25[2,2])

m53 <- glm(test5$event~lp35,family="binomial", x=TRUE,y=TRUE)
Cslope5[3]<-m53$coefficients[2]
s35<-confint(m53)
Cslope5.L[3]<-as.numeric(s35[2,1])
Cslope5.U[3]<-as.numeric(s35[2,2])

m54 <- glm(test5$event~lp45,family="binomial", x=TRUE,y=TRUE)
Cslope5[4]<-m54$coefficients[2]
s45<-confint(m54)
Cslope5.L[4]<-as.numeric(s45[2,1])
Cslope5.U[4]<-as.numeric(s45[2,2])

par(mfrow=c(2,2))
par(pty="s")

# 1 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test5$prob1,breaks=quantile(test5$prob1, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test5,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob1))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1], lwd = 2)
}
h <- hist(test5$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test5$event))~test5$prob1,span=1))
lines_data <- data.frame(test5$prob1,obs_all)
lines_data2 <- lines_data[order(test5$prob1),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 2 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test5$prob1,breaks=quantile(test5$prob2, 
                                          prob = 
                                            c(0,0.1,0.2,0.3,0.4,0.5,
                                              0.6,0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test5,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob2))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[2],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[2], lwd = 2)
}
h <- hist(test5$prob2, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test5$event))~test5$prob2,span=1))
lines_data <- data.frame(test5$prob2,obs_all)
lines_data2 <- lines_data[order(test5$prob2),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),
       lty=c(0,2,1,1),pch=c(16,NA,NA,NA),bty="n")

# 3 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test5$prob3,breaks=quantile(test5$prob3, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,
                                                   0.4,0.5,0.6,
                                                   0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of 
#patients in each risk group 
gpdata <- cbind(test5,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob3))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[3],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[3], lwd = 2)
}
h <- hist(test5$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test5$event))~test5$prob3,span=1))
lines_data <- data.frame(test5$prob3,obs_all)
lines_data2 <- lines_data[order(test5$prob3),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 4 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test5$prob4,breaks=quantile(test5$prob4, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test5,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob4))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[4],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[4], lwd = 2)
}
h <- hist(test5$prob4, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test5$event))~test5$prob4,span=1))
lines_data <- data.frame(test5$prob4,obs_all)
lines_data2 <- lines_data[order(test5$prob4),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[4],"black",cols[4],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

LT5<-c(2.5,2.5,2.5,2.5)
Mod<-c(1,2,3,4)
LT.5<-cbind(LT5, Mod, BS5, c5, c5.L, c5.U, Clge5, 
            Clge5.L, Clge5.U,Cslope5, Cslope5.L,
            Cslope5.U, fit_times, pred_times5)


## ## LT = 3

test6<-read.dta("validation_36binary.dta")
test6$trtment<-factor(test6$trtment)
test6$smoke<-factor(test6$smoke)
test6$interval<-rep(1, nrow(test6))
test6$timevent<-rep(4, nrow(test6))
test6$timevent1<-I(test6$timevent^3)
test6$timevent2<-I(test6$timevent^3*log(test6$timevent))

test6$one<-rep(1, nrow(test6))
test6$roundtime2<-rep(0, nrow(test6))
test6$roundtime3<-rep(0, nrow(test6))
test6$roundtime4<-rep(1, nrow(test6))

#Gee1 - No time since baseline covariate

lpcov16<-test6 %>% select(one, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat16<-as.matrix(lpcov16)
lp16<-gee1coef %*% t(lpmat16)
lp16<-as.vector(t(lp16))

#Gee2 - Categorical time since baseline covariate

lpcov26<-test6 %>% select(one, roundtime2, roundtime3, roundtime4, age, pgen, previous_dmards, 
                          disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat26<-as.matrix(lpcov26)
lp26<-gee2coef %*% t(lpmat26)
lp26<-as.vector(t(lp26))

#Gee3 - Linear continuous time since baseline

lpcov36<-test6 %>% select(one, timevent, age, pgen, previous_dmards, disdur, trt2, trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat36<-as.matrix(lpcov36)
lp36<-gee3coef %*% t(lpmat36)
lp36<-as.vector(t(lp36))


#Gee4 - Nonlinear continuous time since baseline

lpcov46<-test6 %>% select(one, timevent1, timevent2, age, pgen, previous_dmards, disdur, trt2, 
                          trt4, trt1561, trt10060,
                          ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat46<-as.matrix(lpcov46)
lp46<-gee4coef %*% t(lpmat46)
lp46<-as.vector(t(lp46))

## Probability of at least one serious infection per year
p_start <- Sys.time()
test6$prob1<-1-exp(-1*exp(lp16))
p_end <- Sys.time()
p_start2 <- Sys.time()
test6$prob2<-1-exp(-1*exp(lp26))
p_end2 <- Sys.time()
p_start3 <- Sys.time()
test6$prob3<-1-exp(-1*exp(lp36))
p_end3 <- Sys.time()
p_start4 <- Sys.time()
test6$prob4<-1-exp(-1*exp(lp46))
p_end4 <- Sys.time()

p_time6 <- p_end - p_start
p2_time6 <- p2_end - p2_start
p3_time6 <- p3_end - p3_start
p4_time6 <- p4_end - p4_start


BS6<-0
BS6[1]<-BrierScore(resp= test6$event, pred = test6$prob1, 
                   scaled = FALSE)
BS6[2]<-BrierScore(resp= test6$event, pred = test6$prob2, 
                   scaled = FALSE)
BS6[3]<-BrierScore(resp= test6$event, pred = test6$prob3, 
                   scaled = FALSE)
BS6[4]<-BrierScore(resp= test6$event, pred = test6$prob4, 
                   scaled = FALSE)

# C statistic / AUC
c6<-0
c6.L<-0
c6.U<-0

C1 <- roc(test6$event~test6$prob1,ci=TRUE)
c6[1]<-as.numeric(C1$auc)
c6.L[1]<-as.numeric(C1$ci[1])
c6.U[1]<-as.numeric(C1$ci[3])
C2 <- roc(test6event~test6$prob2,ci=TRUE)
c6[2]<-as.numeric(C2$auc)
c6.L[2]<-as.numeric(C2$ci[1])
c6.U[2]<-as.numeric(C2$ci[3])
C3 <- roc(test6$event~test6$prob3,ci=TRUE)
c6[3]<-as.numeric(C3$auc)
c6.L[3]<-as.numeric(C3$ci[1])
c6.U[3]<-as.numeric(C3$ci[3])
C4 <- roc(test6$event~test6$prob4,ci=TRUE)
c6[4]<-as.numeric(C4$auc)
c6.L[4]<-as.numeric(C4$ci[1])
c6.U[4]<-as.numeric(C4$ci[3])


#  Calibration

# Calibration In The Large
set.seed(112)

Clge6<-0
Clge6.L<-0
Clge6.U<-0

m1 <- glm(test6$event~offset(lp16), family="binomial")
Clge6[1] <- as.numeric(m1$coef)
x16<-confint(m1)
Clge6.L[1]<-as.numeric(x16[1])
Clge6.U[1]<-as.numeric(x16[2])

m2 <- glm(test6$event~offset(lp26), family="binomial")
Clge6[2] <- as.numeric(m2$coef)
x26<-confint(m2)
Clge6.L[2]<-as.numeric(x26[1])
Clge6.U[2]<-as.numeric(x26[2])

m3 <- glm(test6$event~offset(lp36), family="binomial")
Clge6[3] <- as.numeric(m3$coef)
x36<-confint(m3)
Clge6.L[3]<-as.numeric(x36[1])
Clge6.U[3]<-as.numeric(x36[2])

m4 <- glm(test6$event~offset(lp46), family="binomial")
Clge6[4] <- as.numeric(m4$coef)
x46<-confint(m4)
Clge6.L[4]<-as.numeric(x46[1])
Clge6.U[4]<-as.numeric(x46[2])


# Calibration slope
set.seed(182)

Cslope6<-0
Cslope6.L<-0
Cslope6.U<-0

m61 <- glm(test6$event~lp16,family="binomial", x=TRUE,y=TRUE)
Cslope6[1]<-m61$coefficients[2]
s16<-confint(m61)
Cslope6.L[1]<-as.numeric(s16[2,1])
Cslope6.U[1]<-as.numeric(s16[2,2])

m62 <- glm(test6$event~lp26,family="binomial", x=TRUE,y=TRUE)
Cslope6[2]<-m62$coefficients[2]
s26<-confint(m62)
Cslope6.L[2]<-as.numeric(s26[2,1])
Cslope6.U[2]<-as.numeric(s26[2,2])

m63 <- glm(test6$event~lp36,family="binomial", x=TRUE,y=TRUE)
Cslope6[3]<-m63$coefficients[2]
s36<-confint(m63)
Cslope6.L[3]<-as.numeric(s36[2,1])
Cslope6.U[3]<-as.numeric(s36[2,2])

m64 <- glm(test6$event~lp46,family="binomial", x=TRUE,y=TRUE)
Cslope6[4]<-m64$coefficients[2]
s46<-confint(m64)
Cslope6.L[4]<-as.numeric(s46[2,1])
Cslope6.U[4]<-as.numeric(s46[2,2])

par(mfrow=c(2,2))
par(pty="s")

# 1 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test6$prob1,breaks=quantile(test6$prob1, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test6,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob1))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[1], lwd = 2)
}
h <- hist(test6$prob1, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test6$event))~test6$prob1,span=1))
lines_data <- data.frame(test6$prob1,obs_all)
lines_data2 <- lines_data[order(test6$prob1),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 2 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test6$prob2,breaks=quantile(test6$prob2, 
                                          prob = 
                                            c(0,0.1,0.2,0.3,0.4,0.5,
                                              0.6,0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test6,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob2))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[2],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[2], lwd = 2)
}
h <- hist(test5$prob2, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test6$event))~test6$prob2,span=1))
lines_data <- data.frame(test6$prob2,obs_all)
lines_data2 <- lines_data[order(test6$prob2),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),
       lty=c(0,2,1,1),pch=c(16,NA,NA,NA),bty="n")

# 3 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test6$prob3,breaks=quantile(test6$prob3, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,
                                                   0.4,0.5,0.6,
                                                   0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of 
#patients in each risk group 
gpdata <- cbind(test6,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob3))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),col=cols[3],
     ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[3], lwd = 2)
}
h <- hist(test6$prob3, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test6$event))~test6$prob3,span=1))
lines_data <- data.frame(test6$prob3,obs_all)
lines_data2 <- lines_data[order(test6$prob3),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

# 4 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test6$prob4,breaks=quantile(test6$prob4, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test6,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob4))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.


plot(obs~exp[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[4],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col=cols[4], lwd = 2)
}
h <- hist(test6$prob4, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                 1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test6$event))~test6$prob4,span=1))
lines_data <- data.frame(test6$prob4,obs_all)
lines_data2 <- lines_data[order(test6$prob4),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend("topleft",c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[4],"black",cols[4],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

pred_times6<-c(p_time6, p2_time6, p3_time6, p4_time6)
LT6<-c(3,3,3,3)
Mod<-c(1,2,3,4)
LT.6<-cbind(LT6, Mod, BS6, c6, c6.L, c6.U, Clge6, 
            Clge6.L, Clge6.U,Cslope6, Cslope6.L,
            Cslope6.U, fit_times, pred_times6)

# Temporal assessment vectors

GEE.TA<-rbind(LT.0, LT.1, LT.2, LT.3, LT.4, LT.5, LT.6)
GEE.TA<-as.data.frame(GEE.TA)
write.table(GEE.TA, file = "geeTA.csv", sep = ",", col.names = NA,
            qmethod = "double")
x<-read.table("geeTA.csv", header = TRUE, sep = ",", row.names = 1)
head(x)

###### Temporal assessment survival edition
library(ipred)
##### LT = 0 #####
LT<-c(0,0.5,1,1.5,2,2.5,3)
Cstat<-0
CstatL<-0
CstatU<-0
BS<-0

# Obtain the predicted probabilities for each subject
test0<-read.dta("baseline_model.dta")
test0$smoke <- as.factor(test0$smoke)
test0$firsttreat <- as.factor(test0$firsttreat)
test0$one<-rep(1, times = nrow(test0))
geecoef <-as.numeric(gee1$coefficients)
lpcov0<-test0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                         trt4, trt1561, trt10060,
                         ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat0<-as.matrix(lpcov0)
pred_LP0<-geecoef %*% t(lpmat0)
test0$pred_LP0<-as.vector(t(pred_LP0))
test0$surv_pred0<-exp(-1*exp(test0$pred_LP0))

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
test0$smoke <- as.factor(test0$smoke)
test0$firsttreat <- as.factor(test0$firsttreat)
test0$one<-rep(1, times = nrow(test0))
geecoef <-as.numeric(gee1$coefficients)
lpcov0<-test0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                         trt4, trt1561, trt10060,
                         ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat0<-as.matrix(lpcov0)
pred_LP0<-geecoef %*% t(lpmat0)
test0$pred_LP0<-as.vector(t(pred_LP0))
test0$surv_pred0<-exp(-1*exp(test0$pred_LP0))

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
              test0$surv_pred0, btime = 1.5)


##### LT = 1 #####

# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_12months.dta")
test0$smoke <- as.factor(test0$smoke)
test0$one<-rep(1, times = nrow(test0))
test0$firsttreat <- as.factor(test0$firsttreat)
geecoef <-as.numeric(gee1$coefficients)
lpcov0<-test0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                         trt4, trt1561, trt10060,
                         ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat0<-as.matrix(lpcov0)
pred_LP0<-geecoef %*% t(lpmat0)
test0$pred_LP0<-as.vector(t(pred_LP0))
test0$surv_pred0<-exp(-1*exp(test0$pred_LP0))

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
              test0$surv_pred0, btime = 2)



##### LT = 1.5 #####

# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_18months.dta")
test0$smoke <- as.factor(test0$smoke)
test0$firsttreat <- as.factor(test0$firsttreat)
test0$one<-rep(1, times = nrow(test0))
geecoef <-as.numeric(gee1$coefficients)
lpcov0<-test0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                         trt4, trt1561, trt10060,
                         ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat0<-as.matrix(lpcov0)
pred_LP0<-geecoef %*% t(lpmat0)
test0$pred_LP0<-as.vector(t(pred_LP0))
test0$surv_pred0<-exp(-1*exp(test0$pred_LP0))

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
              test0$surv_pred0, btime = 2.5)

##### LT = 2 #####


# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_24months.dta")
test0$smoke <- as.factor(test0$smoke)
test0$firsttreat <- as.factor(test0$firsttreat)
test0$one<-rep(1, times = nrow(test0))
geecoef <-as.numeric(gee1$coefficients)
lpcov0<-test0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                         trt4, trt1561, trt10060,
                         ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat0<-as.matrix(lpcov0)
pred_LP0<-geecoef %*% t(lpmat0)
test0$pred_LP0<-as.vector(t(pred_LP0))
test0$surv_pred0<-exp(-1*exp(test0$pred_LP0))

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
              test0$surv_pred0, btime = 3)


##### LT = 2.5 #####


# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_30months.dta")
test0$smoke <- as.factor(test0$smoke)
test0$firsttreat <- as.factor(test0$firsttreat)
test0$one<-rep(1, times = nrow(test0))
geecoef <-as.numeric(gee1$coefficients)
lpcov0<-test0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                         trt4, trt1561, trt10060,
                         ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat0<-as.matrix(lpcov0)
pred_LP0<-geecoef %*% t(lpmat0)
test0$pred_LP0<-as.vector(t(pred_LP0))
test0$surv_pred0<-exp(-1*exp(test0$pred_LP0))

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
              test0$surv_pred0, btime = 3.5)


##### LT = 3 #####


# Obtain the predicted probabilities for each subject
test0<-read.dta("validation_36months.dta")
test0$smoke <- as.factor(test0$smoke)
test0$firsttreat <- as.factor(test0$firsttreat)
test0$one<-rep(1, times = nrow(test0))
geecoef <-as.numeric(gee1$coefficients)
lpcov0<-test0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                         trt4, trt1561, trt10060,
                         ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
lpmat0<-as.matrix(lpcov0)
pred_LP0<-geecoef %*% t(lpmat0)
test0$pred_LP0<-as.vector(t(pred_LP0))
test0$surv_pred0<-exp(-1*exp(test0$pred_LP0))

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
              test0$surv_pred0, btime = 4)
# Performance data frame


GEE.TA_surv<-cbind(LT, Cstat, CstatL, CstatU, BS)
GEE.TA_surv<-as.data.frame(GEE.TA_surv)
write.table(GEE.TA_surv, file = "geeTA_surv.csv", sep = ",", col.names = NA,
            qmethod = "double")

## Survival bootstrap function

bootvadGEE_surv<-function(data, B, seed){
  dat <- data
  dat$studyno <- dat$groupid
  
  indiv <- unique(dat$studyno)
  C<-0
  od.C<-0
  c<-matrix(0, nrow = B, ncol = 7)
  bs<-matrix(0, nrow = B, ncol = 7)
  od.c<-matrix(0, nrow = B, ncol = 7)
  bs.od<-matrix(0, nrow = B, ncol = 7)
  opt.d<-matrix(0, nrow = B, ncol = 7)
  opt.c<-matrix(0, nrow = B, ncol = 7)
  results<-matrix(0, nrow = 2, ncol = 7)
  
  # original data - temporal assessment
  
  od0 <- read.dta("baseline_model.dta")
  od0$trtment<-as.factor(od0$trtment)
  od0$smoke<-as.factor(od0$smoke)
  od1 <- read.dta("validation_6months.dta")
  od1$trtment<-as.factor(od1$trtment)
  od1$smoke<-as.factor(od1$smoke)
  od2 <- read.dta("validation_12months.dta")
  od2$trtment<-as.factor(od2$trtment)
  od2$smoke<-as.factor(od2$smoke)
  od3 <- read.dta("validation_18months.dta")
  od3$trtment<-as.factor(od3$trtment)
  od3$smoke<-as.factor(od3$smoke)
  od4 <- read.dta("validation_24months.dta")
  od4$trtment<-as.factor(od4$trtment)
  od4$smoke<-as.factor(od4$smoke)
  od5 <- read.dta("validation_30months.dta")
  od5$trtment<-as.factor(od5$trtment)
  od5$smoke<-as.factor(od5$smoke)
  od6 <- read.dta("validation_36months.dta")
  od6$trtment<-as.factor(od6$trtment)
  od6$smoke<-as.factor(od6$smoke)
  
  od0$studyno<-od0$groupid
  od1$studyno<-od1$groupid
  od2$studyno<-od2$groupid
  od3$studyno<-od3$groupid
  od4$studyno<-od4$groupid
  od5$studyno<-od5$groupid
  od6$studyno<-od6$groupid
  
  od0$one<-rep(1, times=nrow(od0))
  od1$one<-rep(1, times=nrow(od1))
  od2$one<-rep(1, times=nrow(od2))
  od3$one<-rep(1, times=nrow(od3))
  od4$one<-rep(1, times=nrow(od4))
  od5$one<-rep(1, times=nrow(od5))
  od6$one<-rep(1, times=nrow(od6))
  
  for (j in 1:B){
    seed<-seed+j
    set.seed(seed)
    smp <- sort(sample(indiv, length(indiv), replace=TRUE))
    smp.df <- data.frame(studyno=smp)
    print(length(smp.df$studyno))
    
    # renaming duplicated ID numbers
    
    smp.df<-smp.df %>% dplyr::group_by(studyno) %>% dplyr::mutate(count = dplyr::row_number())
    smp.df$newid<-paste(smp.df$studyno, smp.df$count, sep=".")
    dat.b = merge(smp.df, dat, by = "studyno", all.x=TRUE)
    dat.b$groupid <- dat.b$newid
    dat.b<-dat.b[with(dat.b, order(groupid, start)),]
    print(length(dat.b$studyno))
    
    # fitting GEE (nonlinear time since baseline)
    
    bgee4<-geeglm(serinf~ age + pgen + previous_dmards + 
                    disdur + trtment + ovmean + dascore + steroids +
                    bmi + renal + lung + diabetes + smoke + offset(log(interval)), 
                  family=binomial(link = "cloglog"), data=dat.b, 
                  id=groupid, waves=pyears, corstr = "exchangeable", 
                  std.err="san.se")
    
    # Temporal assessment sets
    
    val0 <- na.omit(merge(smp.df, od0, by = "studyno", all.x=TRUE))
    val1 <- na.omit(merge(smp.df, od1, by = "studyno", all.x=TRUE))
    val2 <- na.omit(merge(smp.df, od2, by = "studyno", all.x=TRUE))
    val3 <- na.omit(merge(smp.df, od3, by = "studyno", all.x=TRUE))
    val4 <- na.omit(merge(smp.df, od4, by = "studyno", all.x=TRUE))
    val5 <- na.omit(merge(smp.df, od5, by = "studyno", all.x=TRUE))
    val6 <- na.omit(merge(smp.df, od6, by = "studyno", all.x=TRUE))
    
    
    # linear predictors
    
    bGeecoef<-as.numeric(bgee4$coefficients)
    
    one0<-rep(1, nrow(val0))
    val0sub<- val0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                              trt4, trt1561, trt10060,
                              ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat0<-as.matrix(val0sub)
    lp0<-bGeecoef %*% t(lpmat0)
    lp0<-as.vector(t(lp0))
    
    one1<-rep(1, nrow(val1))
    val1sub<-val1 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat1<-as.matrix(val1sub)
    lp1<-bGeecoef %*% t(lpmat1)
    lp1<-as.vector(t(lp1))
    
    one2<-rep(1, nrow(val2))
    val2sub<-val2 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat2<-as.matrix(val2sub)
    lp2<-bGeecoef %*% t(lpmat2)
    lp2<-as.vector(t(lp2))
    
    one3<-rep(1, nrow(val3))
    val3sub<-val3 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat3<-as.matrix(val3sub)
    lp3<-bGeecoef %*% t(lpmat3)
    lp3<-as.vector(t(lp3))
    
    one4<-rep(1, nrow(val4))
    val4sub<-val4 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat4<-as.matrix(val4sub)
    lp4<-bGeecoef %*% t(lpmat4)
    lp4<-as.vector(t(lp4))
    
    one5<-rep(1, nrow(val5))
    val5sub<-val5 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat5<-as.matrix(val5sub)
    lp5<-bGeecoef %*% t(lpmat5)
    lp5<-as.vector(t(lp5))
    
    one6<-rep(1, nrow(val6))
    val6sub<- val6 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                              trt4, trt1561, trt10060,
                              ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat6<-as.matrix(val6sub)
    lp6<-bGeecoef %*% t(lpmat6)
    lp6<-as.vector(t(lp6))
    
    
    # predict at least one SI in 12 months
    
    val0$prob<-1-exp(-1*exp(lp0))
    val1$prob<-1-exp(-1*exp(lp1))
    val2$prob<-1-exp(-1*exp(lp2))
    val3$prob<-1-exp(-1*exp(lp3))
    val4$prob<-1-exp(-1*exp(lp4))
    val5$prob<-1-exp(-1*exp(lp5))
    val6$prob<-1-exp(-1*exp(lp6))
    
    val0$survprob <- 1 - val0$prob
    val1$survprob <- 1 - val1$prob
    val2$survprob <- 1 - val2$prob
    val3$survprob <- 1 - val3$prob
    val4$survprob <- 1 - val4$prob
    val5$survprob <- 1 - val5$prob
    val6$survprob <- 1 - val6$prob
    
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
    
    BS0<-sbrier(obj=Surv(val0$timevent, val0$event), val0$survprob, 
                btime=1)
    BS1<-sbrier(obj=Surv(val1$timevent, val1$event), val1$survprob, 
                btime=1.5)
    BS2<-sbrier(obj=Surv(val2$timevent, val2$event), val2$survprob, 
                btime=2)
    BS3<-sbrier(obj=Surv(val3$timevent, val3$event), val3$survprob, 
                btime=2.5)
    BS4<-sbrier(obj=Surv(val4$timevent, val4$event), val4$survprob, 
                btime=3)
    BS5<-sbrier(obj=Surv(val5$timevent, val5$event), val5$survprob, 
                btime=3.5)
    BS6<-sbrier(obj=Surv(val6$timevent, val6$event), val6$survprob, 
                btime=4)
    
    bs[j,]<-c(BS0, BS1, BS2, BS3, BS4, BS5, BS6)
    
    # Survival predictions - original data
    
    one0od<-rep(1, nrow(od0))
    od0sub<-od0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat0od<-as.matrix(od0sub)
    lp0od<-bGeecoef %*% t(lpmat0od)
    lp0od<-as.vector(t(lp0od))
    
    one1od<-rep(1, nrow(od1))
    od1sub<-od1 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat1od<-as.matrix(od1sub)
    lp1od<-bGeecoef %*% t(lpmat1od)
    lp1od<-as.vector(t(lp1od))
    
    one2od<-rep(1, nrow(od2))
    od2sub<-od2 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat2od<-as.matrix(od2sub)
    lp2od<-bGeecoef %*% t(lpmat2od)
    lp2od<-as.vector(t(lp2od))
    
    one3od<-rep(1, nrow(od3))
    od3sub<-od3 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat3od<-as.matrix(od3sub)
    lp3od<-bGeecoef %*% t(lpmat3od)
    lp3od<-as.vector(t(lp3od))
    
    one4od<-rep(1, nrow(od4))
    od4sub<-od4 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat4od<-as.matrix(od4sub)
    lp4od<-bGeecoef %*% t(lpmat4od)
    lp4od<-as.vector(t(lp4od))
    
    one5od<-rep(1, nrow(od5))
    od5sub<-od5 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat5od<-as.matrix(od5sub)
    lp5od<-bGeecoef %*% t(lpmat5od)
    lp5od<-as.vector(t(lp5od))
    
    one6od<-rep(1, nrow(od6))
    od6sub<-od6 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat6od<-as.matrix(od6sub)
    lp6od<-bGeecoef %*% t(lpmat6od)
    lp6od<-as.vector(t(lp6od))
    
    
    # predict at least one SI in 12 months - original data
    
    od0$prob<-1-exp(-1*exp(lp0od))
    od1$prob<-1-exp(-1*exp(lp1od))
    od2$prob<-1-exp(-1*exp(lp2od))
    od3$prob<-1-exp(-1*exp(lp3od))
    od4$prob<-1-exp(-1*exp(lp4od))
    od5$prob<-1-exp(-1*exp(lp5od))
    od6$prob<-1-exp(-1*exp(lp6od))
    
    od0$survprob<-1-od0$prob
    od1$survprob<-1-od1$prob
    od2$survprob<-1-od2$prob
    od3$survprob<-1-od3$prob
    od4$survprob<-1-od4$prob
    od5$survprob<-1-od5$prob
    od6$survprob<-1-od6$prob
    
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
    
    BS0od<-sbrier(obj=Surv(od0$timevent, od0$event), od0$survprob, 
                  btime=1)
    BS1od<-sbrier(obj=Surv(od1$timevent, od1$event), od1$survprob, 
                  btime=1.5)
    BS2od<-sbrier(obj=Surv(od2$timevent, od2$event), od2$survprob, 
                  btime=2)
    BS3od<-sbrier(obj=Surv(od3$timevent, od3$event), od3$survprob, 
                  btime=2.5)
    BS4od<-sbrier(obj=Surv(od4$timevent, od4$event), od4$survprob, 
                  btime=3)
    BS5od<-sbrier(obj=Surv(od5$timevent, od5$event), od5$survprob, 
                  btime=3.5)
    BS6od<-sbrier(obj=Surv(od6$timevent, od6$event), od6$survprob, 
                  btime=4)
    
    
    bs.od[j,]<-c(BS0od, BS1od, BS2od, BS3od, BS4od, BS5od, BS6od)
    print(bs.od)
    print(dim(c))
    print(dim(od.c))
    opt.d[j,]<-as.numeric(c[j,]-od.c[j,])
    opt.c[j,]<-as.numeric(bs.od[j,]-bs[j,])
    print(opt.d)
    print(opt.c)
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


# Bootstrap function

bootvadGEE<-function(data, B, seed){
  dat <- data
  dat$studyno <- dat$groupid
  
  indiv <- unique(dat$studyno)
  C<-0
  od.C<-0
  c<-matrix(0, nrow = B, ncol = 7)
  bs<-matrix(0, nrow = B, ncol = 7)
  od.c<-matrix(0, nrow = B, ncol = 7)
  bs.od<-matrix(0, nrow = B, ncol = 7)
  opt.d<-matrix(0, nrow = B, ncol = 7)
  opt.c<-matrix(0, nrow = B, ncol = 7)
  results<-matrix(0, nrow = 2, ncol = 7)
  
  # original data - temporal assessment
  
  od0 <- test0
  od1 <- test1
  od2 <- test2
  od3 <- test3
  od4 <- test4
  od5 <- test5
  od6 <- test6
  
  od0$studyno<-od0$groupid
  od1$studyno<-od1$groupid
  od2$studyno<-od2$groupid
  od3$studyno<-od3$groupid
  od4$studyno<-od4$groupid
  od5$studyno<-od5$groupid
  od6$studyno<-od6$groupid
  
  od0$one<-rep(1, times=nrow(od0))
  od1$one<-rep(1, times=nrow(od1))
  od2$one<-rep(1, times=nrow(od2))
  od3$one<-rep(1, times=nrow(od3))
  od4$one<-rep(1, times=nrow(od4))
  od5$one<-rep(1, times=nrow(od5))
  od6$one<-rep(1, times=nrow(od6))
  
  for (j in 1:B){
    seed<-seed+j
    set.seed(seed)
    smp <- sort(sample(indiv, length(indiv), replace=TRUE))
    smp.df <- data.frame(studyno=smp)
    print(length(smp.df$studyno))
    
    # renaming duplicated ID numbers
    
    smp.df<-smp.df %>% dplyr::group_by(studyno) %>% dplyr::mutate(count = dplyr::row_number())
    smp.df$newid<-paste(smp.df$studyno, smp.df$count, sep=".")
    dat.b = merge(smp.df, dat, by = "studyno", all.x=TRUE)
    dat.b$groupid <- dat.b$newid
    dat.b<-dat.b[with(dat.b, order(groupid, start)),]
    print(length(dat.b$studyno))
    
    # fitting GEE (nonlinear time since baseline)
    
    bgee4<-geeglm(serinf~age + pgen + previous_dmards + 
                   disdur + trtment + ovmean + dascore + steroids +
                   bmi + renal + lung + diabetes + smoke + offset(log(interval)), 
                 family=binomial(link = "cloglog"), data=dat.b, 
                 id=groupid, waves=pyears, corstr = "exchangeable", 
                 std.err="san.se")
    
    # Temporal assessment sets

    val0 <- na.omit(merge(smp.df, od0, by = "studyno", all.x=TRUE))
    val1 <- na.omit(merge(smp.df, od1, by = "studyno", all.x=TRUE))
    val2 <- na.omit(merge(smp.df, od2, by = "studyno", all.x=TRUE))
    val3 <- na.omit(merge(smp.df, od3, by = "studyno", all.x=TRUE))
    val4 <- na.omit(merge(smp.df, od4, by = "studyno", all.x=TRUE))
    val5 <- na.omit(merge(smp.df, od5, by = "studyno", all.x=TRUE))
    val6 <- na.omit(merge(smp.df, od6, by = "studyno", all.x=TRUE))
    
    
    # linear predictors
    
    bGeecoef<-as.numeric(bgee4$coefficients)
    
    one0<-rep(1, nrow(val0))
    val0sub<- val0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                         trt4, trt1561, trt10060,
                         ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat0<-as.matrix(val0sub)
    lp0<-bGeecoef %*% t(lpmat0)
    lp0<-as.vector(t(lp0))
    
    one1<-rep(1, nrow(val1))
    val1sub<-val1 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat1<-as.matrix(val1sub)
    lp1<-bGeecoef %*% t(lpmat1)
    lp1<-as.vector(t(lp1))
    
    one2<-rep(1, nrow(val2))
    val2sub<-val2 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat2<-as.matrix(val2sub)
    lp2<-bGeecoef %*% t(lpmat2)
    lp2<-as.vector(t(lp2))
    
    one3<-rep(1, nrow(val3))
    val3sub<-val3 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat3<-as.matrix(val3sub)
    lp3<-bGeecoef %*% t(lpmat3)
    lp3<-as.vector(t(lp3))
    
    one4<-rep(1, nrow(val4))
    val4sub<-val4 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat4<-as.matrix(val4sub)
    lp4<-bGeecoef %*% t(lpmat4)
    lp4<-as.vector(t(lp4))
    
    one5<-rep(1, nrow(val5))
    val5sub<-val5 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat5<-as.matrix(val5sub)
    lp5<-bGeecoef %*% t(lpmat5)
    lp5<-as.vector(t(lp5))
    
    one6<-rep(1, nrow(val6))
    val6sub<- val6 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                             trt4, trt1561, trt10060,
                             ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat6<-as.matrix(val6sub)
    lp6<-bGeecoef %*% t(lpmat6)
    lp6<-as.vector(t(lp6))
    

    # predict at least one SI in 12 months
    
    val0$prob<-1-exp(-1*exp(lp0))
    val1$prob<-1-exp(-1*exp(lp1))
    val2$prob<-1-exp(-1*exp(lp2))
    val3$prob<-1-exp(-1*exp(lp3))
    val4$prob<-1-exp(-1*exp(lp4))
    val5$prob<-1-exp(-1*exp(lp5))
    val6$prob<-1-exp(-1*exp(lp6))
    
    # Discrimination measure
    
    C0 <- roc(val0$event~val0$prob,ci=FALSE)
    C[1]<- as.numeric(C0$auc)
    C1 <- roc(val1$event~val1$prob,ci=FALSE)
    C[2]<- as.numeric(C1$auc)
    C2 <- roc(val2$event~val2$prob,ci=FALSE)
    C[3]<- as.numeric(C2$auc)
    C3 <- roc(val3$event~val3$prob,ci=FALSE)
    C[4]<- as.numeric(C3$auc)
    C4 <- roc(val4$event~val4$prob,ci=FALSE)
    C[5]<- as.numeric(C4$auc)
    C5 <- roc(val5$event~val5$prob,ci=FALSE)
    C[6]<- as.numeric(C5$auc)
    C6 <- roc(val6$event~val6$prob,ci=FALSE)
    C[7]<- as.numeric(C6$auc)
    
    c[j,]<-C
    
    # Calibration measure
    
    bs0 <-BrierScore(resp= val0$event, pred = val0$prob, 
                    scaled = FALSE)
    bs1 <-BrierScore(resp= val1$event, pred = val1$prob, 
                    scaled = FALSE)
    bs2 <-BrierScore(resp= val2$event, pred = val2$prob, 
                    scaled = FALSE)
    bs3 <-BrierScore(resp= val3$event, pred = val3$prob, 
                    scaled = FALSE)
    bs4 <-BrierScore(resp= val4$event, pred = val4$prob, 
                    scaled = FALSE)
    bs5 <-BrierScore(resp= val5$event, pred = val5$prob, 
                    scaled = FALSE)
    bs6 <-BrierScore(resp= val6$event, pred = val6$prob, 
                    scaled = FALSE)
    
    bs[j,]<-c(bs0, bs1, bs2, bs3, bs4, bs5, bs6)

    # Survival predictions - original data
    
    one0od<-rep(1, nrow(od0))
    od0sub<-od0 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                       trt4, trt1561, trt10060,
                       ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat0od<-as.matrix(od0sub)
    lp0od<-bGeecoef %*% t(lpmat0od)
    lp0od<-as.vector(t(lp0od))
    
    one1od<-rep(1, nrow(od1))
    od1sub<-od1 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat1od<-as.matrix(od1sub)
    lp1od<-bGeecoef %*% t(lpmat1od)
    lp1od<-as.vector(t(lp1od))
    
    one2od<-rep(1, nrow(od2))
    od2sub<-od2 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat2od<-as.matrix(od2sub)
    lp2od<-bGeecoef %*% t(lpmat2od)
    lp2od<-as.vector(t(lp2od))
    
    one3od<-rep(1, nrow(od3))
    od3sub<-od3 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat3od<-as.matrix(od3sub)
    lp3od<-bGeecoef %*% t(lpmat3od)
    lp3od<-as.vector(t(lp3od))
    
    one4od<-rep(1, nrow(od4))
    od4sub<-od4 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat4od<-as.matrix(od4sub)
    lp4od<-bGeecoef %*% t(lpmat4od)
    lp4od<-as.vector(t(lp4od))
    
    one5od<-rep(1, nrow(od5))
    od5sub<-od5 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat5od<-as.matrix(od5sub)
    lp5od<-bGeecoef %*% t(lpmat5od)
    lp5od<-as.vector(t(lp5od))
    
    one6od<-rep(1, nrow(od6))
    od6sub<-od6 %>% select(one, age, pgen, previous_dmards, disdur, trt2, 
                           trt4, trt1561, trt10060,
                           ovmean, dascore, steroids, bmi, renal, lung, diabetes, smoke1, smoke2)
    lpmat6od<-as.matrix(od6sub)
    lp6od<-bGeecoef %*% t(lpmat6od)
    lp6od<-as.vector(t(lp6od))
    
    
    # predict at least one SI in 12 months - original data
    
    od0$prob<-1-exp(-1*exp(lp0od))
    od1$prob<-1-exp(-1*exp(lp1od))
    od2$prob<-1-exp(-1*exp(lp2od))
    od3$prob<-1-exp(-1*exp(lp3od))
    od4$prob<-1-exp(-1*exp(lp4od))
    od5$prob<-1-exp(-1*exp(lp5od))
    od6$prob<-1-exp(-1*exp(lp6od))
    
    # Discrimination measure
    
    C0od <- roc(od0$event~od0$prob,ci=FALSE)
    od.C[1]<- as.numeric(C0od$auc)
    C1od <- roc(od1$event~od1$prob,ci=FALSE)
    od.C[2]<- as.numeric(C1od$auc)
    C2od <- roc(od2$event~od2$prob,ci=FALSE)
    od.C[3]<- as.numeric(C2od$auc)
    C3od <- roc(od3$event~od3$prob,ci=FALSE)
    od.C[4]<- as.numeric(C3od$auc)
    C4od <- roc(od4$event~od4$prob,ci=FALSE)
    od.C[5]<- as.numeric(C4od$auc)
    C5od <- roc(od5$event~od5$prob,ci=FALSE)
    od.C[6]<- as.numeric(C5od$auc)
    C6od <- roc(od6$event~od6$prob,ci=FALSE)
    od.C[7]<- as.numeric(C6od$auc)
    
    od.c[j,]<-od.C

    
    # Calibration measure
    
    bs0.od<-BrierScore(resp= od0$event, pred = od0$prob, 
                    scaled = FALSE)
    bs1.od<-BrierScore(resp= od1$event, pred = od1$prob, 
                    scaled = FALSE)
    bs2.od<-BrierScore(resp= od2$event, pred = od2$prob, 
                    scaled = FALSE)
    bs3.od<-BrierScore(resp= od3$event, pred = od3$prob, 
                    scaled = FALSE)
    bs4.od<-BrierScore(resp= od4$event, pred = od4$prob, 
                    scaled = FALSE)
    bs5.od<-BrierScore(resp= od5$event, pred = od5$prob, 
                    scaled = FALSE)
    bs6.od<-BrierScore(resp= od6$event, pred = od6$prob, 
                    scaled = FALSE)
    
    bs.od[j,]<-c(bs0.od, bs1.od, bs2.od, bs3.od, bs4.od, bs5.od, bs6.od)
    print(bs.od)
    print(dim(c))
    print(dim(od.c))
    opt.d[j,]<-as.numeric(c[j,]-od.c[j,])
    opt.c[j,]<-as.numeric(bs.od[j,]-bs[j,])
    print(opt.d)
    print(opt.c)
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

#trial <- bootvadGEE(Gdata, B=1, seed=999)
############################ parallelisation ####################################

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1) #always useful to keep one core free
start_time<-Sys.time()
registerDoParallel(cl)
GEE.opt <- foreach(i=1:200, .combine=rbind, .packages = c("geepack", "dplyr", 
                                                     "pROC", "DescTools")) %dopar% {
                                                       seed<-199429+i
                                                       bootvadGEE(Gdata, 1, seed)
                                                     }
stopCluster(cl)
end_time<-Sys.time()
end_time - start_time


GEE.opt.C=GEE.opt[seq(1,399,2),]
GEE.opt.B=GEE.opt[seq(2,400,2),]
GEE.opt2<-matrix(0, nrow=2, ncol=7)
GEE.opt2[1,]=c(mean(GEE.opt.C[,1]), mean(GEE.opt.C[,2]),
                mean(GEE.opt.C[,3]), mean(GEE.opt.C[,4]),
                mean(GEE.opt.C[,5]), mean(GEE.opt.C[,6]),
                mean(GEE.opt.C[,7]))
GEE.opt2[2,]=c(mean(GEE.opt.B[,1]), mean(GEE.opt.B[,2]),
                mean(GEE.opt.B[,3]), mean(GEE.opt.B[,4]),
                mean(GEE.opt.B[,5]), mean(GEE.opt.B[,6]),
                mean(GEE.opt.B[,7]))

GEE.OPT<-data.frame(LT0 = GEE.opt2[,1], LT0.5 = GEE.opt2[,2], LT1 = GEE.opt2[,3], 
                     LT1.5 = GEE.opt2[,4], LT2 = GEE.opt2[,5], LT2.5 = GEE.opt2[,6], 
                     LT3 = GEE.opt2[,7], row.names = c("C statistic (OE)", 
                                                        "Brier Score (OE)"))

write.table(GEE.OPT, file = "geeOPT.csv", sep = ",", col.names = NA,
            qmethod = "double")



########## survival bootstrap ################

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1) 
start_time<-Sys.time()
registerDoParallel(cl)
GEE.opt <- foreach(i=1:200, .combine=rbind, .packages = c("geepack", "dplyr", 
                                                          "pROC", "DescTools","foreign", "ipred", "survival")) %dopar% {
                                                            setwd("R:/BSRBR/Analyses/lucy_bull/Serious_infection/R Datasets")
                                                            seed<-199429+i
                                                            bootvadGEE_surv(Gdata, 1, seed)
                                                          }
stopCluster(cl)
end_time<-Sys.time()
end_time - start_time


GEE.opt.C=GEE.opt[seq(1,399,2),]
GEE.opt.B=GEE.opt[seq(2,400,2),]
GEE.opt2<-matrix(0, nrow=2, ncol=7)
GEE.opt2[1,]=c(mean(GEE.opt.C[,1]), mean(GEE.opt.C[,2]),
               mean(GEE.opt.C[,3]), mean(GEE.opt.C[,4]),
               mean(GEE.opt.C[,5]), mean(GEE.opt.C[,6]),
               mean(GEE.opt.C[,7]))
GEE.opt2[2,]=c(mean(GEE.opt.B[,1]), mean(GEE.opt.B[,2]),
               mean(GEE.opt.B[,3]), mean(GEE.opt.B[,4]),
               mean(GEE.opt.B[,5]), mean(GEE.opt.B[,6]),
               mean(GEE.opt.B[,7]))

GEE.OPT<-data.frame(LT0 = GEE.opt2[,1], LT0.5 = GEE.opt2[,2], LT1 = GEE.opt2[,3], 
                    LT1.5 = GEE.opt2[,4], LT2 = GEE.opt2[,5], LT2.5 = GEE.opt2[,6], 
                    LT3 = GEE.opt2[,7], row.names = c("C statistic (OE)", 
                                                      "Brier Score (OE)"))

write.table(GEE.OPT, file = "geeOPT_surv.csv", sep = ",", col.names = NA,
            qmethod = "double")

