library(foreign)
library(rstpm2)
library(geepack)
library(ROCR)
library(pROC)
library(plyr)

# data structures

pdata<-read.dta("poisson_model.dta")

pdata$smoke<-factor(pdata$smoke, levels=c("0","1","2"), ordered=F)
pdata$firsttreat<-factor(pdata$firsttreat, ordered = F)
pdata$roundtime<-as.factor(pdata$roundtime)

fdata <- read.dta("baseline_model2.dta")
fdata$smoke<-factor(fdata$smoke, levels=c("0","1","2"), ordered=F)
fdata$firsttreat<-factor(fdata$firsttreat, ordered = F)

cdata<-read.dta("cloglog_model.dta")
cdata$smoke<-factor(cdata$smoke, levels=c("0","1","2"), ordered=F)
cdata$firsttreat<-factor(cdata$firsttreat, ordered = F)
cdata$roundtime<-as.factor(cdata$roundtime)

lrdata<-read.dta("baseline_binary.dta")
lrdata$smoke<-factor(lrdata$smoke, levels=c("0","1","2"), ordered=F)
lrdata$firsttreat<-factor(lrdata$firsttreat, ordered = F)

# model development

m1<-glm(countevent~ roundtime + age + pgen + previous_dmards + 
          disdur + firsttreat + ovmean + dascore + steroids +
          bmi + renal + lung + diabetes + smoke + offset(log(timevent)), family = poisson(link="log"), data = pdata)

m2<-stpm2(Surv(timevent, event)~age + pgen + previous_dmards + 
            disdur + firsttreat + ovmean + dascore + steroids +
            bmi + renal + lung + diabetes + smoke, k=2, data = fdata)

m3<-glm(event ~ roundtime + age + pgen + previous_dmards + 
          disdur + firsttreat + ovmean + dascore + steroids +
          bmi + renal + lung + diabetes + smoke + offset(log(timevent)), family = binomial(link="cloglog"), data = cdata)

m3.2<-glm(event ~ age + pgen + previous_dmards + 
            disdur + firsttreat + ovmean + dascore + steroids +
            bmi + renal + lung + diabetes + smoke + offset(log(timevent)), family = binomial(link="cloglog"), data = cdata)

m3.3<-glm((event==0) ~ age + pgen + previous_dmards + 
            disdur + firsttreat + ovmean + dascore + steroids +
            bmi + renal + lung + diabetes + smoke + offset(log(timevent)), family = binomial(link="cloglog"), data = cdata)

m4<-glm(event ~ age + pgen + previous_dmards + 
          disdur + firsttreat + ovmean + dascore + steroids +
          bmi + renal + lung + diabetes + smoke, family = binomial(link="logit"), data = lrdata)


par(pty="s")
base<- data.frame(age=0,pgen=0, disdur=0, firsttreat="1", dascore=0, ovmean=0, steroids =0, firsttreat=0,lung=0, diabetes=0, lung=0, renal=0, previous_dmards=0, bmi=0, smoke="0")
plot(m2,type="hazard", newdata=base, xlab="Event time (years)")

# model prediction data 

testC<-cdata
testC$roundtime<-rep(1,nrow(testC))
testC$roundtime<-as.factor(testC$roundtime)
testC$timevent <-rep(1,nrow(testC))

PA1<-predict(m1, newdata = testC, type="response")
PA2<-predict(m2, newdata = testC, type = "surv", se.fit = TRUE, full = TRUE)
PA2$prediction<- 1-PA2$Estimate
PA3<-predict(m3, newdata = testC, type="response")
PA3.2<-predict(m3.2, newdata = testC, type="response")
PA4<-predict(m4, newdata = testC, type="response")

# probability of at least one event in 
fpred<-PA2$prediction
ppred<-as.numeric(PA1)
cpred<-as.numeric(PA3)
cpred.2<-as.numeric(PA3.2)
lrpred<-as.numeric(PA4)

cor(fpred, lrpred)
plot(fpred, lrpred, xlab="FPSM predictions", ylab="LR predictions", cex.lab=1.2)


cor(fpred, ppred)
plot(fpred,ppred, xlab="FPSM predictions", lab="Poisson predictions")

cor(fpred, cpred)
plot(fpred,cpred, xlab="FPSM predictions", ylab="GLM cloglog predictions")

cor(fpred, cpred.2)
plot(fpred,cpred.2, xlab="FPSM predictions", ylab="GLM cloglog predictions w/o time since baseline")

cor(lrpred, cpred.2)
plot(lrpred,cpred.2, xlab="LR predictions", ylab="GLM cloglog predictions w/o time since baseline")

f.coef<-c(0.03541817, -0.25485802, 0.35703620, 0.25333840, 0.39018778, 0.54989517, 0.03805300)
p.coef<-m1$coefficients[18:24]
c.coef<-m3$coefficients[18:24]
cor(f.coef, c.coef)

# longitudinal poisson model

p.long<-read.dta("long_GEE.dta")
p.long$smoke<-factor(p.long$smoke, levels=c("0","1","2"), ordered=F)
p.long$firsttreat<-factor(p.long$firsttreat, ordered = F)
p.long$roundtime<-as.factor(p.long$roundtime)


p.gee<-geeglm(serinf~age + pgen + previous_dmards + 
                disdur + firsttreat + ovmean + dascore + steroids +
                bmi + renal + lung + diabetes + smoke + offset(log(interval)), 
                   family=poisson(link = "log"), data=p.long, 
                   id=groupid, waves=pyears, corstr = "exchangeable", 
                   std.err="san.se")

p.gee.2<-geeglm(serinf~ roundtime + age + pgen + previous_dmards + 
                  disdur + firsttreat + ovmean + dascore + steroids +
                  bmi + renal + lung + diabetes + smoke + offset(log(interval)), 
              family=poisson(link = "log"), data=p.long, 
              id=studyno, waves=pyears, corstr = "exchangeable", 
              std.err="san.se")

c.gee<-geeglm(serinf~ age + pgen + previous_dmards + 
                disdur + firsttreat + ovmean + dascore + steroids +
                bmi + renal + lung + diabetes + smoke + offset(log(interval)), 
              family=binomial(link = "cloglog"), data=p.long, 
              id=studyno, waves=pyears, corstr = "exchangeable", 
              std.err="san.se")

c.gee.2<-geeglm(serinf~ roundtime + age + pgen + previous_dmards + 
                  disdur + firsttreat + ovmean + dascore + steroids +
                  bmi + renal + lung + diabetes + smoke + offset(log(interval)), 
                family=binomial(link = "cloglog"), data=p.long, 
                id=studyno, waves=pyears, corstr = "exchangeable", 
                std.err="san.se")


cor(summary(p.gee)$coefficients$Estimate[2:length(summary(p.gee)$coefficients$Estimate)], 
    summary(c.gee)$coefficients$Estimate[2:length(summary(c.gee)$coefficients$Estimate)])


cor(summary(p.gee)$coefficients$Estimate[2:length(summary(p.gee)$coefficients$Estimate)], 
    summary(c.gee.2)$coefficients$Estimate[5:length(summary(c.gee.2)$coefficients$Estimate)])

test0<-read.dta("baseline_binary.dta")
test0 <- test0 %>% mutate(trt2 = ifelse(firsttreat == "2", 1, 0), 
                          trt4 = ifelse(firsttreat == "4", 1, 0),
                          trt1561 = ifelse(firsttreat == "1561", 1, 0),
                          trt10060 = ifelse(firsttreat == "10060", 1, 0),
                          smoke1 = ifelse(smoke == "1", 1, 0),
                          smoke2 = ifelse(smoke == "2", 1, 0))
covariate_matrix <- test0 %>% select(age, pgen, previous_dmards, disdur, trt2, 
                                     trt4, trt1561, trt10060, ovmean, dascore, 
                                     steroids, bmi, renal, lung, diabetes, 
                                     smoke1, smoke2)
covariate_matrix <- as.matrix(covariate_matrix)
intercept <- rep(1, times = length(covariate_matrix[,1]))
covariate_matrix <- cbind(intercept, covariate_matrix)

p.gee.lp0 <- covariate_matrix %*% summary(p.gee)$coefficients$Estimate
p.gee.2.lp0 <- covariate_matrix %*% c(summary(p.gee.2)$coefficients$Estimate[1],
                                      summary(p.gee.2)$coefficients$Estimate[5:length(summary(p.gee.2)$
                                                                                        coefficients$
                                                                                        Estimate)])

c.gee.lp0 <-covariate_matrix %*% summary(c.gee)$coefficients$Estimate
c.gee.2.lp0 <- covariate_matrix %*% c(summary(c.gee.2)$coefficients$Estimate[1],
                                   summary(c.gee.2)$coefficients$Estimate[5:length(summary(c.gee.2)$
                                                                                     coefficients$
                                                                                     Estimate)])


ppred.long<-exp(p.gee.lp0)
ppred.long.2<-exp(p.gee.2.lp0)
test0$cpred.long<- 1-exp(-exp(c.gee.lp0))
test0$cpred.long2<- 1-exp(-exp(c.gee.2.lp0))

cor(ppred.long, test0$cpred.long)
cor(ppred.long.2, test0$cpred.long)

plot(ppred.long, test0$cpred.long, 
     xlab="GEE (Poisson link)\n risk predictions at baseline",
     ylab="GEE (Complementary log-log link)\n risk predictions at baseline",
     main="GEE without control for time since baseline", cex.lab = 1.2)
      

plot(ppred.long.2, test0$cpred.long2, 
     xlab="GEE (Poisson link)\n risk predictions at baseline",
     ylab="GEE (Complementary log-log link)\n risk predictions at baseline",
     main="GEE with control for time since baseline", cex.lab = 1.2)

################################################################################

c1 <- roc(test0$event~test0$cpred.long,ci=TRUE)
c1

# Drawing the ROC curve
pred.full<-prediction(test0$cpred.long, test0$event)
perf<- performance(pred.full, "tpr", "fpr")
par(pty="s")
plot(perf)
x<-seq(0,1,length.out=100)
lines(x,x)

m1 <- glm(test0$event~offset(test0$c.gee.lp0), family="binomial")
m1$coef
confint(m1)

m11 <- glm(test0$event~test0$c.gee.lp0,family="binomial", x=TRUE,y=TRUE)
m11
confint(m11)

# 1 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test0$cprob2,breaks=quantile(test0$cprob2, prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test0,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(cprob2))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

par(pty="s")
plot(obs~exp[,2],xlim=c(0,0.2),ylim=c(0,0.2),col="darkgreen",ylab="Observed",xlab="Expected")
lines(c(0,0.2),c(0,0.2),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="darkgreen")
}
h <- hist(test0$cprob2, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test0$event))~test0$cprob2,span=1))
lines_data <- data.frame(test0$cprob2,obs_all)
lines_data2 <- lines_data[order(test0$cprob2),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend(0.08,0.05,c("Risk groups","Reference line","95% CI","Loess"),col=c("darkgreen","black","green","grey"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA),bty="n")


#cross-sectional model performance
test0$on_steroid_at_baseline<-test0$steroids
test0$timevent <- rep(1,nrow(test0))
PA3.2<-predict(m3.2, newdata = test0, type="response")
test0$prob <- as.numeric(PA3.2)
  
gee.blp0<- gee.blp0 <- (-5.82943) + 0.03311*test0$age -0.28882*test0$pgen
gee.blp0 <- gee.blp0 + 0.27837*test0$ovmean +  0.27282*test0$on_steroid_at_baseline 
gee.blp0 <- gee.blp0 + 0.40291*test0$lung +  0.45326*test0$diabetes + 0.09329*test0$previous_dmards
test0$gee.blp0<-gee.blp0
test0$baseprob <- 1-exp(-exp(gee.blp0))
#check they are the same thing

plot(test0$baseprob, test0$prob)

m2 <- glm(test0$event~offset(test0$gee.blp0), family="binomial")
m2$coef
confint(m2)

m11 <- glm(test0$event~test0$gee.blp0,family="binomial", x=TRUE,y=TRUE)
m11
confint(m11)

# 1 ## Visual assessment of calibration by risk groups

# create 10 risk groups
groups <- cut(test0$prob,breaks=quantile(test0$prob, prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(test0,groups)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(event)))[,2]
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob))
attach(gpdata)
obsn <- table(event,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# Calibration plot, created manually.

par(pty="s")
plot(obs~exp[,2],xlim=c(0,0.2),ylim=c(0,0.2),col="darkgreen",ylab="Observed",xlab="Expected")
lines(c(0,0.2),c(0,0.2),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="darkgreen")
}
h <- hist(test0$prob, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(test0$event))~test0$prob,span=1))
lines_data <- data.frame(test0$prob,obs_all)
lines_data2 <- lines_data[order(test0$prob),] 
lines(lines_data2[,1],lines_data2[,2],col="grey")
legend(0.08,0.05,c("Risk groups","Reference line","95% CI","Loess"),col=c("darkgreen","black","green","grey"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA),bty="n")
