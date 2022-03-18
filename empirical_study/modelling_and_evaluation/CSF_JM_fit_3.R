library(survival)
library(JMbayes)
library(splines)
library(ipred)
library(dplyr)
library(foreign)
library(foreach)
library(doParallel)

JMdata <- read.dta("JMlong_das.dta")
JMdata$smoke<-factor(JMdata$smoke)
JMdata$trtment<-factor(JMdata$trtment)

JMsurv<-read.dta("JMsurv.dta")
JMsurv$smoke<-factor(JMsurv$smoke)
JMsurv$trtment<-factor(JMsurv$trtment)

set.seed(310394)
lmeFit<- lme(dascore~ns(pyears, knots=2), data = JMdata, 
             random = ~ ns(pyears,knots=2) | groupid, control = list(opt = "optim"))

coxFit <- coxph(Surv(timevent, event) ~ age + pgen +  
                  disdur + previous_dmards + trtment + lung + 
                  diabetes + bmi + renal + steroids + ovmean + 
                  smoke, data = JMsurv, x = TRUE)

jointFit<-jointModelBayes(lmeFit, coxFit, timeVar = "pyears")
summary(jointFit)
saveRDS(jointFit, "./jointmodelresults/haq3.rds")
