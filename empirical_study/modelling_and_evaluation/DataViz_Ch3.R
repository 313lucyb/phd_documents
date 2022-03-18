library(foreign)
library(ggplot2)
library(tableone)
library(Hmisc)
library(dplyr)
library(patchwork)
library(RColorBrewer)
cols <-brewer.pal(11, "RdBu")
cols_models <- c(cols[1:4], cols[9:11])

# baseline correlations

baseline_data <- read.dta("baseline_model2.dta")
par(pty="s")
attach(baseline_data)
pairs(cbind(ovmean,age, bmi, dascore, previous_dmards, disdur),
      col=cols[10],labels=c("HAQ","Age","BMI","DAS28", "No.\n prev. \n DMARDs", 
                               "Disease\n duration"), cex.labels=1.5)

# Correlation coefficients
round(rcorr(cbind(ovmean,age, bmi, dascore, previous_dmards, 
                  disdur),type="pearson")$r,2)

# Sanple size and events within 12 months over time

v0 <-read.dta("baseline_model.dta")
v1<-read.dta("validation_6months.dta")
v2<-read.dta("validation_12months.dta")
v3<-read.dta("validation_18months.dta")
v4<-read.dta("validation_24months.dta")
v5<-read.dta("validation_30months.dta")
v6<-read.dta("validation_36months.dta")

s1 <- length(v0$studyno)
e1 <- sum(v0$event)
s2 <- length(v1$studyno)
e2 <- sum(v1$event)
s3 <- length(v2$studyno)
e3 <- sum(v2$event)
s4 <- length(v3$studyno)
e4 <- sum(v3$event)
s5 <- length(v4$studyno)
e5 <- sum(v4$event)
s6 <- length(v5$studyno)
e6 <- sum(v5$event)
s7 <- length(v6$studyno)
e7 <- sum(v6$event)


s1.1<-c(s1,s2,s3,s4,s5,s6,s7)
e1.1<-c(e1,e2,e3,e4,e5,e6, e7)
i1.1<-s1.1-e1.1

tabb1<-data.frame(e1.1,i1.1, row.names = c(0,0.5,1,1.5,2,2.5,3))
tabb1<-t(tabb1)

par(pty="s")
barplot(tabb1, col = c(cols[11], cols[9]), 
        xlab = "Time post-baseline (years)", 
        ylab = "Number of exposed subjects at-risk")
legend("topright", c("Event subjects", "Non-event subjects"), 
       fill = c(cols[11], cols[9]))


# Data Tables

noncomplete<-read.dta("long_noncomplete.dta")
noncompletebaseline<-noncomplete[noncomplete$pyears == 0, ]
VarsToFactor <- c("lung", "renal", "diabetes", "on_steroid_at_baseline", "pgen", "trtment", "smoke", "event", "died", "numev")
noncompletebaseline[VarsToFactor]<-lapply(noncompletebaseline[VarsToFactor], factor)

dput(names(noncompletebaseline))
vars = c("timevent", "event", "numev", "died", "age", "pgen", "trtment", "disdur","previous_dmards",
         "on_steroid_at_baseline", "ovmean", "dascore", "bmi", "smoke", "renal", "lung", "diabetes")
nonnormal = c("ovmean", "previous_dmards")
tableOne<-CreateTableOne(vars = vars, strata = c("event"), addOverall = TRUE, includeNA = TRUE, data=noncompletebaseline)
print(tableOne, nonnormal = nonnormal, missing = TRUE)

complete<-read.dta("long_complete.dta")
completebaseline<-complete[complete$pyears == 0, ]
VarsToFactor <- c("lung", "renal", "diabetes", "on_steroid_at_baseline", "pgen", "trtment", "smoke", "event", "died", "numev")
completebaseline[VarsToFactor]<-lapply(completebaseline[VarsToFactor], factor)

dput(names(completebaseline))
vars = c("timevent", "event", "numev", "died", "age", "pgen", "trtment", "disdur","previous_dmards",
         "on_steroid_at_baseline", "ovmean", "dascore", "bmi", "smoke", "renal", "lung", "diabetes")
nonnormal = c("ovmean", "previous_dmards")
tableOne<-CreateTableOne(vars = vars, strata = c("event"), addOverall = TRUE, includeNA = TRUE, data=completebaseline)
print(tableOne, nonnormal = nonnormal, noSpaces = TRUE)

# Temporal assessment datasets

LT = c(0, 0.5, 1, 1.5, 2, 2.5, 3)
LRM.TA<-read.table("lrmTA.csv", header = TRUE, sep = ",", row.names = 1)
FPM1.TA<-read.table("fpmTA_1.csv", header = TRUE, sep = ",", row.names = 1)
FPM2.TA<-read.table("fpmTA_2.csv", header = TRUE, sep = ",", row.names = 1)
TDCM.TA<-read.table("tdcmTA.csv", header = TRUE, sep = ",", row.names = 1)
GEE.TA<-read.table("geeTA.csv", header = TRUE, sep = ",", row.names = 1)
LM.TA<-read.table("lmTA.csv", header = TRUE, sep = ",", row.names = 1)
LM.TA$LM <- rep(LT, times = 3)
SRF.TA<-read.table("srfTA.csv", header = TRUE, sep = ",", row.names = 1)
JM1.TA<-read.table("JM1_TA.csv", header = TRUE, sep = ",", row.names = 1)
JM2.TA<-read.table("JM2_TA.csv", header = TRUE, sep = ",", row.names = 1)
JM3.TA<-read.table("JM3_TA.csv", header = TRUE, sep = ",", row.names = 1)
JM4.TA<-read.table("JM4_TA.csv", header = TRUE, sep = ",", row.names = 1)
RRS.TA<-read.table("rrsTA.csv", header = TRUE, sep = ",", row.names = 1)


# fit times and prediction times

FPM1_fit <- FPM1.TA$fit_times_one[1]/60
FPM2_fit <- FPM2.TA$fit_times_two[1]/60
FPM_fit <- mean(c(FPM1_fit, FPM2_fit))
GEE_fit <- mean(GEE.TA$fit_times)/60
LRM_fit <- LRM.TA$fit_vec[1]/60
LM_fit <- mean(LM.TA$fit_times1/60)
JM_fit <- mean(c(49.9,53.5,57.1,60))
SRF_fit <- SRF.TA$fit_times[1]
models <- c("LRM","GEE","FPSM","SRF","TDCM","LM","JM")
fit_times <- c(LRM_fit, GEE_fit, FPM_fit, SRF_fit, TDCM.TA$fit_times[1]/60, LM_fit, JM_fit)
fit_data <- data.frame(model = models, time = fit_times)

names(FPM2.TA) <- names(FPM1.TA)
FPM.TA <- rbind(FPM1.TA, FPM2.TA)
FPM_pred <- FPM.TA %>% group_by(LT) %>% summarise(pred_time = mean(pred_times_one))
GEE_pred <- GEE.TA %>% group_by(LT0) %>% summarise(pred_time = mean(pred_times0))
LM_pred <- LM.TA %>% group_by(LM) %>% summarise(pred_time = mean(pred_times1))
JM.TA <- rbind(JM1.TA, JM2.TA, JM3.TA, JM4.TA)
JM_pred <- JM.TA %>% group_by(X.2) %>% summarise(pred_time = mean(X.3))

LT_pred = rep(LT, times = 7)
models_pred <- rep(models, each=7)
times_pred <- c(LRM.TA$pred_time/60, GEE_pred$pred_time/60, FPM_pred$pred_time/60, SRF.TA$pred_times/60,
                TDCM.TA$pred_times/60, LM_pred$pred_time/60,  JM_pred$pred_time*60)
pred_time_data<-data.frame(LT = LT_pred, model = models_pred, time = times_pred)


fit_data$model<-factor(fit_data$model, levels= c("LRM","GEE","FPSM","SRF","TDCM","LM","JM"))
fit_times <- ggplot(data=fit_data, aes(x=model, y=time, fill=model)) + geom_bar(stat="identity") +
  ylab("Model fitting time (minutes)") + xlab("Model") + scale_fill_manual("Models", values=cols_models) +
  theme_bw()

pred_no_jm <- pred_time_data[1:42,]
pred_times_no_jm <- ggplot(data=pred_no_jm, aes(x=LT, y=time, fill=model)) + geom_bar(stat="identity", position="dodge") +
  ylab("Prediction times \n (minutes)") + xlab("Landmark time \n (years post-baseline)") + scale_fill_manual("Models", values=cols_models) +
  theme_bw() + theme(legend.position = "none")

pred_time_data$model<-factor(pred_time_data$model, levels= c("LRM","GEE","FPSM","SRF","TDCM","LM","JM"))
pred_times <- ggplot(data=pred_time_data, aes(x=LT, y=time, fill=model)) + geom_bar(stat="identity", position="dodge") +
  ylab("Prediction times \n (minutes)") + xlab("Landmark time \n (years post-baseline)") + scale_fill_manual("Models", values=cols_models) +
  theme_bw()

# IPA calculations
null_binary<-read.table("nullBrier_binary.csv", header = TRUE, sep = ",")
null_survival<-read.table("nullBrier_survival.csv", header = TRUE, sep = ",")

IPA_LRM <- 1 - (as.numeric(LRM.TA$BS) / as.numeric(null_binary[,2]))
LRM.TA$IPA <- 100*IPA_LRM
null_binary_gee <- rep(null_binary[,2], each=4)
IPA_GEE <- 1 - (as.numeric(GEE.TA$BS0) / as.numeric(null_binary_gee))
GEE.TA$IPA <- 100*IPA_GEE

IPA_FPSM1 <- 1- (as.numeric(FPM1.TA$BS) / as.numeric(null_survival[,2]))
FPM1.TA$IPA <- 100*IPA_FPSM1
IPA_FPSM2 <- 1 - (as.numeric(FPM2.TA$BS) / as.numeric(null_survival[,2]))
FPM2.TA$IPA <- 100*IPA_FPSM2
IPA_SRF <- 1- (as.numeric(SRF.TA$BS) / as.numeric(null_survival[,2]))
SRF.TA$IPA <- 100*IPA_SRF
IPA_TDCM <- 1- (as.numeric(TDCM.TA$BS) / as.numeric(null_survival[,2]))
TDCM.TA$IPA <- 100*IPA_TDCM
null_survival_lm <- rep(null_survival[,2], times = 3)
IPA_LM <- 1- (as.numeric(LM.TA$BS1) / as.numeric(null_survival_lm))
LM.TA$IPA <- 100*IPA_LM
null_survival_jm <- rep(null_survival[,2], times = 4)
IPA_JM <- 1- (as.numeric(JM.TA$Brier.score) / as.numeric(null_survival_jm))
JM.TA$IPA <- 100*IPA_JM

#within-method comparison LRM


Cind = ggplot(data=LRM.TA,
                 aes(x = LT,y =Cstat, ymin =CstatL, ymax =CstatU))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("AUC (95% CI)")+
  geom_errorbar(aes(ymin=CstatL, ymax=CstatU),width=0.2,cex=1, col = cols_models[1])+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  geom_pointrange(col = cols_models[1]) +  theme_bw() +
  theme(text = element_text(size=14))
Cind

Cslope = ggplot(data=LRM.TA,
                   aes(x = LT,y =Cslope, ymin =CslopeL, ymax =CslopeU))+
  geom_errorbar(aes(ymin=CslopeL, ymax=CslopeU),width=0.2,cex=1, col = cols_models[1])+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  geom_pointrange(col = cols_models[1])+
  geom_hline(yintercept =1, linetype=2)+ theme_bw() +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Calibration slope \n (95% CI)")+
  theme(aspect.ratio=1, 
        text = element_text(size=14))
Cslope

Clarge = ggplot(data=LRM.TA,
                   aes(x = LT,y =Clarge, ymin =ClargeL, ymax =ClargeU))+
  geom_errorbar(aes(ymin=ClargeL, ymax=ClargeU),width=0.2,cex=1, col = cols_models[1])+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  geom_pointrange(col = cols_models[1])+
  geom_hline(yintercept =0, linetype=2)+
  xlab('Prediction time \n (years post-baseline)')+ ylab("Calibration-in-the-large \n (95% CI)")+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold")) + theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
Clarge

Bs<-ggplot(data=LRM.TA,
              aes(x = LT,y =BS))+
  geom_point(col = cols_models[1])+geom_line(col = cols_models[1]) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) + theme_bw() +
  theme(aspect.ratio=1, text = element_text(size=14))
Bs

ipa<-ggplot(data=LRM.TA,
           aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA), colour = cols_models[1])+
  geom_point(col=cols_models[1], size=2)+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA \n (%)")+ 
  geom_hline(yintercept =0, linetype=2)+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipa


Cind / (Cslope | Clarge) / (Bs | ipa)

#within-method comparison TDCM


Cind = ggplot(data=TDCM.TA,
                 aes(x = LT,y =Cind, ymin =CindL, ymax =CindU))+
  xlab('Prediction time (years post-baseline)')+ ylab("C index (95% CI)")+
  geom_errorbar(aes(ymin=CindL, ymax=CindU),width=0.2,cex=1, col = cols_models[5])+
  geom_pointrange(col = cols_models[5])+ theme_bw() +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(aspect.ratio = 1, text = element_text(size=14))
Cind

Cslope = ggplot(data=TDCM.TA,
                   aes(x = LT,y =Cslope, ymin =CslopeL, ymax =CslopeU))+
  geom_hline(yintercept =1, linetype=2)+
  geom_errorbar(aes(ymin=CslopeL, ymax=CslopeU),width=0.2,cex=1, col = cols_models[5])+
  geom_pointrange(col = cols_models[5])+ theme_bw() +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Calibration slope \n (95% CI)")+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(aspect.ratio=1, 
        text = element_text(size=14))
Cslope


Bs<-ggplot(data=TDCM.TA,
              aes(x = LT,y =BS))+
  geom_point(col = cols_models[5])+ geom_line(col = cols_models[5]) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold")) + theme_bw() +
    theme(aspect.ratio=1, 
          text = element_text(size=14)) 
Bs

ipa<-ggplot(data=TDCM.TA,
            aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA), colour = cols_models[5])+
  geom_point(col=cols_models[5], size=2)+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA (%)")+ 
  geom_hline(yintercept =0, linetype=2)+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipa

(Cind | Cslope) / (Bs | ipa)

#within-method comparison SRF

Cind = ggplot(data=SRF.TA,
                 aes(x = LT,y =Cind, ymin =CindL, ymax =CindU))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("C index (95% CI)")+
  geom_errorbar(aes(ymin=CindL, ymax=CindU),width=0.2,cex=1, col = cols[3])+
  geom_pointrange(col = cols[3])+ theme_bw() +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(text = element_text(size=14))
Cind


Bs<-ggplot(data=SRF.TA,
              aes(x = LT,y =BS))+
  geom_point(col = cols[3])+ geom_line(col = cols[3]) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ 
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold")) + theme_bw() +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(aspect.ratio=1, 
        text = element_text(size=14))
Bs

ipa<-ggplot(data=SRF.TA,
            aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA), colour = cols_models[3])+
  geom_point(col=cols_models[3], size=2)+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA (%)")+ 
  geom_hline(yintercept =0, linetype=2)+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipa

Cind / (Bs | ipa )

# Within method comparisons - GEE

GEE.TA$Mod<-factor(GEE.TA$Mod, levels=c(1,2,3,4))

CindGEE = ggplot(data=GEE.TA,
           aes(x = LT0,y =c0, ymin =c0.L, ymax =c0.U))+
  geom_errorbar(aes(ymin=c0.L, ymax=c0.U, col=Mod),width=0.2,cex=1, position = position_dodge(0.25))+ 
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  scale_color_manual("GEE\nvariations", values=cols[1:4], labels=c("without", "categorical", "linear", "nonlinear")) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("AUC (95% CI)")+ theme_bw() +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(axis.title.x=element_blank(), 
          text = element_text(size=14))
CindGEE

CslopeGEE = ggplot(data=GEE.TA,
                 aes(x = LT0,y =Cslope0, ymin =Cslope0.L, ymax =Cslope0.U))+
  geom_errorbar(aes(ymin=Cslope0.L, ymax=Cslope0.U, col=Mod),width=0.2,cex=1, position = position_dodge(0.25))+
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  geom_hline(yintercept =1, linetype=2)+
  scale_color_manual("GEE\nvariations", values=cols[1:4], labels=c("without", "categorical", "linear", "nonlinear")) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Calibration slope (95% CI)")+ theme_bw() +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(aspect.ratio=1,
        text = element_text(size=14), legend.position = "none")
CslopeGEE

ClargeGEE = ggplot(data=GEE.TA,
                   aes(x = LT0,y =Clge0, ymin =Clge0.L, ymax =Clge0.U))+
  geom_errorbar(aes(ymin=Clge0.L, ymax=Clge0.U, col=Mod),width=0.2,cex=1, position = position_dodge(0.25))+
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  geom_hline(yintercept =0, linetype=2)+
  geom_vline(xintercept =0.25, linetype=2)+
  geom_vline(xintercept =0.75, linetype=2)+
  geom_vline(xintercept =1.25, linetype=2)+
  geom_vline(xintercept =1.75, linetype=2)+
  geom_vline(xintercept =2.25, linetype=2)+
  geom_vline(xintercept =2.75, linetype=2)+
  scale_color_manual("GEE\nvariations", values=cols[1:4], labels=c("without", "categorical", "linear", "nonlinear")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Calibration-in-the-large \n (95% CI)")+ theme_bw() +
  theme(aspect.ratio=1, 
        text = element_text(size=14))
ClargeGEE

BsGEE<-ggplot(data=GEE.TA,
              aes(x = LT0,y =BS0))+
  geom_point(aes(col = Mod), position = position_dodge(0.25))+ geom_line(aes(col = Mod), position = position_dodge(0.25)) +
  scale_color_manual("GEE\nvariations", values=cols[1:4], labels=c("without", "categorical", "linear", "nonlinear")) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ theme_bw() +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(aspect.ratio=1, 
        text = element_text(size=14), legend.position = "none")
BsGEE

ipa<-ggplot(data=GEE.TA,
            aes(x = LT0,y =IPA))+
  geom_linerange(aes(x = LT0, ymin = 0, ymax = IPA, colour=Mod), position = position_dodge(0.25))+
  geom_point(aes(col=Mod), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA\n (%)")+
  scale_color_manual("GEE\nvariations", values=cols[1:4], labels=c("without", "categorical", "linear", "nonlinear")) +
  geom_hline(yintercept =0, linetype=2)+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipa

CindGEE / ( CslopeGEE | ClargeGEE) / (BsGEE | ipa)

# Within-method comparisons - LM

LM.TA$Mod<-c(rep(1,7),rep(2,7),rep(3,7))
LM.TA$Mod<-factor(LM.TA$Mod, levels = c(1,2,3))

CindLM = ggplot(data=LM.TA,
                 aes(x = LM,y =Cind1, ymin =Cind1L, ymax =Cind1U))+
  geom_errorbar(aes(ymin=Cind1L, ymax=Cind1U, col=Mod),width=0.2,cex=1, position = position_dodge(0.25))+ 
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  scale_color_manual("LM\nvariations", values=cols[9:11], labels=c("0.1", "events", "0.5")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  xlab('Prediction time (years post-baseline)')+ ylab("C index (95% CI)")+ theme_bw() +
  theme(text = element_text(size=14), legend.position="none")
CindLM

CslopeLM = ggplot(data=LM.TA,
                   aes(x = LM,y =Cslope1, ymin =Cslope1L, ymax =Cslope1U))+
  geom_errorbar(aes(ymin=Cslope1L, ymax=Cslope1U, col=Mod),width=0.2,cex=1, position = position_dodge(0.25))+
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  geom_hline(yintercept =1, linetype=2)+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Calibration slope \n (95% CI)")+ theme_bw() +
  scale_color_manual("LM\nvariations", values=cols[9:11], labels=c("0.1", "events", "0.5")) +
  theme(aspect.ratio=1,
        text = element_text(size=14))
CslopeLM

BsLM<-ggplot(data=LM.TA,
              aes(x = LM,y =BS1))+
  geom_point(aes(col = Mod), position = position_dodge(0.25))+ geom_line(aes(col = Mod), position = position_dodge(0.25)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ theme_bw() +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  scale_color_manual("LM\nvariations", values=cols[9:11], labels=c("0.1", "events", "0.5")) +
  theme(aspect.ratio=1, 
        text = element_text(size=14), legend.position="none")
BsLM

ipa<-ggplot(data=LM.TA,
            aes(x = LM,y =IPA))+
  geom_linerange(aes(x = LM, ymin = 0, ymax = IPA, colour=Mod), position = position_dodge(0.25))+
  geom_point(aes(col=Mod), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA\n (%)")+
  scale_color_manual("LM\nvariations", values=cols[9:11], labels=c("0.1", "events", "0.5")) +
  geom_hline(yintercept =0, linetype=2)+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipa

(CindLM | CslopeLM) / (BsLM | ipa)

# Within-method comparisons - FPM
FPM.TA <- rbind(FPM1.TA, FPM2.TA)
FPM.TA$Mod<-c(rep(1,7), rep(2,7))
FPM.TA$Mod<-factor(FPM.TA$Mod, levels=c(1,2))

CindFPM = ggplot(data=FPM.TA,
                aes(x = LT,y =Cstat, ymin =CstatL, ymax =CstatU))+
  geom_errorbar(aes(ymin=CstatL, ymax=CstatU, col=Mod),width=0.2,cex=1, position = position_dodge(0.25))+
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  scale_color_manual("Year\ncensored", values=c(cols_models[2],cols_models[6]), labels=c("1", "4")) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("C index (95% CI)")+ theme_bw() +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(text = element_text(size=14), legend.position="none")

CindFPM


CslopeFPM = ggplot(data=FPM.TA,
                  aes(x = LT,y =Cslope, ymin =CslopeL, ymax =CslopeU))+
  geom_errorbar(aes(ymin=CslopeL, ymax=CslopeU, col=Mod),width=0.2,cex=1, position = position_dodge(0.25))+
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  geom_hline(yintercept =1, linetype=2)+
  xlab('Prediction time \n (years post-baseline)')+ ylab("Calibration slope (95% CI)")+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  scale_color_manual("Year\ncensored", values=c(cols_models[2],cols_models[6]), labels=c("1", "4")) + theme_bw() +
  theme(aspect.ratio = 1,  text = element_text(size=14))
CslopeFPM

BsFPM<-ggplot(data=FPM.TA,
             aes(x = LT,y =BS))+
  geom_point(aes(col = Mod), position = position_dodge(0.25))+ geom_line(aes(col = Mod), position = position_dodge(0.25)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  scale_color_manual("Year\ncensored", values=c(cols_models[2],cols_models[6]), labels=c("1", "4")) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size=14), legend.position="none")
BsFPM

ipa<-ggplot(data=FPM.TA,
            aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA, colour=Mod), position = position_dodge(0.25))+
  geom_point(aes(col=Mod), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA (%)")+
  scale_color_manual("Year\ncensored", values=c(cols_models[2],cols_models[6]), labels=c("1", "4")) +
  geom_hline(yintercept =0, linetype=2)+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipa

(CindFPM | CslopeFPM) / ( BsFPM | ipa)

# Within-method comparisons - JM

JM.TA$Model<-rep(c(1,2,3,4), each = 7)
JM.TA$Model <- factor(JM.TA$Model, levels=c(1,2,3,4))

Cind = ggplot(data=JM.TA,
                 aes(x = X.2,y =C, ymin =C.1, ymax =C.2))+
  geom_errorbar(aes(ymin=C.1, ymax=C.2, col=Model),width=0.2,cex=1, position = position_dodge(0.25))+
  geom_pointrange(aes(col=Model), position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("C index (95% CI)")+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  scale_color_manual("JM\nvariations", values=c(cols_models[1:2],cols_models[6:7]), labels=c("HAQ (3,3)", "HAQ (3, int)", "DAS (2,2)", "DAS (2, int)")) +
  theme_bw() +
  theme(text = element_text(size=14))
Cind


Bs<-ggplot(data=JM.TA,
              aes(x = X.2,y =Brier.score))+
  geom_point(aes(col = Model), position = position_dodge(0.25))+ geom_line(aes(col = Model), position = position_dodge(0.25)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  scale_color_manual("JM\nvariations", values=c(cols_models[1:2],cols_models[6:7]), labels=c("HAQ (3,3)", "HAQ (3, int)", "DAS (2,2)", "DAS (2, int)")) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size=14), legend.position="none")
Bs

ipa<-ggplot(data=JM.TA,
            aes(x = X.2,y =IPA))+
  geom_linerange(aes(x = X.2, ymin = 0, ymax = IPA, colour=Model), position = position_dodge(0.25))+
  geom_point(aes(col=Model), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA (%)")+
  geom_hline(yintercept =0, linetype=2)+
  scale_color_manual("JM\nvariations", values=c(cols_models[1:2],cols_models[6:7]), labels=c("HAQ (3,3)", "HAQ (3, int)", "DAS (2,2)", "DAS (2, int)")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipa

Cind / ( Bs | ipa )

## RRS within-method comparison

RRS.TA$Model<-factor(RRS.TA$Scenario, levels=c(1,2,3))
null_brier_rrs <- rep(null_survival[,2], times = 3)
RRS.TA$IPA <- 100*(1-(RRS.TA$Brier/null_brier_rrs))

CindGEE = ggplot(data=RRS.TA,
                 aes(x = LT,y =C, ymin =C.low, ymax =C.high))+
  geom_errorbar(aes(ymin=C.low, ymax=C.high, col=Model),width=0.2,cex=1, position = position_dodge(0.25))+ 
  geom_pointrange(aes(col=Model), position = position_dodge(0.25))+
  scale_color_manual("Steroid\nassumptions", values=c(cols[1],cols[3],cols[9]), labels=c("low", "high", "mixed")) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("AUC (95% CI)")+ theme_bw() +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(text = element_text(size=14), legend.position = "none")
CindGEE


ClargeGEE = ggplot(data=RRS.TA,
                   aes(x = LT,y =Clarge, ymin =Clarge.low, ymax =Clarge.high))+
  geom_errorbar(aes(ymin=Clarge.low, ymax=Clarge.high, col=Model),width=0.2,cex=1, position = position_dodge(0.25))+
  geom_pointrange(aes(col=Model), position = position_dodge(0.25))+
  geom_hline(yintercept =0, linetype=2)+
  geom_vline(xintercept =0.25, linetype=2)+
  geom_vline(xintercept =0.75, linetype=2)+
  geom_vline(xintercept =1.25, linetype=2)+
  geom_vline(xintercept =1.75, linetype=2)+
  geom_vline(xintercept =2.25, linetype=2)+
  geom_vline(xintercept =2.75, linetype=2)+
  scale_color_manual("Steroid\nassumptions", values=c(cols[1],cols[3],cols[9]), labels=c("low", "high", "mixed")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Calibration-in-the-large \n (95% CI)")+ theme_bw() +
  theme(aspect.ratio=1, 
        text = element_text(size=14))
ClargeGEE

BsGEE<-ggplot(data=RRS.TA,
              aes(x = LT,y =Brier))+
  geom_point(aes(col = Model), position = position_dodge(0.25))+ geom_line(aes(col = Model), position = position_dodge(0.25)) +
  scale_color_manual("Steroid\nassumptions", values=c(cols[1],cols[3],cols[9]), labels=c("low", "high", "mixed")) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ theme_bw() +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(aspect.ratio=1, 
        text = element_text(size=14), legend.position = "none")
BsGEE

ipa<-ggplot(data=RRS.TA,
            aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA, colour=Model), position = position_dodge(0.25))+
  geom_point(aes(col=Model), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA (%)")+
  geom_hline(yintercept =0, linetype=2)+
  scale_color_manual("Steroid\nassumptions", values=c(cols[1],cols[3],cols[9]), labels=c("low", "high", "mixed")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipa

(CindGEE | ClargeGEE) / (BsGEE | ipa)

# Between-method comparisons

LRComp<-LRM.TA[, c(1,c(8:11,14))]
FPMComp<-FPM.TA[(8:14), c(1,(5:8), 11)]
TDCMComp<-TDCM.TA[, c(1,(5:8),11)]
GEE.TA.2 <- GEE.TA %>% filter(Mod == 1)
GEEComp<-GEE.TA.2[, c(1, (4:6), 3, 15)]
LM.TA.2 <- LM.TA %>% filter(Mod == 3)
LMComp<-LM.TA.2[, c(10, 4:7,11)]
SRFComp<-SRF.TA[, c(1:5, 8)]
JM.TA.2 <- JM.TA %>% filter(Model == 2)
JMComp<- JM.TA.2[,c(2:6, 8)]

Binary_Comp<-rbind(as.matrix(LRComp), as.matrix(GEEComp))
Survival_Comp<-rbind(as.matrix(FPMComp), as.matrix(SRFComp), as.matrix(TDCMComp), 
              as.matrix(LMComp), as.matrix(JMComp))

Binary_Comp<-as.data.frame(Binary_Comp)
Survival_Comp<-as.data.frame(Survival_Comp)

Mod_bin<-c(rep(1,7), rep(2,7))
Mod_surv <- c(rep(3,7), rep(4,7), rep(5,7), rep(6,7), rep(7,7))

Survival_Comp$Mod <-Mod_surv
Binary_Comp$Mod <- Mod_bin

# Adding in RRS data

RRS.TA.2 <- RRS.TA %>% filter(Scenario == 1)
RRSComp <-RRS.TA.2[, c(3:6, 10,12)]
RRSComp$Mod <- rep(8, times = 7)
names(RRSComp) <- names(Binary_Comp)
Binary_Comp<-rbind(Binary_Comp, RRSComp)
Binary_Comp$Mod <- factor(Binary_Comp$Mod, levels <- c(1,2,8))
Survival_Comp$Mod <- factor(Survival_Comp$Mod, levels <- c(3,4,5,6,7))

#---------------------------

Cindcomp = ggplot(data=Binary_Comp,
                  aes(x = LT,y =Cstat, ymin =CstatL, ymax =CstatU))+
  geom_errorbar(aes(ymin=CstatL, ymax=CstatU, col=Mod),width=0.2,cex=1, position = position_dodge(0.25))+ 
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("AUC (95% CI)")+
 scale_colour_manual(name=NULL, values=c(cols_models[1],cols_models[2],"black"), labels=c("LRM","GEE", "RRS")) +theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size=14), legend.position = "none")
Cindcomp


Cindcomp2 = ggplot(data=Survival_Comp,
                  aes(x = LT,y =Cstat, ymin =CstatL, ymax =CstatU))+
  geom_errorbar(aes(ymin=CstatL, ymax=CstatU, col=Mod),width=0.2, cex=1, position = position_dodge(0.25))+ 
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Harrell's C-index (95% CI)")+
  scale_colour_manual(name=NULL, values=cols_models[3:7], labels=c("FPSM","SRF", "TDCM", "LM", "JM")) +theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size=14), legend.position = "none")
Cindcomp2

Cindcomp / Cindcomp2

BsComp<-ggplot(data=Binary_Comp,
              aes(x = LT,y =BS))+
  geom_point(aes(col = Mod), position = position_dodge(0.25))+ geom_line(aes(col = Mod), position = position_dodge(0.25)) +
  scale_colour_manual(name=NULL, values=c(cols_models[1],cols_models[2],"black"), labels=c("LRM","GEE", "RRS")) +theme_bw() +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(aspect.ratio = 1, text = element_text(size=14), legend.position = "none")
BsComp

BsComp2<-ggplot(data=Survival_Comp,
               aes(x = LT,y =BS))+
  geom_point(aes(col = Mod), position = position_dodge(0.25))+ geom_line(aes(col = Mod), position = position_dodge(0.25)) +
  scale_colour_manual(name=NULL, values=cols_models[3:7], labels=c("FPSM","SRF", "TDCM", "LM", "JM")) +theme_bw() +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(text = element_text(size=14), legend.position = "none")
BsComp2

ipacomp<-ggplot(data=Binary_Comp,
            aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA, colour=Mod), position = position_dodge(0.25))+
  geom_point(aes(col=Mod), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA (%)")+
  geom_hline(yintercept =0, linetype=2)+
  scale_colour_manual(name=NULL, values=c(cols_models[1],cols_models[2],"black"), labels=c("LRM","GEE", "RRS")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipacomp

ipacomp2<-ggplot(data=Survival_Comp,
                aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA, colour=Mod), position = position_dodge(0.25))+
  geom_point(aes(col=Mod), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA (%)")+
  geom_hline(yintercept =0, linetype=2)+
  scale_colour_manual(name=NULL, values=cols_models[3:7], labels=c("FPSM","SRF", "TDCM", "LM", "JM")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipacomp2
Cindcomp | BsComp | ipacomp
Cindcomp2 |BsComp2 | ipacomp2



# Optimism estimate

LRM.opt<-read.table("lrmOPT.csv", header = TRUE, sep = ",", row.names = 1)
FPM.opt<-read.table("fpm1OPT.csv", header = TRUE, sep = ",", row.names = 1)

SRF.opt<-read.table("srfOPT.csv", header = TRUE, sep = ",", row.names = 1)
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


TDCM.opt<-read.table("tdcmOPT.csv", header = TRUE, sep = ",", row.names = 1)
GEE.opt<-read.table("geeOPT.csv", header = TRUE, sep = ",", row.names = 1)
LM.opt<-read.table("lmOPT.csv", header = TRUE, sep = ",", row.names = 1)

JM1.opt<-read.table("JM_optimism_results_updated_1.csv", header = TRUE, sep = ",")
JM2.opt<-read.table("JM_optimism_results_updated_2.csv", header = TRUE, sep = ",")
JM3.opt<-read.table("JM_optimism_results_updated_3.csv", header = TRUE, sep = ",")
JM4.opt<-read.table("JM_optimism_results_updated_4.csv", header = TRUE, sep = ",")
JM5.opt<-read.table("JM_optimism_results_updated_5.csv", header = TRUE, sep = ",")
JM6.opt<-read.table("JM_optimism_results_updated_6.csv", header = TRUE, sep = ",")
JM7.opt<-read.table("JM_optimism_results_updated_7.csv", header = TRUE, sep = ",")
JM8.opt<-read.table("JM_optimism_results_updated_8.csv", header = TRUE, sep = ",")
JM9.opt<-read.table("JM_optimism_results_updated_9.csv", header = TRUE, sep = ",")
JM10.opt<-read.table("JM_optimism_results_updated_10.csv", header = TRUE, sep = ",")

JM.opt<-rbind(JM1.opt, 
              JM2.opt,
              JM3.opt,
              JM4.opt,
              JM5.opt,
              JM6.opt,
              JM7.opt,
              JM8.opt,
              JM9.opt,
              JM10.opt)
names(JM.opt) <- c("row_num","iteration", "LT", "C", "C_low", "C_high", "Brier")
JM.opt2 <- JM.opt %>% group_by(LT) %>% summarise(C_mean = mean(C), Brier_mean = mean(Brier))

JM.opt2 <- as.data.frame(JM.opt2)
JM.OPT<-data.frame(LT0 = as.numeric(JM.opt2[1,(2:3)]), LT0.5 = as.numeric(JM.opt2[2,(2:3)]), LT1 = as.numeric(JM.opt2[3,(2:3)]), 
                    LT1.5 = as.numeric(JM.opt2[4,(2:3)]), LT2 = as.numeric(JM.opt2[5,(2:3)]), LT2.5 = as.numeric(JM.opt2[6,(2:3)]), 
                    LT3 = as.numeric(JM.opt2[7,(2:3)]), row.names = c("C statistic (OE)", 
                                                      "Brier Score (OE)"))

optimism<-rbind(LRM.opt, GEE.opt, FPM.opt, SRF.opt, TDCM.opt, LM.opt, JM.OPT)
optimism.C<-optimism[seq(1,13,2),]
optimism.B<-optimism[seq(2,14,2),]
model<-c(1,2,3,4,5,6,7)
optimism.C<-cbind(model, optimism.C)
optimism.B<-cbind(model, optimism.B)


LT0.C<-as.matrix(optimism.C[,1:2])
LT0.5.C<-as.matrix(optimism.C[,(c(1,3))])
LT1.C<-as.matrix(optimism.C[,(c(1,4))])
LT1.5.C<-as.matrix(optimism.C[,(c(1,5))])
LT2.C<-as.matrix(optimism.C[,(c(1,6))])
LT2.5.C<-as.matrix(optimism.C[,(c(1,7))])
LT3.C<-as.matrix(optimism.C[,(c(1,8))])
LT<-c(rep(0,7),rep(0.5,7),rep(1,7),rep(1.5,7),rep(2,7),rep(2.5,7),rep(3,7))
optC<-rbind(LT0.C, LT0.5.C, LT1.C, LT1.5.C, LT2.C, LT2.5.C, LT3.C)
optC<-as.data.frame(optC)
optC<-cbind(LT,optC)
#optC$LT<-as.factor(optC$LT)
optC$model<-as.factor(optC$model)


optC_bin <- optC %>% filter(model == 1 | model == 2)
OptcompC<-ggplot(data=optC_bin,
                 aes(x = LT,y =LT0))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = LT0, colour = model), 
                 position = position_dodge(0.25))+
  geom_hline(yintercept=0, linetype=2) +
  geom_point(aes(col=model),position = position_dodge(0.25), size=2)+ 
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism estimate (AUC)")+ 
  scale_colour_manual(name=NULL, values=c(cols_models[1],cols_models[2]), labels=c("LRM","GEE")) + theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size=14), legend.position="none") 
OptcompC


optC_surv <- optC %>% filter(model == 3 | model == 4| model == 5| model == 6 | model == 7 )
OptcompC2<-ggplot(data=optC_surv,
                 aes(x = LT,y =LT0))+
  geom_hline(yintercept=0, linetype=2) +
  geom_linerange(aes(x = LT, ymin = 0, ymax = LT0, colour = model), 
                 position = position_dodge(0.25))+
  geom_point(aes(col=model),position = position_dodge(0.25), size=2)+ 
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism estimate \n(Harrell's C index)")+ 
  scale_colour_manual(name=NULL, values=cols_models[3:7], labels=c("FPSM","SRF", "TDCM", "LM", "JM")) + theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size=14), legend.position="none") 
OptcompC2


LT0.B<-as.matrix(optimism.B[,1:2])
LT0.5.B<-as.matrix(optimism.B[,(c(1,3))])
LT1.B<-as.matrix(optimism.B[,(c(1,4))])
LT1.5.B<-as.matrix(optimism.B[,(c(1,5))])
LT2.B<-as.matrix(optimism.B[,(c(1,6))])
LT2.5.B<-as.matrix(optimism.B[,(c(1,7))])
LT3.B<-as.matrix(optimism.B[,(c(1,8))])
LT<-c(rep(0,7),rep(0.5,7),rep(1,7),rep(1.5,7),rep(2,7),rep(2.5,7),rep(3,7))
optB<-rbind(LT0.B, LT0.5.B, LT1.B, LT1.5.B, LT2.B, LT2.5.B, LT3.B)
optB<-as.data.frame(optB)
optB<-cbind(LT,optB)

optB$model<-as.factor(optB$model)

optB_bin <- optB %>% filter(model == 1 | model == 2)
optB_surv <- optB %>% filter(model == 3 | model == 4| model == 5| model == 6 | model == 7)

OptcompB<-ggplot(data=optB_bin,
                 aes(x = LT,y =LT0))+
  geom_hline(yintercept=0, linetype=2) +
  geom_linerange(aes(x = LT, ymin = 0, ymax = LT0, colour = model), 
                 position = position_dodge(0.25))+
  geom_point(aes(col=model),position = position_dodge(0.25), size = 2) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism estimate \n(Brier)")+ 
  scale_colour_manual(name=NULL, values=c(cols_models[1],cols_models[2]), labels=c("LRM","GEE")) + theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size=14)) 
OptcompB

OptcompB2<-ggplot(data=optB_surv,
                  aes(x = LT,y =LT0))+
  geom_hline(yintercept=0, linetype=2) +
  geom_linerange(aes(x = LT, ymin = 0, ymax = LT0, colour = model), 
                 position = position_dodge(0.25))+
  geom_point(aes(col=model),position = position_dodge(0.25), size = 2) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism estimate \n(Brier)")+ 
  scale_colour_manual(name=NULL, values=cols_models[3:7], labels=c("FPSM","SRF", "TDCM", "LM", "JM")) + theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size=14)) 
OptcompB2

OptcompC | OptcompB 
OptcompC2 | OptcompB2



OptB<-optB[with(optB, order(model, LT)), ]
OptC<-optC[with(optC, order(model, LT)), ]
OptB_bin <- filter(OptB, model == 1 | model == 2)
OptB_surv <- OptB %>% filter(model == 3 | model == 4| model == 5| model == 6 | model == 7)
OptC_bin <- filter(OptC, model == 1 | model == 2)
OptC_surv <- OptC %>% filter(model == 3 | model == 4| model == 5| model == 6 | model == 7)

bin_opt<-Binary_Comp[1:14,]
surv_opt<-Survival_Comp[1:35,]

bin_opt$BS<-bin_opt$BS+sqrt(OptB_bin$LT0^2)
bin_opt$Cstat<-bin_opt$Cstat-sqrt(OptC_bin$LT0^2)
bin_opt$CstatU<-bin_opt$CstatU-sqrt(OptC_bin$LT0^2)
bin_opt$CstatL<-bin_opt$CstatL-sqrt(OptC_bin$LT0^2)
null_bin <- rep(null_binary[,2], times = 2)
bin_opt$IPA <- 100*(1-(bin_opt$BS/null_bin))
bin_opt_RRS<-rbind(bin_opt,Binary_Comp[15:21,])

null_surv <- rep(null_survival[,2], times=5)
surv_opt$BS<-surv_opt$BS+sqrt(OptB_surv$LT0^2)
surv_opt$Cstat<-surv_opt$Cstat-sqrt(OptC_surv$LT0^2)
surv_opt$CstatU<-surv_opt$CstatU-sqrt(OptC_surv$LT0^2)
surv_opt$CstatL<-surv_opt$CstatL-sqrt(OptC_surv$LT0^2)
surv_opt$IPA <- 100*(1-(surv_opt$BS/null_surv))

Cindcomp = ggplot(data=bin_opt_RRS,
                  aes(x = LT,y =Cstat, ymin =CstatL, ymax =CstatU))+
  geom_errorbar(aes(ymin=CstatL, ymax=CstatU, col=Mod),width=0.2,cex=1, position = position_dodge(0.25))+ 
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism-adjusted \n AUC (95% CI)")+
  scale_colour_manual(name=NULL, values=c(cols_models[1],cols_models[2],"black"), labels=c("LRM","GEE", "RRS")) +theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size=14), legend.position = "none")
Cindcomp

Cindcomp2 = ggplot(data=surv_opt,
                   aes(x = LT,y =Cstat, ymin =CstatL, ymax =CstatU))+
  geom_errorbar(aes(ymin=CstatL, ymax=CstatU, col=Mod),width=0.2, cex=1, position = position_dodge(0.25))+ 
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism-adjusted \n Harrell's C-index (95% CI)")+
  scale_colour_manual(name=NULL, values=cols_models[3:7], labels=c("FPSM","SRF", "TDCM", "LM", "JM")) +theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size=14), legend.position = "none")
Cindcomp2

Cindcomp / Cindcomp2

BsComp<-ggplot(data=bin_opt_RRS,
               aes(x = LT,y =BS))+
  geom_point(aes(col = Mod), position = position_dodge(0.25))+ geom_line(aes(col = Mod), position = position_dodge(0.25)) +
  scale_colour_manual(name=NULL, values=c(cols_models[1],cols_models[2],"black"), labels=c("LRM","GEE", "RRS")) +theme_bw() +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism-adjusted \n Brier score")+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(aspect.ratio = 1, text = element_text(size=14), legend.position = "none")
BsComp

BsComp2<-ggplot(data=surv_opt,
                aes(x = LT,y =BS))+
  geom_point(aes(col = Mod), position = position_dodge(0.25))+ geom_line(aes(col = Mod), position = position_dodge(0.25)) +
  scale_colour_manual(name=NULL, values=cols_models[3:7], labels=c("FPSM","SRF", "TDCM", "LM", "JM")) +theme_bw() +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism-adjusted \n Brier score")+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(text = element_text(size=14), legend.position = "none")
BsComp2

ipacomp<-ggplot(data=bin_opt_RRS,
                aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA, colour=Mod), position = position_dodge(0.25))+
  geom_point(aes(col=Mod), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism-adjusted \n IPA (%)")+
  geom_hline(yintercept =0, linetype=2)+
  scale_colour_manual(name=NULL, values=c(cols_models[1],cols_models[2],"black"), labels=c("LRM","GEE", "RRS")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipacomp

ipacomp2<-ggplot(data=surv_opt,
                 aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA, colour=Mod), position = position_dodge(0.25))+
  geom_point(aes(col=Mod), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism-adjusted \n IPA (%)")+
  geom_hline(yintercept =0, linetype=2)+
  scale_colour_manual(name=NULL, values=cols_models[3:7], labels=c("FPSM","SRF", "TDCM", "LM", "JM")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipacomp2

Cindcomp | BsComp | ipacomp
Cindcomp2 |BsComp2 | ipacomp2


## Between all methods

LRM.TA.surv <- read.table("lrmTA_surv.csv", header = TRUE, sep = ",", row.names = 1)
LRM.opt.surv<- read.table("lrmOPT_surv.csv", header = TRUE, sep = ",", row.names = 1)
GEE.TA.surv <- read.table("geeTA_surv.csv", header = TRUE, sep = ",", row.names = 1)
GEE.opt.surv<- read.table("geeOPT_surv.csv", header = TRUE, sep = ",", row.names = 1)
RRS.TA.surv <- read.table("rrsTA_surv.csv", header = TRUE, sep = ",", row.names = 1)

# current datasets to compare apparent and optimism-adjusted survival measures
Survival_Comp
surv_opt

LRM.TA.surv$IPA <- 100*(1 - (LRM.TA.surv$BS/null_survival[,2]))
GEE.TA.surv$IPA <- 100*(1 - (GEE.TA.surv$BS/null_survival[,2]))
RRS.TA.surv$IPA <- 100*(1 - (RRS.TA.surv$BS/null_survival[,2]))

LRM.TA.surv$Mod <- rep(1, times = 7)
GEE.TA.surv$Mod <- rep(2, times = 7)
RRS.TA.surv$Mod <- rep(8, times = 7)
Survival_Comp$Mod <-factor(Survival_Comp$Mod, levels=c(1,2,3,4,5,6,7,8))
Survival_Comp_all <- rbind(Survival_Comp, LRM.TA.surv, GEE.TA.surv, RRS.TA.surv)

LRM.C.opt <- as.numeric(LRM.opt.surv[1,])
LRM.BS.opt <- as.numeric(LRM.opt.surv[2,])
LRM_surv_opt <- LRM.TA.surv
LRM_surv_opt$Cstat <- LRM_surv_opt$Cstat - sqrt(LRM.C.opt^2)
LRM_surv_opt$CstatL <- LRM_surv_opt$CstatL - sqrt(LRM.C.opt^2)
LRM_surv_opt$CstatU <- LRM_surv_opt$CstatU - sqrt(LRM.C.opt^2)
LRM_surv_opt$BS <- LRM_surv_opt$BS + sqrt(LRM.BS.opt^2)
LRM_surv_opt$IPA <- 100*(1-(LRM_surv_opt$BS/null_survival[,2]))

GEE.C.opt <- as.numeric(GEE.opt.surv[1,])
GEE.BS.opt <- as.numeric(GEE.opt.surv[2,])
GEE_surv_opt <- GEE.TA.surv
GEE_surv_opt$Cstat <- GEE_surv_opt$Cstat - sqrt(GEE.C.opt^2)
GEE_surv_opt$CstatL <- GEE_surv_opt$CstatL - sqrt(GEE.C.opt^2)
GEE_surv_opt$CstatU <- GEE_surv_opt$CstatU - sqrt(GEE.C.opt^2)
GEE_surv_opt$BS <- GEE_surv_opt$BS + sqrt(GEE.BS.opt^2)
GEE_surv_opt$IPA <- 100*(1-(GEE_surv_opt$BS/null_survival[,2]))

surv_opt$Mod <- factor(surv_opt$Mod, levels =c(1,2,3,4,5,6,7,8))
surv_opt_all <- rbind(surv_opt, LRM_surv_opt, GEE_surv_opt, RRS.TA.surv)

Cindcomp = ggplot(data=surv_opt_all,
                   aes(x = LT,y =Cstat, ymin =CstatL, ymax =CstatU))+
  geom_errorbar(aes(ymin=CstatL, ymax=CstatU, col=Mod),width=0.2, cex=1, position = position_dodge(0.25))+ 
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism-adjusted \n Harrell's C-index (95% CI)")+
  scale_colour_manual(name=NULL, values=c(cols_models, "black"), labels=c("LRM", "GEE", "FPSM","SRF", "TDCM", "LM", "JM","RRS")) +theme_bw() +
  theme(text = element_text(size=14))
Cindcomp


BsComp<-ggplot(data=surv_opt_all,
                aes(x = LT,y =BS))+
  geom_point(aes(col = Mod), position = position_dodge(0.25))+ geom_line(aes(col = Mod), position = position_dodge(0.25)) +
  scale_colour_manual(name=NULL, values=c(cols_models, "black"), labels=c("LRM", "GEE", "FPSM","SRF", "TDCM", "LM", "JM","RRS")) + theme_bw() +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism-adjusted \n Brier score")+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(aspect.ratio=1, text = element_text(size=14), legend.position = "none")
BsComp


ipacomp<-ggplot(data=surv_opt_all,
                 aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA, colour=Mod), position = position_dodge(0.25))+
  geom_point(aes(col=Mod), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("Optimism-adjusted \n IPA (%)")+
  geom_hline(yintercept =0, linetype=2)+
  scale_colour_manual(name=NULL, values=c(cols_models, "black"), labels=c("LRM", "GEE", "FPSM","SRF", "TDCM", "LM", "JM","RRS")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipacomp

###  Apparent
Cindcomp2 = ggplot(data=Survival_Comp_all,
                   aes(x = LT,y =Cstat, ymin =CstatL, ymax =CstatU))+
  geom_errorbar(aes(ymin=CstatL, ymax=CstatU, col=Mod),width=0.2, cex=1, position = position_dodge(0.25))+ 
  geom_pointrange(aes(col=Mod), position = position_dodge(0.25))+
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Harrell's C-index (95% CI)")+
  scale_colour_manual(name=NULL, values=c(cols_models, "black"), labels=c("LRM", "GEE", "FPSM","SRF", "TDCM", "LM", "JM","RRS")) +theme_bw() +
  theme( text = element_text(size=14))
Cindcomp2


BsComp2<-ggplot(data=Survival_Comp_all,
                aes(x = LT,y =BS))+
  geom_point(aes(col = Mod), position = position_dodge(0.25))+ geom_line(aes(col = Mod), position = position_dodge(0.25)) +
  scale_colour_manual(name=NULL, values=c(cols_models, "black"), labels=c("LRM", "GEE", "FPSM","SRF", "TDCM", "LM", "JM","RRS")) + theme_bw() +
  xlab('Prediction time \n (years post-baseline)')+ ylab("Brier score")+ 
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme(text = element_text(size=14), legend.position = "none")
BsComp2


ipacomp2<-ggplot(data=Survival_Comp_all,
                 aes(x = LT,y =IPA))+
  geom_linerange(aes(x = LT, ymin = 0, ymax = IPA, colour=Mod), position = position_dodge(0.25))+
  geom_point(aes(col=Mod), size=2, position = position_dodge(0.25))+
  xlab('Prediction time \n (years post-baseline)')+ ylab("IPA (%)")+
  geom_hline(yintercept =0, linetype=2)+
  scale_colour_manual(name=NULL, values=c(cols_models, "black"), labels=c("LRM", "GEE", "FPSM","SRF", "TDCM", "LM", "JM","RRS")) +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  theme_bw() + theme(aspect.ratio=1, text = element_text(size=14))
ipacomp2

Cindcomp / (BsComp | ipacomp)
Cindcomp2 / (BsComp2 | ipacomp2)
