library(foreign)
library(dplyr)
library(Hmisc)
library(MASS)
library(randomForestSRC)
library(ipred)
library(pec)

bootvadSRF<-function(data, B, seed, test0, test1, test2, test3, test4, test5, test6){
  dat <- data
  dat$studyno <- dat$groupid
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
    
    # tuning and fitting SRFs
    
    rsf.tune<-tune(formula = Surv(timevent, event) ~ age+pgen+dascore+ovmean
                   +firsttreat+bmi+as.factor(smoke)+as.factor(steroids)+lung+renal+
                     diabetes+disdur+previous_dmards, data = boot, 
                   ntreeTry = 500, nodesizeTry = seq(10, 300, by = 10), 
                   nsplit = 5, mtryStart = 3, trace = TRUE, doBest = TRUE)
    
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

SRF <- read.dta("baseline_model2_CSF.dta")
SRF$smoke<-factor(SRF$smoke)
SRF$firsttreat<-factor(SRF$firsttreat)
SRF$pgen<-factor(SRF$pgen)
SRF$renal<-factor(SRF$renal)
SRF$lung<-factor(SRF$lung)
SRF$diabetes<-factor(SRF$diabetes)
SRF$steroids<-factor(SRF$steroids)

test0<-SRF

test1<-read.dta("validation_6months_CSF.dta")
test1$smoke<-factor(test1$smoke)
test1$firsttreat<-factor(test1$firsttreat)
test1$pgen<-factor(test1$pgen)
test1$renal<-factor(test1$renal)
test1$lung<-factor(test1$lung)
test1$diabetes<-factor(test1$diabetes)
test1$steroids<-factor(test1$steroids)

test2<-read.dta("validation_12months_CSF.dta")
test2$smoke<-factor(test2$smoke)
test2$firsttreat<-factor(test2$firsttreat)
test2$pgen<-factor(test2$pgen)
test2$renal<-factor(test2$renal)
test2$lung<-factor(test2$lung)
test2$diabetes<-factor(test2$diabetes)
test2$steroids<-factor(test2$steroids)

test3<-read.dta("validation_18months_CSF.dta")
test3$smoke<-factor(test3$smoke)
test3$firsttreat<-factor(test3$firsttreat)
test3$pgen<-factor(test3$pgen)
test3$renal<-factor(test3$renal)
test3$lung<-factor(test3$lung)
test3$diabetes<-factor(test3$diabetes)
test3$steroids<-factor(test3$steroids)

test3<-read.dta("validation_24months_CSF.dta")
test4$smoke<-factor(test4$smoke)
test4$firsttreat<-factor(test4$firsttreat)
test4$pgen<-factor(test4$pgen)
test4$renal<-factor(test4$renal)
test4$lung<-factor(test4$lung)
test4$diabetes<-factor(test4$diabetes)
test4$steroids<-factor(test4$steroids)

test5<-read.dta("validation_30months_CSF.dta")
test5$smoke<-factor(test5$smoke)
test5$firsttreat<-factor(test5$firsttreat)
test5$pgen<-factor(test5$pgen)
test5$renal<-factor(test5$renal)
test5$lung<-factor(test5$lung)
test5$diabetes<-factor(test5$diabetes)
test5$steroids<-factor(test5$steroids)

test6<-read.dta("validation_36months_CSF.dta")
test6$smoke<-factor(test6$smoke)
test6$firsttreat<-factor(test6$firsttreat)
test6$pgen<-factor(test6$pgen)
test6$renal<-factor(test6$renal)
test6$lung<-factor(test6$lung)
test6$diabetes<-factor(test6$diabetes)
test6$steroids<-factor(test6$steroids)

cl<-makeCluster(20, outfile="") 
registerDoParallel(cl)
start_time <- Sys.time()
SRF.opt <- foreach(i=1:200, .combine=rbind, .packages = c("ipred", "MASS", "randomForestSRC", 
                                                          "pec", "Hmisc", "dplyr")) %dopar% {
                                                            seed<-199429+i
                                                            bootvadSRF(data = SRF, B=1, seed=seed, test0=test0, test1=test1, test2=test2, test3=test3, test4=test4, test5=test5, test6=test6)
                                                          }
stopCluster(cl)
end_time<-Sys.time()
bootcomptimeSRF <- end_time - start_time

print(SRF.opt)
print(dim(SRF.opt))
print(bootcomptimeSRF)

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
