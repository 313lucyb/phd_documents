library(survival)
library(JMbayes)
library(splines)
library(ipred)
library(dplyr)
library(foreign)
library(foreach)
library(doParallel)
library(tidyr)
library(zoo)

long_data<-read.dta("val_long_JM.dta")
long_data$smoke<-factor(long_data$smoke)
long_data$trtment<-factor(long_data$trtment)
long_data <- long_data %>% group_by(groupid) %>% mutate(dascore = na.locf(dascore))

# fitted model

jointFit <- readRDS("~/jointmodelresults/haq1.rds")

file.names<-c("cal10.csv","cal11.csv","cal12.csv","cal13.csv", "cal14.csv", "cal15.csv", "cal16.csv")

JMperf<-function(LT, HT, obj, long_data, i, long_pred, seed){
  
  JM_val_new <- long_data %>% filter(pyears <= LT & LMvadtime > LT) %>% drop_na(long_pred)
  JM_ev_new <- JM_val_new %>% group_by(groupid) %>% top_n(1, pyears) %>%
    mutate(val_event = case_when(timevent < HT ~ 1, TRUE ~ 0)) %>%
    mutate(val_timevent = case_when(timevent < HT ~ timevent, TRUE ~ HT)) %>%
    select(groupid, val_event, val_timevent)
  set.seed(seed)
  start_time <- Sys.time()
  s <- survfitJM(object = obj, newdata = JM_val_new, 
                 idVar = "groupid", survTimes = HT, last.time = LT, 
                 type = "SurvProb", simulate = TRUE)
  end_time <- Sys.time()
  pred_time <- end_time - start_time
  output<-matrix(unlist(s[[1]]), ncol = 5, byrow = TRUE)
  JM_ev_new$pred<-output[,2]
  m<-coxph(Surv(val_timevent,val_event)~pred, data=JM_ev_new)
  c<-summary(m)$concordance[1]
  cL<-summary(m)$concordance[1]-(1.96*summary(m)$concordance[2])
  cU<-summary(m)$concordance[1]+(1.96*summary(m)$concordance[2])
  BS <- sbrier(obj= Surv(JM_ev_new$val_timevent, JM_ev_new$val_event), 
               pred = JM_ev_new$pred, btime = HT)
  
  caldata<-cbind(JM_ev_new$pred, JM_ev_new$val_timevent, JM_ev_new$val_event)
  caldata<-as.data.frame(caldata)
  names(caldata)<-c("pred", "survtime", "event")
  write.table(caldata, file = file.names[i], sep= ",", col.names= NA, qmethod = "double")
  perf<-c(1,LT,c,cL,cU, BS, pred_time)
  print(perf)
  return(perf)
}

LT<-c(0,0.5,1,1.5,2,2.5,3)
HT<-LT+1

cl<-makeCluster(7, outfile="") 
registerDoParallel(cl)
JMTA<- foreach(i=1:7, .combine=rbind, .packages = c("ipred", "splines", 
                                                    "survival", "JMbayes", "dplyr", "tidyr")) %dopar% {
                                                      JMperf(LT[i], HT[i], jointFit, long_data, i, "ovmean", 310394+i)
                                                    }
stopCluster(cl)
print(JMTA)
write.table(JMTA, file = "./JM1_TA.csv", sep = ",", col.names = NA,
            qmethod = "double")
