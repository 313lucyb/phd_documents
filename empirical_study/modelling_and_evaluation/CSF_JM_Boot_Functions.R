JMperf<-function(LT, HT, obj, long_data, i, long_pred, seed){
  
  JM_val_new <- long_data %>% filter(pyears <= LT & LMvadtime > LT) %>% drop_na(long_pred)
  JM_val_new <- JM_val_new %>% group_by(groupid) %>% mutate(bmi = na.locf(bmi))
  JM_ev_new <- JM_val_new %>% group_by(groupid) %>% top_n(1, pyears) %>%
    mutate(val_event = case_when(timevent < HT ~ 1, TRUE ~ 0)) %>%
    mutate(val_timevent = case_when(timevent < HT ~ timevent, TRUE ~ HT)) %>%
    select(groupid, val_event, val_timevent)
  set.seed(seed)
  start_time <- Sys.time()
  print(head(JM_val_new))
  s <- tryCatch(
    {
      survfitJM(object = obj, newdata = JM_val_new, 
                idVar = "groupid", survTimes = HT, last.time = LT, 
                type = "SurvProb", simulate = TRUE)
    },
    error=function(cond){
      print(cond)
      return(NA)
    }
  )
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
  perf<-c(1,LT,c,cL,cU, BS, pred_time)
  print(perf)
  return(perf)
}

# TODO: make sure the long_data have min event times, censoring at 4 years and LMvadtime

bootvadJM<-function(long_data, event_data, B, seed, s){
  
  dat <- long_data
  sdat <- event_data
  indiv <- unique(dat$groupid)
  output <-list()  
  
  for (j in 1:B){
    
    # generate bootstrap dataset, and renaming duplicate IDs
    
    seed<-seed+j
    set.seed(seed)
    
    smp <- sort(sample(indiv, length(indiv), replace=TRUE))
    smp.df <- data.frame(groupid=smp)
    smp.df<-smp.df %>% group_by(groupid) %>% mutate(count = row_number())
    smp.df$newid<-paste(smp.df$groupid, smp.df$count, sep=".")
    
    boot_data = merge(smp.df, dat, by = "groupid", all.x=TRUE)
    boot_surv_data = merge(smp.df, sdat, by = "groupid", all.x=TRUE)
    boot_data$groupid <- boot_data$newid
    boot_data<-boot_data[with(boot_data, order(groupid, pyears)),]
    boot_surv_data$groupid <- boot_surv_data$newid
    boot_surv_data<-boot_surv_data[with(boot_surv_data, order(groupid, pyears)),]
    # model fit to bootstrap data
    
    lmeFit<- lme(ovmean~ns(pyears, 3), data = boot_data, 
                 random = ~ ns(pyears, 2) | groupid, control = list(opt = "optim"))
    print(summary(lmeFit))
    survFit <- coxph(Surv(timevent, event) ~ age + pgen + as.factor(trtment) + 
                       disdur + previous_dmards + lung + 
                       diabetes + bmi + renal + steroids + dascore + 
                       as.factor(smoke), data = boot_surv_data, x = TRUE)
    print(summary(survFit))
    jointFit <- jointModelBayes(lmeFit, survFit, timeVar = "pyears")
    
    
    LT<-c(0,0.5,1,1.5,2,2.5,3)
    HT<-LT+1
    
    cl<-makeCluster(cores[1]-1) #always useful to keep one core free
    registerDoParallel(cl)
    JMTA_apparent<- foreach(i=1:7, .combine=rbind, .packages = c("ipred", "splines", 
                                                                 "survival", "JMbayes", "dplyr")) %do% {
                                                                   JMperf(LT[i],HT[i], jointFit, boot_data, i, "ovmean", seed)
                                                                 }
    stopCluster(cl)
    
    cl<-makeCluster(cores[1]-1) #always useful to keep one core free
    registerDoParallel(cl)
    JMTA_original<- foreach(i=1:7, .combine=rbind, .packages = c("ipred", "splines", 
                                                                 "survival", "JMbayes", "dplyr")) %do% {
                                                                   JMperf(LT[i],HT[i], jointFit, long_data, i, "ovmean", seed)
                                                                 }
    stopCluster(cl)
    optimism <- matrix(NA, nrow = 7, ncol=6)
    optimism[,1] <- rep(1, time=7)
    optimism[,2] <- LT
    optimism[, (3:6)] = JMTA_apparent[, (3:6)] - JMTA_original[,(3:6)]
    optimism < -as.data.frame(optimism)
    names(optimism)<- c("mod", "LT", "C", "C_low", "C_high", "Brier")
    output[[j]] <- optimism
  }
  
  output<-as.data.frame(do.call(rbind, output))
  write.table(output, file = paste("./JM_optimism_results_",s, ".csv", sep=""), sep = ",", col.names = NA,
              qmethod = "double")
  print(output)
  return(output)
}
