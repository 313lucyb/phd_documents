library(rstpm2)
library(parallel)

singlerun_sim <- function(dgm, dimension, probability_missing, n, j, weibull_shape, 
                        weibull_scale, validation_data, 
                        prediction_window, landmark_times, i){
  
  data <- simulate_data(dgm, dimension, probability_missing, n, j, weibull_shape, weibull_scale, i)
  survival_data <- dplyr::filter(data, time == 0)
  validation_survival_data <- validation_data %>% dplyr::filter(time == 0)
  print("singlerun_sim :: model fitting data simulated")
  
  events <- survival_data %>% filter(event == 1, time == 0)
  event_prop <- (nrow(events)/n)*100
  print(paste("Simulated dataset has sample size:", i))
  print(n)
  print(paste("Simulated dataset has event prevalence:", i))
  print(event_prop)
  print("%")
  
  ### adjust data structure for different models
  
  data_fpsm <- to_baseline(data)
  data_tdcm <- to_tdcm(data)
  data_lm_agg <- rbindlist(to_lm_agg(data, prediction_window, max_time = 4))
  jm_data<-data
  jm_survival_data<-survival_data
  
  #data_lm_long <- to_lm_long(data, prediction_window, max_time = 4)
  
  print(paste("singlerun_sim :: data transformed",i))
  ### model fitting
  
  #### flexi-parametric survival model
  start <- Sys.time()
  fpsm <- tryCatch(
    {
      stpm2(Surv(Time, event)~ age + sex + 
              disease_duration + clinician_disease + 
              patient_disease + bmi, data = data_fpsm)
    },
    error = function(error_message){
      print(error_message)
      return(NULL)
    }
  )
  end <- Sys.time()
  fpsm_fit_time <- as.numeric(end - start, units = "mins")
  print(paste("singlerun_sim :: fpsm fitted",i))
  #### time-dependent Cox model
  
  start <- Sys.time()
  tdcm <- tryCatch(
    {
      coxph(Surv(start, stop, event) ~ age + sex + 
              disease_duration + clinician_disease + 
              patient_disease + bmi, data = data_tdcm)
    },
    error = function(error_message){
      print(error_message)
      return(NULL)
    }
  )
  end <- Sys.time()
  tdcm_fit_time <- as.numeric(end - start, units = "mins")
  print(paste("singlerun_sim :: tdcm fitted",i))
  #### aggregated landmark model (LOCF)
  
  start <- Sys.time()
  lm_locf <- tryCatch(
    {
      coxph(Surv(LM, Time, event) ~ age + sex + 
              disease_duration + clinician_disease + 
              patient_disease + bmi + strata(LM) + 
              cluster(id), 
            data=data_lm_agg, method="breslow")
    },
    error = function(error_message){
      print(error_message)
      return(NULL)
    }
  )
  end <- Sys.time()  
  lm_locf_fit_time <- as.numeric(end - start, units = "mins")
  print(paste("singlerun_sim :: lm_locf fitted", i))
  #### models dependent on dimension (lm variations and jm)
  
  if (dimension == 1){
    
    start <- Sys.time()
    lm_median <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + median_clinician_disease + 
                patient_disease + bmi + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("lm median error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_median_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_median fitted", i))
    #### Aggregated landmark model (Min)
    
    start <- Sys.time()
    lm_min <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + min_clinician_disease + 
                patient_disease + bmi + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("lm min error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_min_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_min fitted", i))
    #### Aggregated landmark model (Max)
    
    start <- Sys.time()
    lm_max <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + max_clinician_disease + 
                patient_disease + bmi + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("lm max error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_max_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_max fitted", i))
    #### Aggregated landmark model (LOCF + SD)
    
    start <- Sys.time()
    lm_locf_sd <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + clinician_disease + 
                patient_disease + bmi +
                sd_clinician_disease + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("locf sd lm error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_locf_sd_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_locf_sd fitted",i))
    #### Joint model (current value)
    
    start<-Sys.time()
    survival_model <- tryCatch(
      {
        coxph(Surv(Time,event) ~ age + sex + 
                disease_duration + patient_disease + bmi, data = jm_survival_data)
      },
      error = function(error_message){
        print("baseline survival for joint model")
        print(error_message)
        return(NULL)
      }
    )
    print(paste("singlerun_sim :: survival model fitted for joint model", i))
    longitudinal_model_1 <- tryCatch(
      {
        lme(clinician_disease ~ sex*time, random = ~ time | id, data = jm_data, control = list(opt = "optim"), na.action = na.exclude)
      },
      error = function(error_message){
        print("clinician disease longitudinal model")
        print(error_message)
        return(NULL)
      }
    )
    print(paste("singlerun_sim :: ME model fitted for joint model", i))
    if (!is.null(survival_model) & !is.null(longitudinal_model_1)){
      joint_model <- tryCatch(
        {
          jm(survival_model, list(longitudinal_model_1), time_var = "time", id_var = "id", cores = 1)
        },
        error = function(error_message){
          print("joint model")
          print(error_message)
          return(NULL)
        })
    } else {
      joint_model <- NULL
    }
    end <- Sys.time()
    jm_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: joint model fitted", i))
  }
  
  if (dimension == 2){
    
    start <- Sys.time()
    lm_median <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + median_clinician_disease + 
                median_patient_disease + bmi + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("lm median error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_median_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_median fitted", i))
    #### Aggregated landmark model (Min)
    
    start <- Sys.time()
    lm_min <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + min_clinician_disease + 
                min_patient_disease + bmi + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("lm min error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_min_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_min fitted", i))
    #### Aggregated landmark model (Max)
    
    start <- Sys.time()
    lm_max <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + min_clinician_disease + 
                min_patient_disease + bmi + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("lm max error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_max_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("lm_locf", i))
    
    print(paste("singlerun_sim :: lm_max fitted", i))
    #### Aggregated landmark model (LOCF + SD)
    
    start <- Sys.time()
    lm_locf_sd <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + clinician_disease + 
                patient_disease + bmi + 
                sd_clinician_disease + sd_patient_disease + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("lm locf error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_locf_sd_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_locf_sd fitted", i))
    #### Joint model (current value)
    
    start<-Sys.time()
    survival_model <- tryCatch(
      {
        coxph(Surv(Time,event) ~ age + sex + 
                disease_duration + bmi, data = jm_survival_data)
      },
      error = function(error_message){
        print("baseline survival for joint model")
        print(error_message)
        return(NULL)
      }
    )
    print(paste("singlerun_sim :: survival model for joint model fitted", i))
    longitudinal_model_1 <- tryCatch(
      {
        lme(clinician_disease ~ sex*time, random = ~ time | id, data = jm_data, control = list(opt = "optim"), na.action = na.exclude)
      },
      error = function(error_message){
        print("clinician disease longitudinal model")
        print(error_message)
        return(NULL)
      }
    )
    print(paste("singlerun_sim :: ME 1 for joint model fitted", i))
    longitudinal_model_2 <- tryCatch(
      {
        lme(patient_disease ~ sex*time, random = ~ time | id, data = jm_data, control = list(opt = "optim"), na.action = na.exclude)
      },
      error = function(error_message){
        print("patient disease longitudinal model")
        print(error_message)
        return(NULL)
      }
    )
    print(paste("singlerun_sim :: ME model 2 for joint model fitted", i))
    if (!is.null(survival_model) & !is.null(longitudinal_model_1) & !is.null(longitudinal_model_2)){

      joint_model <- tryCatch(
        {
          jm(survival_model, list(longitudinal_model_1,longitudinal_model_2), time_var = "time", id_var = "id", cores = 1)
        },
        error = function(error_message){
          print("joint model")
          print(error_message)
          return(NULL)
        })
    } else {
      joint_model <- NULL
    }
    end <- Sys.time()
    jm_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: joint model fitted", i))
  }
  
  if (dimension == 3){
    #### Aggregated landmark model (median)
    
    start <- Sys.time()
    lm_median <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + median_clinician_disease + 
                median_patient_disease + median_bmi + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print(error_message)
        print("lm median error")
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_median_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_median fitted", i))
    #### Aggregated landmark model (Min)
    
    start <- Sys.time()
    lm_min <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + min_clinician_disease +
                min_patient_disease + min_bmi + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("lm min error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_min_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_min fitted", i))
    #### Aggregated landmark model (Max)
    
    start <- Sys.time()
    lm_max <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + max_clinician_disease +
                max_patient_disease + max_bmi + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("lm max error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_max_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_max fitted", i))
    #### Aggregated landmark model (LOCF + SD)
    
    start <- Sys.time()
    lm_locf_sd <- tryCatch(
      {
        coxph(Surv(LM, Time, event) ~ age + sex + 
                disease_duration + clinician_disease +
                patient_disease + bmi + sd_clinician_disease +
                sd_patient_disease + sd_bmi + strata(LM) + 
                cluster(id), 
              data=data_lm_agg, method="breslow")
      },
      error = function(error_message){
        print("lm sd error")
        print(error_message)
        return(NULL)
      }
    )
    end <- Sys.time()
    lm_locf_sd_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: lm_locf_sd fitted", i))
    #### Joint model (current value)
    
    start<-Sys.time()
    survival_model <- tryCatch(
      {
        coxph(Surv(Time,event) ~ age + sex + disease_duration, data = jm_survival_data)
      },
      error = function(error_message){
        print("baseline survival for joint model")
        print(error_message)
        return(NULL)
      }
    )
    print(paste("singlerun_sim :: survival model for joint model fitted", i))
    longitudinal_model_1 <- tryCatch(
      {
        lme(clinician_disease ~ sex*time, random = ~ time | id, data = jm_data, control = list(opt = "optim"), na.action = na.exclude)
      },
      error = function(error_message){
        print("clinician disease longitudinal model")
        print(error_message)
        return(NULL)
      }
    )
    print(paste("singlerun_sim :: ME 1 for joint model fitted", i))
    longitudinal_model_2 <- tryCatch(
      {
        lme(patient_disease ~ sex*time, random = ~ time | id, data = jm_data, control = list(opt = "optim"), na.action = na.exclude)
      },
      error = function(error_message){
        print("patient disease longitudinal model")
        print(error_message)
        return(NULL)
      }
    )
    print(paste("singlerun_sim :: ME 2 for joint model fitted", i))
    longitudinal_model_3 <- tryCatch(
      {
        lme(bmi ~ sex*time, random = ~ time | id, data = data, control = list(opt = "optim"), na.action = na.exclude)
      },
      error = function(error_message){
        print("bmi longitudinal model")
        print(error_message)
        return(NULL)
      }
    )
    print(paste("singlerun_sim :: ME 3 for joint model fitted", i))
    if (!is.null(survival_model) & !is.null(longitudinal_model_1) & 
        !is.null(longitudinal_model_2) & !is.null(longitudinal_model_3)){
      joint_model <- tryCatch(
        {
          jm(survival_model, list(longitudinal_model_1, 
                                  longitudinal_model_2, 
                                  longitudinal_model_3), time_var = "time", id_var = "id", cores = 1)
        },
        error = function(error_message){
          print("joint model")
          print(error_message)
          return(NULL)
        })
    } else {
      joint_model <- NULL
    }
    end <- Sys.time()
    jm_fit_time <- as.numeric(end - start, units = "mins")
    print(paste("singlerun_sim :: joint model fitted", i))
  }
  
  ### Generate temporal assessment data
  val_dev_data <- prediction_datasets(data, prediction_window, landmark_times, type = "short")
  val_new_data <- prediction_datasets(validation_data, prediction_window, landmark_times, type = "short")
  print(paste("singlerun_sim :: prediction datasets created", i))
  ### model predictions and performance metrics
  
  if(!is.null(fpsm)) {
    if (!any(is.na(as.numeric(summary(fpsm)@coef[2:7,1])))){
      fpsm_perf <- model_output(model = fpsm, fit_time = fpsm_fit_time, 
                                val_dev_data, val_new_data, 
                                pred_var = "fpsm_pred", 
                                lp_var = "fpsm_lp", model_id = 1, iter = i, landmark_times,
                                prediction_window, model_type = "fpsm")
      fpsm_output <- fpsm_perf[[1]]
      fpsm_stability <- fpsm_perf[[2]]
    }
    else{
      fpsm_output <- matrix(NA, nrow=1, ncol=6)
      fpsm_output[1,]<- c(i, 1, NA, NA, "failed_fit", 1)
      fpsm_stability<- matrix(NA, nrow=1, ncol=6)
      fpsm_stability[1,]<- c(i, NA, 1, NA, NA, NA)
      fpsm_stability<-data.frame(fpsm_stability)
      names(fpsm_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                "prediction_window", "prediction")
    }
  } else {
    fpsm_output <- matrix(NA, nrow=1, ncol=6)
    fpsm_output[1,]<- c(i, 1, NA, NA, "failed_fit", 1)
    fpsm_stability<- matrix(NA, nrow=1, ncol=6)
    fpsm_stability[1,]<- c(i, NA, 1, NA, NA, NA)
    fpsm_stability<-data.frame(fpsm_stability)
    names(fpsm_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                              "prediction_window", "prediction")
    }
  print(paste("singlerun_sim :: fpsm output generated", i))
  
  if(!is.null(tdcm)){
    if (!any(is.na(tdcm$coefficients))){
      tdcm_perf <- model_output(model = tdcm, fit_time = tdcm_fit_time, 
                                val_dev_data, val_new_data, 
                                pred_var = "tdcm_pred", 
                                lp_var = "tdcm_lp", model_id = 2, iter = i, landmark_times,
                                prediction_window, model_type = "tdcm")
      tdcm_output <- tdcm_perf[[1]]
      tdcm_stability <- tdcm_perf[[2]]
    } else {
        tdcm_output <- matrix(NA, nrow=1, ncol=6)
        tdcm_output[1,]<- c(i, 2, NA, NA, "failed_fit", 1)
        tdcm_stability<- matrix(NA, nrow=1, ncol=6)
        tdcm_stability[1,]<- c(i, NA, 2, NA, NA, NA)
        tdcm_stability<-data.frame(tdcm_stability)
        names(tdcm_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                  "prediction_window", "prediction")
    }
  } else {
    tdcm_output <- matrix(NA, nrow=1, ncol=6)
    tdcm_output[1,]<- c(i, 2, NA, NA, "failed_fit", 1)
    tdcm_stability<- matrix(NA, nrow=1, ncol=6)
    tdcm_stability[1,]<- c(i, NA, 2, NA, NA, NA)
    tdcm_stability<-data.frame(tdcm_stability)
    names(tdcm_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                 "prediction_window", "prediction")
  }
  print(paste("singlerun_sim :: tdcm output generated", i))
  
  if (!is.null(lm_locf)){
    if (!any(is.na(lm_locf$coefficients))){
      lm <- lm_model_output(model = lm_locf, fit_time = lm_locf_fit_time, 
                            development_data = data_lm_agg, 
                            val_dev_data, val_new_data, pred_var = "lm_pred", 
                            lp_var = "lm_lp", model_id = 3, iter = i, landmark_times,
                            prediction_window)
      lm_locf_output <- lm[[1]]
      lm_locf_stability <- lm[[2]]
    } else {
      lm_locf_output <- matrix(NA, nrow=1, ncol=6)
      lm_locf_output[1,]<- c(i, 3, NA, NA, "failed_fit", 1)
      
      lm_locf_stability<- matrix(NA, nrow=1, ncol=6)
      lm_locf_stability[1,]<- c(i, NA, 3, NA, NA, NA)
      lm_locf_stability<-data.frame(lm_locf_stability)
      names(lm_locf_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                   "prediction_window", "prediction")
    }
  } else {
    lm_locf_output <- matrix(NA, nrow=1, ncol=6)
    lm_locf_output[1,]<- c(i, 3, NA, NA, "failed_fit", 1)
    
    lm_locf_stability<- matrix(NA, nrow=1, ncol=6)
    lm_locf_stability[1,]<- c(i, NA, 3, NA, NA, NA)
    lm_locf_stability<-data.frame(lm_locf_stability)
    names(lm_locf_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                "prediction_window", "prediction")
  }
  print(paste("singlerun_sim :: lm_locf output generated", i))
  
  if (!is.null(lm_median)){
    if (!any(is.na(lm_median$coefficients))){
    
    lm <- lm_model_output(model = lm_median, fit_time = lm_median_fit_time, 
                          development_data = data_lm_agg, 
                          val_dev_data, val_new_data, pred_var = "lm_pred", 
                          lp_var = "lm_lp", model_id = 4, iter = i, landmark_times,
                          prediction_window)
    lm_median_output <- lm[[1]]
    lm_median_stability <- lm[[2]]
    } else {
      lm_median_output <- matrix(NA, nrow=1, ncol=6)
      lm_median_output[1,]<- c(i, 4, NA, NA, "failed_fit", 1)
      
      lm_median_stability<- matrix(NA, nrow=1, ncol=6)
      lm_median_stability[1,]<- c(i, NA, 4, NA, NA, NA)
      lm_median_stability<-data.frame(lm_median_stability)
      names(lm_median_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                     "prediction_window", "prediction")
    }
  } else {
    lm_median_output <- matrix(NA, nrow=1, ncol=6)
    lm_median_output[1,]<- c(i, 4, NA, NA, "failed_fit", 1)
    
    lm_median_stability<- matrix(NA, nrow=1, ncol=6)
    lm_median_stability[1,]<- c(i, NA, 4, NA, NA, NA)
    lm_median_stability<-data.frame(lm_median_stability)
    names(lm_median_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                "prediction_window", "prediction")
  }
  print(paste("singlerun_sim :: lm_median output generated", i))
  
  if (!is.null(lm_min)){
    if (!any(is.na(lm_min$coefficients))){
      lm <- lm_model_output(model = lm_min, fit_time = lm_min_fit_time, 
                            development_data = data_lm_agg, 
                            val_dev_data, val_new_data, pred_var = "lm_pred", 
                            lp_var = "lm_lp", model_id = 5, iter = i, landmark_times,
                            prediction_window)
      lm_min_output <- lm[[1]]
      lm_min_stability <- lm[[2]]
    } else {
      lm_min_output <- matrix(NA, nrow=1, ncol=6)
      lm_min_output[1,]<- c(i, 5, NA, NA, "failed_fit", 1)
      
      lm_min_stability<- matrix(NA, nrow=1, ncol=6)
      lm_min_stability[1,]<- c(i, NA, 5, NA, NA, NA)
      lm_min_stability<-data.frame(lm_min_stability)
      names(lm_min_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                  "prediction_window", "prediction")
    }
  } else {
    lm_min_output <- matrix(NA, nrow=1, ncol=6)
    lm_min_output[1,]<- c(i, 5, NA, NA, "failed_fit", 1)
    
    lm_min_stability<- matrix(NA, nrow=1, ncol=6)
    lm_min_stability[1,]<- c(i, NA, 5, NA, NA, NA)
    lm_min_stability<-data.frame(lm_min_stability)
    names(lm_min_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                "prediction_window", "prediction")
  }
  print(paste("singlerun_sim :: lm_min output generated", i))
  if (!is.null(lm_max)){
    if (!any(is.na(lm_max$coefficients))){
    
    lm <- lm_model_output(model = lm_max, fit_time = lm_max_fit_time, 
                          development_data = data_lm_agg, 
                          val_dev_data, val_new_data, pred_var = "lm_pred", 
                          lp_var = "lm_lp", model_id = 6, iter = i, landmark_times,
                          prediction_window)
    lm_max_output <- lm[[1]]
    lm_max_stability <- lm[[2]]
    } else {
      lm_max_output <- matrix(NA, nrow=1, ncol=6)
      lm_max_output[1,]<- c(i, 6, NA, NA, "failed_fit", 1)
      
      lm_max_stability<- matrix(NA, nrow=1, ncol=6)
      lm_max_stability[1,]<- c(i, NA, 6, NA, NA, NA)
      lm_max_stability<-data.frame(lm_max_stability)
      names(lm_max_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                  "prediction_window", "prediction")
    }
  } else {
    lm_max_output <- matrix(NA, nrow=1, ncol=6)
    lm_max_output[1,]<- c(i, 6, NA, NA, "failed_fit", 1)
    
    lm_max_stability<- matrix(NA, nrow=1, ncol=6)
    lm_max_stability[1,]<- c(i, NA, 6, NA, NA, NA)
    lm_max_stability<-data.frame(lm_max_stability)
    names(lm_max_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                    "prediction_window", "prediction")
  }
  print(paste("singlerun_sim :: lm_max output generated", i))
  
  if (!is.null(lm_locf_sd)){
    if (!any(is.na(lm_locf_sd$coefficients))){
      lm <- lm_model_output(model = lm_locf_sd, fit_time = lm_locf_sd_fit_time, 
                            development_data = data_lm_agg, 
                            val_dev_data, val_new_data, pred_var = "lm_pred", 
                            lp_var = "lm_lp", model_id = 7, iter = i, landmark_times,
                            prediction_window)
      lm_locf_sd_output <- lm[[1]]
      lm_locf_sd_stability <- lm[[2]]
    } else {
      lm_locf_sd_output <- matrix(NA, nrow=1, ncol=6)
      lm_locf_sd_output[1,]<- c(i, 7, NA, NA, "failed_fit", 1)
      
      lm_locf_sd_stability<- matrix(NA, nrow=1, ncol=6)
      lm_locf_sd_stability[1,]<- c(i, NA, 7, NA, NA, NA)
      lm_locf_sd_stability<-data.frame(lm_locf_sd_stability)
      names(lm_locf_sd_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                                      "prediction_window", "prediction")
    }
  } else {
    lm_locf_sd_output <- matrix(NA, nrow=1, ncol=6)
    lm_locf_sd_output[1,]<- c(i, 7, NA, NA, "failed_fit", 1)
    
    lm_locf_sd_stability<- matrix(NA, nrow=1, ncol=6)
    lm_locf_sd_stability[1,]<- c(i, NA, 7, NA, NA, NA)
    lm_locf_sd_stability<-data.frame(lm_locf_sd_stability)
    names(lm_locf_sd_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                            "prediction_window", "prediction")
  }
  print(paste("singlerun_sim :: lm_locf_sd output generated", i))

  if (!is.null(joint_model)){
    
    ## extract model fitting times
    
    fit_times<-matrix(NA, nrow=1, ncol=6)
    fit_times[1,] <- c(i, 8, NA, NA, "fit_time", jm_fit_time)
    
    ## generate predictions and prediction times from model
    
    val_dev_JM_output <- predictions_jm(joint_model, data, survival_data, landmark_times, prediction_window)
    val_dev_preds_JM <- val_dev_JM_output[1:length(landmark_times)]
    val_new_JM_output <- predictions_jm(joint_model, validation_data, 
                                        validation_survival_data, 
                                        landmark_times, prediction_window)
    val_new_preds_JM <- val_new_JM_output[1:length(landmark_times)]
    val_new_preds_time_JM <- val_new_JM_output[[(length(landmark_times)+1)]]
    
    ## put prediction times into consistent structure across models
    
    pred_times<-prediction_time_matrix(landmark_times, val_new_preds_time_JM, 8, prediction_window, i)
    
    ## extract predictions on validation data by model to assess stability
    
    stability <- stability_data(val_new_preds_JM, "pred", landmark_times, 8, i)
    stability$PW<-rep(prediction_window, length(stability$SubjID))
    stability <- stability %>%
      dplyr::select(SimID, SubjID, Mod, LT, PW, Pred)
    stability<-as.data.frame(stability)
    names(stability)<-c("simulation_id", "subject_id", "model", "landmark_time", 
                        "prediction_window", "prediction")
    stability<-as.data.frame(stability)
    
    ## predictive performance on validation and development data

    validation_data<-performance_assessment(val_new_preds_JM, "pred", "pred", "Time", "event", 
                                            landmark_times, prediction_window, i, 8)
    apparent_data<-performance_assessment(val_dev_preds_JM, "pred", "pred", "Time", "event", 
                                          landmark_times, prediction_window, i, 8)
    optimism_values <-as.numeric(apparent_data[,6])-as.numeric(validation_data[,6])
    optimism<-matrix(NA, nrow = length(optimism_values), ncol = 6)
    optimism[,1]<-rep(i, times = length(optimism_values))
    optimism[,2]<-rep(8, length(optimism_values))
    optimism[,3]<-c(rep(landmark_times, each = 4))
    optimism[,4]<-rep(prediction_window, times = length(optimism_values))
    optimism[,5]<-rep(c("brier_adjusted", "c_index_adjusted", 
                        "c_index_lower_adjusted", "c_index_higher_adjusted"), times=length(landmark_times))
    optimism[,6]<-optimism_values
    
    jm_output<- rbind(fit_times, pred_times, validation_data, optimism)
    jm_stability<-stability
    
  } else {
    jm_output <- matrix(NA, nrow=1, ncol=6)
    jm_output[1,]<- c(i, 8, NA, NA, "failed_fit", 1)
    
    jm_stability<- matrix(NA, nrow=1, ncol=6)
    jm_stability[1,]<- c(i, NA, 8, NA, NA, NA)
    jm_stability<-data.frame(jm_stability)
    names(jm_stability)<- c("simulation_id", "subject_id", "model", "landmark_time", 
                            "prediction_window", "prediction")
  }
  print(paste("singlerun_sim :: jm output generated", i))
  
  simulations<-rbind(fpsm_stability, 
                     tdcm_stability, 
                     lm_locf_stability,
                     lm_median_stability,
                     lm_min_stability,
                     lm_max_stability,
                     lm_locf_sd_stability,
                     jm_stability)
  names(simulations)<-c("simulation_id", "subject_id", "model", "landmark_time", 
                        "prediction_window", "prediction")
  print(paste("singlerun_sim :: stability data merged", i))
  
  ## null brier score
  
  null_briers = list()
  for (k in 1:length(val_new_data)){
    dataset = val_new_data[[k]]
    null_brier_val = as.numeric(null_brier(data = dataset, 
                            btime = landmark_times[k]+prediction_window)[1])
    null_briers <- c(null_briers, null_brier_val)
  }
  null_brier_mat = matrix(NA, nrow=length(landmark_times), ncol=6)
  null_brier_mat[,1]<-rep(i, times = length(landmark_times))
  null_brier_mat[,2]<-rep(NA, times = length(landmark_times))
  null_brier_mat[,3]<-landmark_times
  null_brier_mat[,4]<-rep(prediction_window, times = length(landmark_times))
  null_brier_mat[,5]<-rep("null_brier", times = length(landmark_times))
  null_brier_mat[,6]<-unlist(null_briers)
  print(paste("singlerun_sim :: null briers generated", i))
  ### output reporting
  
  results<-rbind(fpsm_output, 
                 tdcm_output, 
                 lm_locf_output,
                 lm_median_output,
                 lm_min_output,
                 lm_max_output,
                 lm_locf_sd_output,
                 jm_output,
                 null_brier_mat)
  results<-as.data.frame(results)
  
  names(results)<-c("simulation_id", "model_id", "landmark_time", 
                    "prediction_window", "statistic", "value")
  
  print(paste("singlerun_sim :: performance data merged", i))
  
  output<-list(results, simulations)
  
  return(output)
}

across_sim_func_parallel<-function(dgm, n_simulation, 
                                   validation_sample_size, dimension,
                                   n, j, weibull_shape, weibull_scale, probability_missing,
                                   prediction_window, landmark_times, scenario){
  
  validation_data <- simulate_data(dgm, dimension, probability_missing, 
                                   n = validation_sample_size, j, 
                                   weibull_shape, weibull_scale, i = 454 + scenario)
  val_events <- validation_data %>% filter(time == 0, event == 1)
  val_ev_prop <- (nrow(val_events)/validation_sample_size)*100
  print("Validation dataset has event prevalence of:")
  print(val_ev_prop)
  
  print("across_sim_func_parallel :: validation data simulated")
  numCoresAllowed <- as.integer(Sys.getenv("NSLOTS", unset=1))
  cl<-makeCluster(4, outfile="") 
  registerDoParallel(cl)
  output <- foreach(i=1:n_simulation, .verbose=TRUE, .errorhandling="remove", .combine=rbind, .packages = c("ipred", "MASS", "dynpred", "JMbayes2",
                                                                    "rms","survival", "Hmisc", "dplyr", "data.table", "rstpm2")) %dopar% {
                                                                      print("test print")
                                                                      set.seed(656543 + scenario + i)
                                                                      source("./data_generation_functions_new.R", local=TRUE)
                                                                      source("./prediction_functions_new.R", local=TRUE)
                                                                      source("./performance_functions.R", local=TRUE)
                                                                      source("./simulation_functions_new.R", local=TRUE)
                                                                      singlerun_sim(dgm, dimension, 
                                                                                    probability_missing, 
                                                                                    n, j, weibull_shape, 
                                                                                    weibull_scale, 
                                                                                    validation_data, 
                                                                                    prediction_window, 
                                                                                    landmark_times, i)
                                                                    }
  stopCluster(cl)
  print("across_sim_func_parallel :: simulations performed for scenario")
  results_output<-as.data.frame(rbindlist(output[,1]))
  print("across_sim_func_parallel :: results_output")
  stability_output<-as.data.frame(rbindlist(output[,2]))
  print("across_sim_func_parallel :: output")
  names(stability_output)<-c("simulation_id", "subject_id", "model", "landmark_time", 
                             "prediction_window", "prediction")
  print("across_sim_func_parallel :: names")
  stability<- stability_output %>%
    dplyr::group_by(subject_id, model, landmark_time, prediction_window) %>%
    summarise(range=quantile(prediction, 0.95, na.rm = TRUE)-quantile(prediction, 0.05, na.rm = TRUE))
  stability<-as.data.frame(stability)
  print("across_sim_func_parallel :: stability summary")
  output<-list(results_output, stability)
  return(output)
}



