# predictions from FPSM

predictions_fpsm<-function(prediction_datasets, model, landmark_times, covariates, PW){
  data = list()
  time = c()
  model_coefs<-as.numeric(summary(model)@coef[2:7,1])
  for (i in 1:length(prediction_datasets)){
    dataset = data.frame(prediction_datasets[[i]])
    lp_cov <- dataset %>% dplyr::select(all_of(covariates))
    lp_mat <- as.matrix(lp_cov)
    lp <- model_coefs %*% t(lp_mat)
    dataset$fpsm_lp <- as.vector(t(lp))
    print("predictions_fpsm :: linear predictor calculated")
    pred_dataset <- dataset
    pred_dataset$Time <- rep(landmark_times[i] + PW, times = nrow(dataset))
    start_time <- Sys.time()
    predictions <- rstpm2::predict(model, newdata = pred_dataset, type = "surv", se.fit=TRUE)
    print("predictions_fpsm :: predictions generated")
    dataset$fpsm_pred <- as.numeric(predictions$Estimate)
    print("predictions_fpsm :: predictions entered into dataset")
    end_time <- Sys.time()
    time_taken = as.numeric(end_time - start_time, units = "mins")
    time = c(time, time_taken)
    data[[i]] = dataset
  }
  data[[length(data)+1]] <- time
  return(data)
}

# predictions from TDCM

risk = function(model, newdata, time) {
  as.numeric(1-summary(survfit(model, newdata = newdata, 
                               se.fit = F, conf.int = F), 
                       times = time, extend = TRUE)$surv)
}

predictions_tdcm<-function(prediction_datasets, model, landmark_times, covariates, PW){
  data = list()
  time = c()
  model_coefs<-as.numeric(model$coefficients)
  for (i in 1:length(prediction_datasets)){
    dataset = data.frame(prediction_datasets[[i]])
    lp_cov <- dataset %>% dplyr::select(all_of(covariates))
    lp_mat <- as.matrix(lp_cov)
    lp <- model_coefs %*% t(lp_mat)
    dataset$tdcm_lp <- as.vector(t(lp))
    print("predictions_tdcm :: linear predictor generated")
    start_time <- Sys.time()
    tdcm_pred1 <- 1- risk(
      model, newdata = dataset, 
      time = landmark_times[i] + PW
    )
    print("predictions_tdcm :: predictions generated")
    dataset$tdcm_pred <- tdcm_pred1
    end_time <- Sys.time()
    time_taken = as.numeric(end_time - start_time, units = "mins")
    time = c(time, time_taken)
    data[[i]] = dataset
  }
  data[[length(data)+1]] <- time
  return(data)
}


# predictions from LM

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))]}

predictions_lm<-function(dev_data, prediction_datasets, model, landmark_times, PW){
  data = list()
  time = c()
  for (i in 1:length(prediction_datasets)){
    dataset = as.data.frame(prediction_datasets[[i]])
    dataset_LM = dataset
    dataset_LM$LM<-min(unique(closest(dev_data$LM, landmark_times[i])))
    start_time = Sys.time()
    predictions<-predict(model, newdata = dataset_LM, 
                         type = "surv", se.fit =TRUE)
    end_time = Sys.time()
    linear_predictor<-predict(model, newdata = dataset_LM, 
                              type = "lp", se.fit =TRUE)
    dataset$lm_pred<-predictions$fit
    dataset$lm_lp<-linear_predictor$fit
    time_taken = as.numeric(end_time - start_time, units = "mins")
    data[[i]] = dataset
    time = c(time, time_taken)
  }
  data[[length(data)+1]] <- time
  return(data)
}


## predictions from a joint model fitted using JMbayes2

predictions_jm<-function(model, long_data, surv_data, landmark_times, 
                          prediction_window){
  data = list()
  time = c()
  for (i in (1:length(landmark_times))){
    assessment_data <- long_data %>% filter(time <= landmark_times[i])
    assessment_data$Time <- landmark_times[i]
    horizon_time = landmark_times[i] + prediction_window
    start_time <- Sys.time()
    predictions <- JMbayes2:::predict.jm(model, newdata = assessment_data, 
                                         times = horizon_time, process="event", return_newdata=FALSE, cores = 1)
    pred_output <- data.frame(id = unlist(predictions$id),
                              times = unlist(predictions$times), 
                              preds = unlist(predictions$pred)) %>% 
      filter(times == landmark_times[i]+prediction_window)
    event_data <- surv_data
    event_data$pred <-1-pred_output$preds
    end_time <- Sys.time()
    event_data <- event_data %>% filter(Time > landmark_times[i])
    time_taken = as.numeric(end_time - start_time, units = "mins")
    data[[i]] <- event_data
    time <- c(time, time_taken)
  }
  data[[length(data)+1]] <- time
  return(data)
}

## prediction timing function

prediction_time_matrix<-function(landmark_times, prediction_times, model_id, predicion_window, i){
  time_matrix<-matrix(0, nrow = length(landmark_times), ncol = 6)
  
  for (j in 1:length(landmark_times)){
    time_matrix[j,]<-c(i, model_id, landmark_times[j], prediction_window, "pred_time", prediction_times[j])
  }
  return(time_matrix)  
}

## add stability prediction data from iteration

stability_data <-function(prediction_data, pred_var, landmark_times, model_id, i){
  data = list()
  for (j in 1:length(landmark_times)){
    model_preds <- data.frame(prediction_data[[j]])[, c("id", pred_var)]
    model_preds$simulation_id <-rep(i, times=length(model_preds$id))
    model_preds$model <- rep(model_id, times=length(model_preds$id))
    model_preds$landmark_time <- rep(landmark_times[j], times=length(model_preds$id))
    data[[j]] <- model_preds
  }
  output = do.call(rbind, data)
  names(output)<-c("SubjID", "Pred", "SimID", "Mod", "LT")
  return(output)
}
