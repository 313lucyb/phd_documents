# performance metric calculation

performance_assessment<-function(prediction_datasets, pred_var, lp_var, survival_time, event_status, 
                                 landmark_times, prediction_window, iter,
                                 model_num){
  results = list()
  for (i in 1:length(prediction_datasets)){
    dataset = data.frame(prediction_datasets[[i]])
    brier <- sbrier(obj=Surv(as.numeric(dplyr::select(dataset, all_of(survival_time))[,1]), dplyr::select(dataset, all_of(event_status))[,1]), pred=dplyr::select(dataset, all_of(pred_var))[,1],
                    btime=landmark_times[i] + prediction_window)
    cox <- coxph(Surv(dplyr::select(dataset, all_of(survival_time))[,1], dplyr::select(dataset, all_of(event_status))[,1])~dplyr::select(dataset, all_of(lp_var))[,1])
    c_index <-summary(cox)$concordance[1]
    c_index_lower <- summary(cox)$concordance[1]-
      (1.96*summary(cox)$concordance[2])
    c_index_upper <- summary(cox)$concordance[1]+
      (1.96*summary(cox)$concordance[2])
    values <- c(brier, c_index, c_index_lower, c_index_upper)
    #4 = number of metric values
    mini_matrix <- matrix(NA, nrow=4, ncol=6)
    mini_matrix[,1]<-rep(iter, times = 4)
    mini_matrix[,2]<-rep(model_num, times=4)
    mini_matrix[,3]<-rep(landmark_times[i], times = 4)
    mini_matrix[,4]<-rep(prediction_window, times=4)
    mini_matrix[,5]<-c("brier", "c_index", "c_index_lower", "c_index_upper")
    mini_matrix[,6]<-values
    results[[i]] <- mini_matrix
  }
  output<- do.call(rbind, results)
  return(output)
}


null_brier<-function(data, btime){
  mod<-coxph(Surv(Time, event) ~ 1, data)
  pred.mat <-predict(mod, type = "survival", data = data)
  data$survpred<-as.numeric(pred.mat)
  BS<-sbrier(obj = Surv(data$Time,data$event), 
             data$survpred, btime = btime)
  return(BS)
}

lm_model_output <- function(model, fit_time, development_data, 
                            val_dev_data, val_new_data, pred_var, 
                            lp_var, model_id, iter, landmark_times,
                            prediction_window){
  
  fit_times<-matrix(NA, nrow=1, ncol=6)
  fit_times[1,] <- c(iter, model_id, NA, NA, "fit_time", fit_time)
  
  ## generate predictions and prediction times from model
  
  dev_output <- predictions_lm(development_data, val_dev_data, 
                               model, landmark_times, 
                               prediction_window)
  dev_preds <- dev_output[1:length(landmark_times)]
  
  new_output <- predictions_lm(development_data, val_new_data, 
                               model, landmark_times, 
                               prediction_window)
  
  new_preds <- new_output[1:length(landmark_times)]
  new_preds_time <- new_output[[(length(landmark_times)+1)]]
  
  ## put prediction times into consistent structure across models
  
  pred_times<-prediction_time_matrix(landmark_times, new_preds_time, model_id, prediction_window, iter)
  
  ## extract predictions on validation data by model to assess stability
  
  stability <- stability_data(new_preds, pred_var, landmark_times, model_id, iter)
  stability$PW<-rep(prediction_window, length(stability$SubjID))
  stability <- stability %>%
    dplyr::select(SimID, SubjID, Mod, LT, PW, Pred)
  stability<-as.data.frame(stability)
  names(stability)<-c("simulation_id", "subject_id", "model", "landmark_time", 
                      "prediction_window", "prediction")
  stability<-as.data.frame(stability)
  
  ## predictive performance on validation and development data

  validation_data<-performance_assessment(new_preds, pred_var, lp_var, "Time", "event", 
                                          landmark_times, prediction_window, i, model_id)
  apparent_data<-performance_assessment(dev_preds, pred_var, lp_var, "Time", "event", 
                                        landmark_times, prediction_window, i, 7)
  optimism_values <-as.numeric(apparent_data[,6])-as.numeric(validation_data[,6])
  optimism<-matrix(NA, nrow = length(optimism_values), ncol = 6)
  optimism[,1]<-rep(iter, times = length(optimism_values))
  optimism[,2]<-rep(model_id, length(optimism_values))
  optimism[,3]<-c(rep(landmark_times, each = 4))
  optimism[,4]<-rep(prediction_window, times = length(optimism_values))
  optimism[,5]<-rep(c("brier_adjusted", "c_index_adjusted", 
                      "c_index_lower_adjusted", "c_index_higher_adjusted"), times=length(landmark_times))
  optimism[,6]<-optimism_values
  
  lm_output<- rbind(fit_times, pred_times, validation_data, optimism)
  lm_stability<-stability
  output_list <- list(lm_output, lm_stability)
  return(output_list)
}


model_output <- function(model, fit_time, 
                          val_dev_data, val_new_data, pred_var, 
                          lp_var, model_id, iter, landmark_times,
                          prediction_window, model_type){
  
  covariates <- c("age", "sex", "disease_duration", "clinician_disease",  
                    "patient_disease", "bmi")
  
  fit_times<-matrix(NA, nrow=1, ncol=6)
  fit_times[1,] <- c(iter, model_id, NA, NA, "fit_time", fit_time)
  
  ## generate predictions and prediction times from model
  if (model_type == "tdcm"){
    dev_output <- predictions_tdcm(val_dev_data, model, 
                                   landmark_times, 
                                   covariates, 
                                   prediction_window)
    new_output <- predictions_tdcm(val_new_data, model,
                                   landmark_times, 
                                   covariates, 
                                   prediction_window)
  } else {
    dev_output <- predictions_fpsm(val_dev_data, model, 
                                   landmark_times, 
                                   covariates, 
                                   prediction_window)
    new_output <- predictions_fpsm(val_new_data, model,
                                   landmark_times, 
                                   covariates, 
                                   prediction_window)
  }
  dev_preds <- dev_output[1:length(landmark_times)]
  print("model_output :: predictions generated on dev data")
  new_preds <- new_output[1:length(landmark_times)]
  new_preds_time <- new_output[[(length(landmark_times)+1)]]
  print("model_output :: predictions generated on new data data")
  ## put prediction times into consistent structure across models
  
  pred_times<-prediction_time_matrix(landmark_times, new_preds_time, model_id, prediction_window, iter)
  ## extract predictions on validation data by model to assess stability
  
  stability <- stability_data(new_preds, pred_var, landmark_times, model_id, iter)
  stability$PW<-rep(prediction_window, length(stability$SubjID))
  stability <- stability %>%
    dplyr::select(SimID, SubjID, Mod, LT, PW, Pred)
  stability<-as.data.frame(stability)
  names(stability)<-c("simulation_id", "subject_id", "model", "landmark_time", 
                      "prediction_window", "prediction")
  stability<-as.data.frame(stability)
  
  ## predictive performance on validation and development data
  
  validation_data<-performance_assessment(new_preds, pred_var, lp_var, "Time", "event", 
                                          landmark_times, prediction_window, i, model_id)
  apparent_data<-performance_assessment(dev_preds, pred_var, lp_var, "Time", "event", 
                                        landmark_times, prediction_window, i, 7)
  optimism_values <-as.numeric(apparent_data[,6])-as.numeric(validation_data[,6])
  optimism<-matrix(NA, nrow = length(optimism_values), ncol = 6)
  optimism[,1]<-rep(iter, times = length(optimism_values))
  optimism[,2]<-rep(model_id, length(optimism_values))
  optimism[,3]<-c(rep(landmark_times, each = 4))
  optimism[,4]<-rep(prediction_window, times = length(optimism_values))
  optimism[,5]<-rep(c("brier_adjusted", "c_index_adjusted", 
                      "c_index_lower_adjusted", "c_index_higher_adjusted"), times=length(landmark_times))
  optimism[,6]<-optimism_values
  
  output<- rbind(fit_times, pred_times, validation_data, optimism)
  stability<-stability
  output_list <- list(output, stability)
  return(output_list)
}
 
