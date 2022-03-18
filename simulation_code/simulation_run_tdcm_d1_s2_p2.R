
source("./data_generation_functions_new.R")
source("./prediction_functions_new.R")
source("./performance_functions.R")
source("./simulation_functions_new.R")

require(dplyr)
require(dynpred)
require(nlme)
require(ipred)
require(foreach)
require(doParallel)
require(tidyr)
require(data.table)
require(JMbayes2)
require(rstpm2)
require(rms)
require(splus2R)
require(truncnorm)

start <- Sys.time()

sims_parameters <- crossing(
  dgm = "tdcm",
  validation_sample_size = 10000,
  dimension = 1,
  n = 1000,
  j = 8,
  weibull_shape = 0.9,
  weibull_scale = c(0.002152, 0.01076),
  probability_missing = 0.8,
  prediction_window = 1,
  landmark_times = list(c(0, 1, 2))
)

args <- commandArgs(trailingOnly = T) # pull in args from qsub file
s <- as.numeric(args[1]) # specify class of argument

n_reps = 100

#set.seed(643546577 + s)
print("s")
print(s)
simulation_output <- across_sim_func_parallel(n_simulation = n_reps,
                                              dgm = sims_parameters$dgm[s],
                                              validation_sample_size = sims_parameters$validation_sample_size[s], 
                                              dimension = sims_parameters$dimension[s],
                                              n = sims_parameters$n[s], 
                                              j = sims_parameters$j[s], 
                                              weibull_shape = sims_parameters$weibull_shape[s],
                                              weibull_scale = sims_parameters$weibull_scale[s],
                                              probability_missing = sims_parameters$probability_missing[s],
                                              prediction_window = sims_parameters$prediction_window[s], 
                                              landmark_times = sims_parameters$landmark_times[s][[1]],
                                              scenario=s)

performance_results <-simulation_output[[1]]
stability_results <- simulation_output[[2]]


performance_results <- performance_results %>%
  mutate("simulation_scenario" = s,
         "dgm" = sims_parameters$dgm[s],
         "validation_sample_size" = sims_parameters$validation_sample_size[s],
         "dimension" = sims_parameters$dimension[s],
         "n" = sims_parameters$n[s], 
         "j" = sims_parameters$j[s], 
         "weibull_shape" = sims_parameters$weibull_shape[s],
         "weibull_scale" = sims_parameters$weibull_scale[s],
         "prediction_window" = sims_parameters$prediction_window[s], 
         "landmark_times" = toString(sims_parameters$landmark_times[s][[1]]),
         .after = "model_id")

stability_results <- stability_results %>%
  mutate("simulation_scenario" = s,
         "dgm" = sims_parameters$dgm[s],
         "validation_sample_size" = sims_parameters$validation_sample_size[s],
         "dimension" = sims_parameters$dimension[s],
         "n" = sims_parameters$n[s], 
         "j" = sims_parameters$j[s], 
         "weibull_shape" = sims_parameters$weibull_shape[s],
         "weibull_scale" = sims_parameters$weibull_scale[s],
         "prediction_window" = sims_parameters$prediction_window[s], 
         "landmark_times" = toString(sims_parameters$landmark_times[s][[1]]),
         .after = "model")

write.table(performance_results, file = paste("./performance_results_tdcm_d1_s2_p2_",s, ".csv", sep=""), sep = ",", col.names = NA,
            qmethod = "double")
write.table(stability_results, file = paste("./stability_results_tdcm_d1_s2_p2_",s, ".csv", sep=""), sep = ",", col.names = NA,
            qmethod = "double")

warnings()
end <- Sys.time()
time_taken <- as.numeric(end - start, units = "mins")
print("time_taken")
print(time_taken)
