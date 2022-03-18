setwd("R:/BSRBR/Analyses/lucy_bull/PhD work/Simulation study/Simulation_output_Feb22/")
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(glue)
cols <-brewer.pal(11, "RdBu")
cols_models <- c(cols[1], cols[3], cols[9], cols[11])


#######################
## Performance data ###
#######################

performance_data = read.csv(file ="performance_results_all.csv")

performance_data <- performance_data %>% 
  mutate(probability_missing = ifelse(stringr::str_detect(file.name, "p1"), 0.2, 0.8),
         weibull_scale = round(weibull_scale, 10),
         event_prevalence = ifelse(weibull_scale == round(exp(-13.331), 10), 10, 
                                   ifelse(weibull_scale == round(exp(-11.14453), 10), 40, 
                                          ifelse(weibull_scale == round(exp(-6.283), 10), 10, 
                                                 ifelse(weibull_scale == round(exp(-4.652), 10), 40, 
                                                        ifelse(weibull_scale == 0.002152, 10, 
                                                               ifelse(weibull_scale == 0.01076, 40, 
                                                                      ifelse(weibull_scale == 0.00257, 10, 40)))))))) %>%
  filter(model_id %in% c(1,2,3,8, NA))


performance_data$n <- factor(performance_data$n, levels=c(200, 1000),
                             labels = c("n = 200", "n = 1000"))
performance_data$event_prevalence <- factor(performance_data$event_prevalence, levels = c(10, 40),
                                            labels = c("prev ~ 10%", "prev ~ 40%"))
performance_data$probability_missing <- factor(performance_data$probability_missing, levels = c(0.2, 0.8),
                                               labels = c("20% missing", "80% missing"))
performance_data$dimension <- factor(performance_data$dimension, levels = c(1, 3),
                                     labels = c("1 longitudinal", "3 longitudinal"))
performance_data$dgm <- factor(performance_data$dgm, levels=c("jm", "tdcm"),
                               labels = c("DGM = JM", "DGM = TDCM"))

performance_data <- performance_data %>% select(simulation_id, 
                                                model_id, 
                                                dgm, 
                                                dimension, 
                                                probability_missing, 
                                                n, 
                                                weibull_shape, 
                                                weibull_scale, 
                                                landmark_time, 
                                                prediction_window, 
                                                statistic, 
                                                value, 
                                                event_prevalence)


null_brier_table <- performance_data %>% filter(statistic == "null_brier") %>% 
  group_by(dgm, n, event_prevalence, dimension, landmark_time) %>% 
  summarise(nullbrier = mean(value))

brier_table <- performance_data %>% filter(statistic == "brier")
brier_ipa <- left_join(brier_table, null_brier_table, by= c("dgm", "n", "event_prevalence",
                                                               "dimension", "landmark_time")) %>%
                mutate(ipa = (1 - value/nullbrier)*100, statistic = "ipa", value = as.numeric(ipa))

#######################
## Stability data   ###
#######################

stability_data = read.csv(file ="stability_results_all.csv")

stability_data <- stability_data %>% 
  mutate(probability_missing = ifelse(stringr::str_detect(file.name, "p1"), 0.2, 0.8),
         weibull_scale = round(weibull_scale, 10),
         event_prevalence = ifelse(weibull_scale == round(exp(-13.331), 10), 10, 
                                   ifelse(weibull_scale == round(exp(-11.14453), 10), 40, 
                                          ifelse(weibull_scale == round(exp(-6.283), 10), 10, 
                                                 ifelse(weibull_scale == round(exp(-4.652), 10), 40, 
                                                        ifelse(weibull_scale == 0.002152, 10, 
                                                               ifelse(weibull_scale == 0.01076, 40, 
                                                                      ifelse(weibull_scale == 0.00257, 10, 40))))))))
stability_data$n <- factor(stability_data$n, levels=c(200, 1000),
                           labels = c("n = 200", "n = 1000"))
stability_data$event_prevalence <- factor(stability_data$event_prevalence, levels = c(10, 40),
                                          labels = c("prev ~ 10%", "prev ~ 40%"))
stability_data$probability_missing <- factor(stability_data$probability_missing, levels = c(0.2, 0.8),
                                             labels = c("20% missing", "80% missing"))
stability_data$dimension <- factor(stability_data$dimension, levels = c(1, 3),
                                   labels = c("1 longitudinal", "3 longitudinal"))
stability_data$dgm <- factor(stability_data$dgm, levels=c("jm", "tdcm"),
                               labels = c("DGM = JM", "DGM = TDCM"))
stability_data <- stability_data %>% filter(model %in% c(1,2,3,8))

source("./plot_functions.R")

#################################
# Single longitudinal predictor #
#################################


##########################
### Comparing extremes ###
##########################
#************************#

#######################
####### C-index #######
#######################


high_extreme_c <- performance_plot_dgm(data = performance_data,
                                       stat = "c_index", 
                                       ylimits = c(0.5,0.9), 
                                       ylabel= " ", 
                                       dim = "1 longitudinal", 
                                       prob = "20% missing", 
                                       sample = "n = 1000", 
                                       prev = "prev ~ 40%",
                                       legend = "right")
low_extreme_c <- performance_plot_dgm(data = performance_data, 
                                      stat = "c_index", 
                                      ylimits = c(0.5,0.9), 
                                      ylabel= "Harrell's C-index (95% quantile range)", 
                                      dim = "1 longitudinal", 
                                      prob = "80% missing", 
                                      sample = "n = 200", 
                                      prev = "prev ~ 10%",
                                      legend = "none")

low_extreme_c | high_extreme_c

#######################
####### Brier   #######
#######################
high_extreme_brier <- performance_plot_dgm(data = performance_data,
                                           stat = "brier", 
                                           ylimits = c(0, 0.175), 
                                           ylabel= " ", 
                                           dim = "1 longitudinal", 
                                           prob = "20% missing", 
                                           sample = "n = 1000", 
                                           prev = "prev ~ 40%",
                                           legend = "right")
low_extreme_brier <- performance_plot_dgm(data = performance_data,
                                          stat = "brier", 
                                          ylimits = c(0, 0.175), 
                                          ylabel= "Survival Brier score (95% quantile range)", 
                                          dim = "1 longitudinal", 
                                          prob = "80% missing", 
                                          sample = "n = 200", 
                                          prev = "prev ~ 10%",
                                          legend = "none")

low_extreme_brier | high_extreme_brier

#######################
####### IPA     #######
#######################

high_extreme_ipa <- performance_plot_dgm(data = brier_ipa,
                                           stat = "ipa", 
                                           ylimits = c(-25, 75), 
                                           ylabel= " ", 
                                           dim = "1 longitudinal", 
                                           prob = "20% missing", 
                                           sample = "n = 1000", 
                                           prev = "prev ~ 40%",
                                           legend = "right")
low_extreme_ipa <- performance_plot_dgm(data = brier_ipa,
                                          stat = "ipa", 
                                          ylimits = c(-25, 75), 
                                          ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                          dim = "1 longitudinal", 
                                          prob = "80% missing", 
                                          sample = "n = 200", 
                                          prev = "prev ~ 10%",
                                          legend = "none")

low_extreme_ipa | high_extreme_ipa

##########################
### Optimism - C-index ###
##########################

high_extreme_optC <- performance_plot_dgm(data = performance_data,
                                         stat = "c_index_adjusted", 
                                         ylimits = c(-0.25, 0.5), 
                                         ylabel= " ", 
                                         dim = "1 longitudinal", 
                                         prob = "20% missing", 
                                         sample = "n = 1000", 
                                         prev = "prev ~ 40%",
                                         legend = "right")

low_extreme_optC <- performance_plot_dgm(data = performance_data,
                                        stat = "c_index_adjusted", 
                                        ylimits = c(-0.25, 0.5), 
                                        ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                        dim = "1 longitudinal", 
                                        prob = "80% missing", 
                                        sample = "n = 200", 
                                        prev = "prev ~ 10%",
                                        legend = "none")

low_extreme_optC | high_extreme_optC

##########################
### Optimism - Brier   ###
##########################

high_extreme_optB <- performance_plot_dgm(data = performance_data,
                                          stat = "brier_adjusted", 
                                          ylimits = c(-0.03, 0.03), 
                                          ylabel= " ", 
                                          dim = "1 longitudinal", 
                                          prob = "20% missing", 
                                          sample = "n = 1000", 
                                          prev = "prev ~ 40%",
                                          legend = "right")

low_extreme_optB <- performance_plot_dgm(data = performance_data,
                                         stat = "brier_adjusted", 
                                         ylimits = c(-0.03, 0.03), 
                                         ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                         dim = "1 longitudinal", 
                                         prob = "80% missing", 
                                         sample = "n = 200", 
                                         prev = "prev ~ 10%",
                                         legend = "none")

low_extreme_optB | high_extreme_optB

#####################
#### Stability   ####
#####################

high_extreme_stability <- stability_plot_dgm(stability_data, 
                                            ylimits = c(0, 0.9), 
                                            ylabel = " ", 
                                            dim = "1 longitudinal", 
                                            prob = "20% missing", 
                                            sample = "n = 1000", 
                                            prev = "prev ~ 40%", 
                                            legend = "right")

low_extreme_stability <- stability_plot_dgm(stability_data, 
                                            ylimits = c(0, 0.9), 
                                            ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                            dim = "1 longitudinal", 
                                            prob = "20% missing", 
                                            sample = "n = 200", 
                                            prev = "prev ~ 10%", 
                                            legend = "none")

low_extreme_stability | high_extreme_stability


###################
### Consistency ###
###################

high_extreme_consistencyC <- consistency_plot_dgm(data = performance_data,
                                               stat = "c_index",
                                               dim = "1 longitudinal",
                                               prob = "20% missing",
                                               sample = "n = 1000",
                                               prev = "prev ~ 40%",
                                               ylabel = " ",
                                               ylimits = c(0, 0.07),
                                               legend = "none")
low_extreme_consistencyC <- consistency_plot_dgm(data = performance_data,
                                               stat = "c_index",
                                               dim = "1 longitudinal",
                                               prob = "80% missing",
                                               sample = "n = 200",
                                               prev = "prev ~ 10%",
                                               ylabel = "Harrell's C-index range\n across landmark times",
                                               ylimits = c(0, 0.07),
                                               legend = "none")

low_extreme_consistencyC | high_extreme_consistencyC



high_extreme_consistencyB <- consistency_plot_dgm(data = performance_data,
                                                  stat = "brier",
                                                  dim = "1 longitudinal",
                                                  prob = "20% missing",
                                                  sample = "n = 1000",
                                                  prev = "prev ~ 40%",
                                                  ylabel = " ",
                                                  ylimits = c(0, 0.1),
                                                  legend = "none")
low_extreme_consistencyB <- consistency_plot_dgm(data = performance_data,
                                                 stat = "brier",
                                                 dim = "1 longitudinal",
                                                 prob = "80% missing",
                                                 sample = "n = 200",
                                                 prev = "prev ~ 10%",
                                                 ylabel = "Harrell's C-index range\n across landmark times",
                                                 ylimits = c(0, 0.1),
                                                 legend = "none")

low_extreme_consistencyB | high_extreme_consistencyB

#####################
### Fitting times ###
#####################

high_extreme_fitting <- fitting_time_dgm(data = performance_data,
                                       stat = "fit_time", 
                                       ylimits = c(0,20), 
                                       ylabel= " ", 
                                       dim = "3 longitudinal", 
                                       prob = "20% missing", 
                                       sample = "n = 1000", 
                                       prev = "prev ~ 40%",
                                       legend = "none")
low_extreme_fitting <- fitting_time_dgm(data = performance_data,
                                        stat = "fit_time", 
                                        ylimits = c(0,20), 
                                        ylabel= "Model fitting time (minutes)", 
                                        dim = "1 longitudinal", 
                                        prob = "80% missing", 
                                        sample = "n = 200", 
                                        prev = "prev ~ 10%",
                                        legend = "none")

low_extreme_fitting | high_extreme_fitting

########################
### Prediction times ###
########################

high_extreme_pred_time <- performance_plot_dgm(data = performance_data,
                                       stat = "pred_time", 
                                       ylimits = c(0,50), 
                                       ylabel= " ", 
                                       dim = "1 longitudinal", 
                                       prob = "20% missing", 
                                       sample = "n = 1000", 
                                       prev = "prev ~ 40%",
                                       legend = "right")

low_extreme_pred_time <- performance_plot_dgm(data = performance_data,
                                               stat = "pred_time", 
                                               ylimits = c(0,50), 
                                               ylabel= "Prediction time (minutes)", 
                                               dim = "1 longitudinal", 
                                               prob = "80% missing", 
                                               sample = "n = 200", 
                                               prev = "prev ~ 10%",
                                               legend = "none")
low_extreme_pred_time | high_extreme_pred_time

#########################
### Fitting failures  ###
#########################


single_failed <- failed_dgm_single(data = performance_data, dim = "1 longitudinal", col = 4)
single_failed
multiple_failed <- failed_dgm_multiple(data = performance_data, dim = "3 longitudinal", col = 3)
multiple_failed


#############################
### Sample size influence ###
#############################
#***************************#

#######################
####### C-index #######
#######################


high_c_sample_jm <- performance_plot_sample(data = performance_data,
                                       stat = "c_index", 
                                       ylimits = c(0.5,0.9), 
                                       ylabel= "Harrell's C-index (95% quantile range)", 
                                       dim = "1 longitudinal", 
                                       prob = "20% missing", 
                                       dgm_arg = "DGM = JM", 
                                       prev = "prev ~ 40%",
                                       legend = "right")
high_c_sample_jm

low_c_sample_jm <- performance_plot_sample(data = performance_data, 
                                      stat = "c_index", 
                                      ylimits = c(0.5,0.9), 
                                      ylabel= "Harrell's C-index (95% quantile range)", 
                                      dim = "1 longitudinal", 
                                      prob = "80% missing", 
                                      dgm_arg = "DGM = JM",
                                      prev = "prev ~ 10%",
                                      legend = "right")
low_c_sample_jm

high_c_sample_tdcm <- performance_plot_sample(data = performance_data,
                                            stat = "c_index", 
                                            ylimits = c(0.5,0.9), 
                                            ylabel= "Harrell's C-index (95% quantile range)", 
                                            dim = "1 longitudinal", 
                                            prob = "20% missing", 
                                            dgm_arg = "DGM = TDCM", 
                                            prev = "prev ~ 40%",
                                            legend = "right")
high_c_sample_tdcm

low_c_sample_tdcm <- performance_plot_sample(data = performance_data, 
                                           stat = "c_index", 
                                           ylimits = c(0.5,0.9), 
                                           ylabel= "Harrell's C-index (95% quantile range)", 
                                           dim = "1 longitudinal", 
                                           prob = "80% missing", 
                                           dgm_arg = "DGM = TDCM",
                                           prev = "prev ~ 10%",
                                           legend = "right")
low_c_sample_tdcm

#######################
####### Brier   #######
#######################

high_brier_sample_jm <- performance_plot_sample(data = performance_data,
                                           stat = "brier", 
                                           ylimits = c(0, 0.175), 
                                           ylabel= "Survival Brier score (95% quantile range)", 
                                           dim = "1 longitudinal", 
                                           prob = "20% missing", 
                                           dgm_arg = "DGM = JM",
                                           prev = "prev ~ 40%",
                                           legend = "right")
high_brier_sample_jm

low_brier_sample_jm <- performance_plot_sample(data = performance_data,
                                          stat = "brier", 
                                          ylimits = c(0, 0.05), 
                                          ylabel= "Survival Brier score (95% quantile range)", 
                                          dim = "1 longitudinal", 
                                          prob = "80% missing", 
                                          dgm_arg = "DGM = JM",
                                          prev = "prev ~ 10%",
                                          legend = "right")
low_brier_sample_jm

high_brier_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                stat = "brier", 
                                                ylimits = c(0, 0.175), 
                                                ylabel= "Survival Brier score (95% quantile range)", 
                                                dim = "1 longitudinal", 
                                                prob = "20% missing", 
                                                dgm_arg = "DGM = TDCM",
                                                prev = "prev ~ 40%",
                                                legend = "right")
high_brier_sample_tdcm

low_brier_sample_tdcm <- performance_plot_sample(data = performance_data,
                                               stat = "brier", 
                                               ylimits = c(0, 0.05), 
                                               ylabel= "Survival Brier score (95% quantile range)", 
                                               dim = "1 longitudinal", 
                                               prob = "80% missing", 
                                               dgm_arg = "DGM = TDCM",
                                               prev = "prev ~ 10%",
                                               legend = "right")
low_brier_sample_tdcm


#######################
####### IPA     #######
#######################

high_ipa_sample_jm <- performance_plot_sample(data = brier_ipa,
                                         stat = "ipa", 
                                         ylimits = c(-25, 75), 
                                         ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                         dim = "1 longitudinal", 
                                         prob = "20% missing", 
                                         dgm_arg = "DGM = JM", 
                                         prev = "prev ~ 40%",
                                         legend = "right")
high_ipa_sample_jm
low_ipa_sample_jm <- performance_plot_sample(data = brier_ipa,
                                        stat = "ipa", 
                                        ylimits = c(-25, 75), 
                                        ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                        dim = "1 longitudinal", 
                                        prob = "80% missing", 
                                        dgm_arg = "DGM = JM",  
                                        prev = "prev ~ 10%",
                                        legend = "right")
low_ipa_sample_jm
high_ipa_sample_tdcm <- performance_plot_sample(data = brier_ipa,
                                              stat = "ipa", 
                                              ylimits = c(-25, 75), 
                                              ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                              dim = "1 longitudinal", 
                                              prob = "20% missing", 
                                              dgm_arg = "DGM = TDCM", 
                                              prev = "prev ~ 40%",
                                              legend = "right")
high_ipa_sample_tdcm
low_ipa_sample_tdcm <- performance_plot_sample(data = brier_ipa,
                                             stat = "ipa", 
                                             ylimits = c(-25, 75), 
                                             ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                             dim = "1 longitudinal", 
                                             prob = "80% missing", 
                                             dgm_arg = "DGM = TDCM",  
                                             prev = "prev ~ 10%",
                                             legend = "right")
low_ipa_sample_tdcm


##########################
### Optimism - C-index ###
##########################

high_optC_sample_jm <- performance_plot_sample(data = performance_data,
                                          stat = "c_index_adjusted", 
                                          ylimits = c(-0.2, 0.2), 
                                          ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                          dim = "1 longitudinal", 
                                          prob = "20% missing", 
                                          dgm_arg = "DGM = JM", 
                                          prev = "prev ~ 40%",
                                          legend = "right")
high_optC_sample_jm
low_optC_sample_jm <- performance_plot_sample(data = performance_data,
                                         stat = "c_index_adjusted", 
                                         ylimits = c(-0.25, 0.5), 
                                         ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                         dim = "1 longitudinal", 
                                         prob = "80% missing", 
                                         dgm_arg = "DGM = JM",  
                                         prev = "prev ~ 10%",
                                         legend = "right")
low_optC_sample_jm
high_optC_sample_tdcm <- performance_plot_sample(data = performance_data,
                                               stat = "c_index_adjusted", 
                                               ylimits = c(-0.2, 0.2), 
                                               ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                               dim = "1 longitudinal", 
                                               prob = "20% missing", 
                                               dgm_arg = "DGM = TDCM", 
                                               prev = "prev ~ 40%",
                                               legend = "right")
high_optC_sample_tdcm
low_optC_sample_tdcm <- performance_plot_sample(data = performance_data,
                                              stat = "c_index_adjusted", 
                                              ylimits = c(-0.25, 0.5), 
                                              ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                              dim = "1 longitudinal", 
                                              prob = "80% missing", 
                                              dgm_arg = "DGM = TDCM",  
                                              prev = "prev ~ 10%",
                                              legend = "right")
low_optC_sample_tdcm 

##########################
### Optimism - Brier   ###
##########################

high_optB_sample_jm <- performance_plot_sample(data = performance_data,
                                          stat = "brier_adjusted", 
                                          ylimits = c(-0.07, 0.07), 
                                          ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                          dim = "1 longitudinal", 
                                          prob = "20% missing", 
                                          dgm_arg = "DGM = JM",  
                                          prev = "prev ~ 40%",
                                          legend = "right")
high_optB_sample_jm
low_optB_sample_jm <- performance_plot_sample(data = performance_data,
                                         stat = "brier_adjusted", 
                                         ylimits = c(-0.07, 0.07), 
                                         ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                         dim = "1 longitudinal", 
                                         prob = "80% missing", 
                                         dgm_arg = "DGM = JM", 
                                         prev = "prev ~ 10%",
                                         legend = "right")

low_optB_sample_jm
high_optB_sample_tdcm <- performance_plot_sample(data = performance_data,
                                               stat = "brier_adjusted", 
                                               ylimits = c(-0.07, 0.07), 
                                               ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                               dim = "1 longitudinal", 
                                               prob = "20% missing", 
                                               dgm_arg = "DGM = TDCM",  
                                               prev = "prev ~ 40%",
                                               legend = "right")
high_optB_sample_tdcm
low_optB_sample_tdcm <- performance_plot_sample(data = performance_data,
                                              stat = "brier_adjusted", 
                                              ylimits = c(-0.07, 0.07), 
                                              ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                              dim = "1 longitudinal", 
                                              prob = "80% missing", 
                                              dgm_arg = "DGM = TDCM", 
                                              prev = "prev ~ 10%",
                                              legend = "right")

low_optB_sample_tdcm

#####################
#### Stability   ####
#####################

high_stability_sample_jm <- stability_plot_sample(stability_data, 
                                             ylimits = c(0, 0.9), 
                                             ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                             dim = "1 longitudinal", 
                                             prob = "20% missing", 
                                             dgm_arg = "DGM = JM", 
                                             prev = "prev ~ 40%", 
                                             legend = "right")
high_stability_sample_jm
low_stability_sample_jm <- stability_plot_sample(stability_data, 
                                            ylimits = c(0, 0.9), 
                                            ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                            dim = "1 longitudinal", 
                                            prob = "80% missing", 
                                            dgm_arg = "DGM = JM", 
                                            prev = "prev ~ 10%", 
                                            legend = "right")
low_stability_sample_jm
high_stability_sample_tdcm <- stability_plot_sample(stability_data, 
                                                  ylimits = c(0, 0.9), 
                                                  ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                  dim = "1 longitudinal", 
                                                  prob = "20% missing", 
                                                  dgm_arg = "DGM = TDCM", 
                                                  prev = "prev ~ 40%", 
                                                  legend = "right")
high_stability_sample_tdcm 
low_stability_sample_tdcm <- stability_plot_sample(stability_data, 
                                                 ylimits = c(0, 0.9), 
                                                 ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                 dim = "1 longitudinal", 
                                                 prob = "80% missing", 
                                                 dgm_arg = "DGM = TDCM", 
                                                 prev = "prev ~ 10%", 
                                                 legend = "right")
low_stability_sample_tdcm



###################
### Consistency ###
###################

high_consistencyC_sample_jm <- consistency_plot_sample(data = performance_data,
                                                  stat = "c_index",
                                                  dim = "1 longitudinal",
                                                  prob = "20% missing",
                                                  dgm_arg = "DGM = JM",
                                                  prev = "prev ~ 40%",
                                                  ylabel = "Harrell's C-index range\n across landmark times",
                                                  ylimits = c(0, 0.1),
                                                  legend = "right")
high_consistencyC_sample_jm
low_consistencyC_sample_jm <- consistency_plot_sample(data = performance_data,
                                                 stat = "c_index",
                                                 dim = "1 longitudinal",
                                                 prob = "80% missing",
                                                 dgm_arg = "DGM = JM",
                                                 prev = "prev ~ 10%",
                                                 ylabel = "Harrell's C-index range\n across landmark times",
                                                 ylimits = c(0, 0.07),
                                                 legend = "right")
low_consistencyC_sample_jm
high_consistencyC_sample_tdcm <- consistency_plot_sample(data = performance_data,
                                                       stat = "c_index",
                                                       dim = "1 longitudinal",
                                                       prob = "20% missing",
                                                       dgm_arg = "DGM = TDCM",
                                                       prev = "prev ~ 40%",
                                                       ylabel = "Harrell's C-index range\n across landmark times",
                                                       ylimits = c(0, 0.07),
                                                       legend = "right")
high_consistencyC_sample_tdcm
low_consistencyC_sample_tdcm <- consistency_plot_sample(data = performance_data,
                                                      stat = "c_index",
                                                      dim = "1 longitudinal",
                                                      prob = "80% missing",
                                                      dgm_arg = "DGM = TDCM",
                                                      prev = "prev ~ 10%",
                                                      ylabel = "Harrell's C-index range\n across landmark times",
                                                      ylimits = c(0, 0.07),
                                                      legend = "right")
low_consistencyC_sample_tdcm



high_consistencyB_sample_jm <- consistency_plot_sample(data = performance_data,
                                                  stat = "brier",
                                                  dim = "1 longitudinal",
                                                  prob = "20% missing",
                                                  dgm_arg = "DGM = JM",
                                                  prev = "prev ~ 40%",
                                                  ylabel = "Brier score range\n across landmark times",
                                                  ylimits = c(0, 0.1),
                                                  legend = "right")
high_consistencyB_sample_jm
low_consistencyB_sample_jm <- consistency_plot_sample(data = performance_data,
                                                 stat = "brier",
                                                 dim = "1 longitudinal",
                                                 prob = "80% missing",
                                                 dgm_arg = "DGM = JM",
                                                 prev = "prev ~ 10%",
                                                 ylabel = "Brier score range\n across landmark times",
                                                 ylimits = c(0, 0.1),
                                                 legend = "right")
low_consistencyB_sample_jm
high_consistencyB_sample_tdcm <- consistency_plot_sample(data = performance_data,
                                                       stat = "brier",
                                                       dim = "1 longitudinal",
                                                       prob = "20% missing",
                                                       dgm_arg = "DGM = TDCM",
                                                       prev = "prev ~ 40%",
                                                       ylabel = "Brier score range\n across landmark times",
                                                       ylimits = c(0, 0.1),
                                                       legend = "right")
high_consistencyB_sample_tdcm
low_consistencyB_sample_tdcm <- consistency_plot_sample(data = performance_data,
                                                      stat = "brier",
                                                      dim = "1 longitudinal",
                                                      prob = "80% missing",
                                                      dgm_arg = "DGM = TDCM",
                                                      prev = "prev ~ 10%",
                                                      ylabel = "Brier score range\n across landmark times",
                                                      ylimits = c(0, 0.1),
                                                      legend = "right")
low_consistencyB_sample_tdcm

#####################
### Fitting times ###
#####################

high_fitting_sample_jm <- fitting_time_sample(data = performance_data,
                                         stat = "fit_time", 
                                         ylimits = c(0,3), 
                                         ylabel= "Model fitting time (minutes)", 
                                         dim = "1 longitudinal", 
                                         prob = "20% missing", 
                                         dgm_arg = "DGM = JM", 
                                         prev = "prev ~ 40%",
                                         legend = "none")
high_fitting_sample_jm
low_fitting_sample_jm <- fitting_time_sample(data = performance_data,
                                        stat = "fit_time", 
                                        ylimits = c(0,3), 
                                        ylabel= "Model fitting time (minutes)", 
                                        dim = "1 longitudinal", 
                                        prob = "80% missing",
                                        dgm_arg = "DGM = JM",
                                        prev = "prev ~ 10%",
                                        legend = "right")
low_fitting_sample_jm
high_fitting_sample_tdcm <- fitting_time_sample(data = performance_data,
                                              stat = "fit_time", 
                                              ylimits = c(0,3), 
                                              ylabel= "Model fitting time (minutes)", 
                                              dim = "1 longitudinal", 
                                              prob = "20% missing", 
                                              dgm_arg = "DGM = TDCM", 
                                              prev = "prev ~ 40%",
                                              legend = "none")
high_fitting_sample_tdcm
low_fitting_sample_tdcm <- fitting_time_sample(data = performance_data,
                                             stat = "fit_time", 
                                             ylimits = c(0,20), 
                                             ylabel= "Model fitting time (minutes)", 
                                             dim = "1 longitudinal", 
                                             prob = "80% missing",
                                             dgm_arg = "DGM = TDCM",
                                             prev = "prev ~ 10%",
                                             legend = "right")
low_fitting_sample_tdcm


########################
### Prediction times ###
########################


high_pred_time_sample_jm <- performance_plot_sample(data = performance_data,
                                               stat = "pred_time", 
                                               ylimits = c(0,10), 
                                               ylabel= "Prediction time (minutes)", 
                                               dim = "1 longitudinal", 
                                               prob = "20% missing", 
                                               dgm_arg = "DGM = JM", 
                                               prev = "prev ~ 40%",
                                               legend = "right")
high_pred_time_sample_jm
low_pred_time_sample_jm <- performance_plot_sample(data = performance_data,
                                              stat = "pred_time", 
                                              ylimits = c(0,10), 
                                              ylabel= "Prediction time (minutes)", 
                                              dim = "1 longitudinal", 
                                              prob = "80% missing", 
                                              dgm_arg = "DGM = JM", 
                                              prev = "prev ~ 10%",
                                              legend = "right")
low_pred_time_sample_jm 
high_pred_time_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                    stat = "pred_time", 
                                                    ylimits = c(0,20), 
                                                    ylabel= "Prediction time (minutes)", 
                                                    dim = "1 longitudinal", 
                                                    prob = "20% missing", 
                                                    dgm_arg = "DGM = TDCM", 
                                                    prev = "prev ~ 40%",
                                                    legend = "right")
high_pred_time_sample_tdcm
low_pred_time_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                   stat = "pred_time", 
                                                   ylimits = c(0,30), 
                                                   ylabel= "Prediction time (minutes)", 
                                                   dim = "1 longitudinal", 
                                                   prob = "80% missing", 
                                                   dgm_arg = "DGM = TDCM", 
                                                   prev = "prev ~ 10%",
                                                   legend = "right")
low_pred_time_sample_tdcm 

##################################
### Event prevalence influence ###
##################################
#********************************#

#######################
####### C-index #######
#######################


high_c_prev_jm <- performance_plot_prev(data = performance_data,
                                            stat = "c_index", 
                                            ylimits = c(0.5,0.9), 
                                            ylabel= "Harrell's C-index (95% quantile range)", 
                                            dim = "1 longitudinal", 
                                            prob = "20% missing", 
                                            dgm_arg = "DGM = JM", 
                                            sample = "n = 1000",
                                            legend = "right")
high_c_prev_jm

low_c_prev_jm <- performance_plot_prev(data = performance_data, 
                                           stat = "c_index", 
                                           ylimits = c(0.5,0.9), 
                                           ylabel= "Harrell's C-index (95% quantile range)", 
                                           dim = "1 longitudinal", 
                                           prob = "80% missing", 
                                           dgm_arg = "DGM = JM",
                                          sample = "n = 200",
                                           legend = "right")
low_c_prev_jm

high_c_prev_tdcm <- performance_plot_prev(data = performance_data,
                                              stat = "c_index", 
                                              ylimits = c(0.5,0.9), 
                                              ylabel= "Harrell's C-index (95% quantile range)", 
                                              dim = "1 longitudinal", 
                                              prob = "20% missing", 
                                              dgm_arg = "DGM = TDCM", 
                                              sample = "n = 1000",
                                              legend = "right")
high_c_prev_tdcm

low_c_prev_tdcm <- performance_plot_prev(data = performance_data, 
                                             stat = "c_index", 
                                             ylimits = c(0.5,0.9), 
                                             ylabel= "Harrell's C-index (95% quantile range)", 
                                             dim = "1 longitudinal", 
                                             prob = "80% missing", 
                                             dgm_arg = "DGM = TDCM",
                                             sample = "n = 200",
                                             legend = "right")
low_c_prev_tdcm

#######################
####### Brier   #######
#######################

high_brier_prev_jm <- performance_plot_prev(data = performance_data,
                                                stat = "brier", 
                                                ylimits = c(0, 0.175), 
                                                ylabel= "Survival Brier score (95% quantile range)", 
                                                dim = "1 longitudinal", 
                                                prob = "20% missing", 
                                                dgm_arg = "DGM = JM",
                                                sample = "n = 1000",
                                                legend = "right")
high_brier_prev_jm

low_brier_prev_jm <- performance_plot_prev(data = performance_data,
                                               stat = "brier", 
                                               ylimits = c(0, 0.175), 
                                               ylabel= "Survival Brier score (95% quantile range)", 
                                               dim = "1 longitudinal", 
                                               prob = "80% missing", 
                                               dgm_arg = "DGM = JM",
                                               sample = "n = 200",
                                               legend = "right")
low_brier_prev_jm

high_brier_prev_tdcm <- performance_plot_prev(data = performance_data,
                                                  stat = "brier", 
                                                  ylimits = c(0, 0.175), 
                                                  ylabel= "Survival Brier score (95% quantile range)", 
                                                  dim = "1 longitudinal", 
                                                  prob = "20% missing", 
                                                  dgm_arg = "DGM = TDCM",
                                                  sample = "n = 1000",
                                                  legend = "right")
high_brier_prev_tdcm

low_brier_prev_tdcm <- performance_plot_prev(data = performance_data,
                                                 stat = "brier", 
                                                 ylimits = c(0, 0.175), 
                                                 ylabel= "Survival Brier score (95% quantile range)", 
                                                 dim = "1 longitudinal", 
                                                 prob = "80% missing", 
                                                 dgm_arg = "DGM = TDCM",
                                                 sample = "n = 200",
                                                 legend = "right")
low_brier_prev_tdcm


#######################
####### IPA     #######
#######################

high_ipa_prev_jm <- performance_plot_prev(data = brier_ipa,
                                              stat = "ipa", 
                                              ylimits = c(-25, 75), 
                                              ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                              dim = "1 longitudinal", 
                                              prob = "20% missing", 
                                              dgm_arg = "DGM = JM", 
                                              sample = "n = 1000",
                                              legend = "right")
high_ipa_prev_jm
low_ipa_prev_jm <- performance_plot_prev(data = brier_ipa,
                                             stat = "ipa", 
                                             ylimits = c(-25, 100), 
                                             ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                             dim = "1 longitudinal", 
                                             prob = "80% missing", 
                                             dgm_arg = "DGM = JM",  
                                             sample = "n = 200",
                                             legend = "right")
low_ipa_prev_jm
high_ipa_prev_tdcm <- performance_plot_prev(data = brier_ipa,
                                                stat = "ipa", 
                                                ylimits = c(-25, 75), 
                                                ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                                dim = "1 longitudinal", 
                                                prob = "20% missing", 
                                                dgm_arg = "DGM = TDCM", 
                                                sample = "n = 1000",
                                                legend = "right")
high_ipa_prev_tdcm
low_ipa_prev_tdcm <- performance_plot_prev(data = brier_ipa,
                                               stat = "ipa", 
                                               ylimits = c(-25, 75), 
                                               ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                               dim = "1 longitudinal", 
                                               prob = "80% missing", 
                                               dgm_arg = "DGM = TDCM",  
                                               sample = "n = 200",
                                               legend = "right")
low_ipa_prev_tdcm


##########################
### Optimism - C-index ###
##########################

high_optC_prev_jm <- performance_plot_prev(data = performance_data,
                                               stat = "c_index_adjusted", 
                                               ylimits = c(-0.2, 0.2), 
                                               ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                               dim = "1 longitudinal", 
                                               prob = "20% missing", 
                                               dgm_arg = "DGM = JM", 
                                               sample = "n = 1000",
                                               legend = "right")
high_optC_prev_jm
low_optC_prev_jm <- performance_plot_prev(data = performance_data,
                                              stat = "c_index_adjusted", 
                                              ylimits = c(-0.25, 0.5), 
                                              ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                              dim = "1 longitudinal", 
                                              prob = "80% missing", 
                                              dgm_arg = "DGM = JM",  
                                              sample = "n = 200",
                                              legend = "right")
low_optC_prev_jm
high_optC_prev_tdcm <- performance_plot_prev(data = performance_data,
                                                 stat = "c_index_adjusted", 
                                                 ylimits = c(-0.2, 0.2), 
                                                 ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                                 dim = "1 longitudinal", 
                                                 prob = "20% missing", 
                                                 dgm_arg = "DGM = TDCM", 
                                                 sample = "n = 1000",
                                                 legend = "right")
high_optC_prev_tdcm
low_optC_prev_tdcm <- performance_plot_prev(data = performance_data,
                                                stat = "c_index_adjusted", 
                                                ylimits = c(-0.25, 0.5), 
                                                ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                                dim = "1 longitudinal", 
                                                prob = "80% missing", 
                                                dgm_arg = "DGM = TDCM",  
                                                sample = "n = 200",
                                                legend = "right")
low_optC_prev_tdcm 

##########################
### Optimism - Brier   ###
##########################

high_optB_prev_jm <- performance_plot_prev(data = performance_data,
                                               stat = "brier_adjusted", 
                                               ylimits = c(-0.07, 0.07), 
                                               ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                               dim = "1 longitudinal", 
                                               prob = "20% missing", 
                                               dgm_arg = "DGM = JM",  
                                               sample = "n = 1000",
                                               legend = "right")
high_optB_prev_jm
low_optB_prev_jm <- performance_plot_prev(data = performance_data,
                                              stat = "brier_adjusted", 
                                              ylimits = c(-0.07, 0.07), 
                                              ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                              dim = "1 longitudinal", 
                                              prob = "80% missing", 
                                              dgm_arg = "DGM = JM", 
                                              sample = "n = 200",
                                              legend = "right")

low_optB_prev_jm
high_optB_prev_tdcm <- performance_plot_prev(data = performance_data,
                                           stat = "brier_adjusted", 
                                           ylimits = c(-0.07, 0.07), 
                                           ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                           dim = "1 longitudinal", 
                                           prob = "20% missing", 
                                           dgm_arg = "DGM = TDCM",  
                                           sample = "n = 1000",
                                           legend = "right")
high_optB_prev_tdcm
low_optB_prev_tdcm <- performance_plot_prev(data = performance_data,
                                          stat = "brier_adjusted", 
                                          ylimits = c(-0.07, 0.07), 
                                          ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                          dim = "1 longitudinal", 
                                          prob = "80% missing", 
                                          dgm_arg = "DGM = TDCM", 
                                          sample = "n = 200",
                                          legend = "right")

low_optB_prev_tdcm


#####################
#### Stability   ####
#####################

high_stability_prev_jm <- stability_plot_prev(stability_data, 
                                                  ylimits = c(0, 0.9), 
                                                  ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                  dim = "1 longitudinal", 
                                                  prob = "20% missing", 
                                                  dgm_arg = "DGM = JM", 
                                                  sample = "n = 1000", 
                                                  legend = "right")
high_stability_prev_jm
low_stability_prev_jm <- stability_plot_prev(stability_data, 
                                                 ylimits = c(0, 0.9), 
                                                 ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                 dim = "1 longitudinal", 
                                                 prob = "80% missing", 
                                                 dgm_arg = "DGM = JM", 
                                                  sample = "n = 200", 
                                                 legend = "right")
low_stability_prev_jm
high_stability_prev_tdcm <- stability_plot_prev(stability_data, 
                                                    ylimits = c(0, 0.9), 
                                                    ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                    dim = "1 longitudinal", 
                                                    prob = "20% missing", 
                                                    dgm_arg = "DGM = TDCM", 
                                                    sample = "n = 1000", 
                                                    legend = "right")
high_stability_prev_tdcm 
low_stability_prev_tdcm <- stability_plot_prev(stability_data, 
                                                   ylimits = c(0, 0.9), 
                                                   ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                   dim = "1 longitudinal", 
                                                   prob = "80% missing", 
                                                   dgm_arg = "DGM = TDCM", 
                                                   sample = "n = 200", 
                                                   legend = "right")
low_stability_prev_tdcm



###################
### Consistency ###
###################

high_consistencyC_prev_jm <- consistency_plot_prev(data = performance_data,
                                                       stat = "c_index",
                                                       dim = "1 longitudinal",
                                                       prob = "20% missing",
                                                       dgm_arg = "DGM = JM",
                                                       sample = "n = 1000",
                                                       ylabel = "Harrell's C-index range\n across landmark times",
                                                       ylimits = c(0, 0.1),
                                                       legend = "right")
high_consistencyC_prev_jm
low_consistencyC_prev_jm <- consistency_plot_prev(data = performance_data,
                                                      stat = "c_index",
                                                      dim = "1 longitudinal",
                                                      prob = "80% missing",
                                                      dgm_arg = "DGM = JM",
                                                      sample = "n = 200",
                                                      ylabel = "Harrell's C-index range\n across landmark times",
                                                      ylimits = c(0, 0.07),
                                                      legend = "right")
low_consistencyC_prev_jm
high_consistencyC_prev_tdcm <- consistency_plot_prev(data = performance_data,
                                                         stat = "c_index",
                                                         dim = "1 longitudinal",
                                                         prob = "20% missing",
                                                         dgm_arg = "DGM = TDCM",
                                                         sample = "n = 1000",
                                                         ylabel = "Harrell's C-index range\n across landmark times",
                                                         ylimits = c(0, 0.07),
                                                         legend = "right")
high_consistencyC_prev_tdcm
low_consistencyC_prev_tdcm <- consistency_plot_prev(data = performance_data,
                                                        stat = "c_index",
                                                        dim = "1 longitudinal",
                                                        prob = "80% missing",
                                                        dgm_arg = "DGM = TDCM",
                                                        sample = "n = 200",
                                                        ylabel = "Harrell's C-index range\n across landmark times",
                                                        ylimits = c(0, 0.07),
                                                        legend = "right")
low_consistencyC_prev_tdcm



high_consistencyB_prev_jm <- consistency_plot_prev(data = performance_data,
                                                       stat = "brier",
                                                       dim = "1 longitudinal",
                                                       prob = "20% missing",
                                                       dgm_arg = "DGM = JM",
                                                       sample = "n = 1000",
                                                       ylabel = "Brier score range\n across landmark times",
                                                       ylimits = c(0, 0.1),
                                                       legend = "right")
high_consistencyB_prev_jm
low_consistencyB_prev_jm <- consistency_plot_prev(data = performance_data,
                                                      stat = "brier",
                                                      dim = "1 longitudinal",
                                                      prob = "80% missing",
                                                      dgm_arg = "DGM = JM",
                                                      sample = "n = 200",
                                                      ylabel = "Brier score range\n across landmark times",
                                                      ylimits = c(0, 0.1),
                                                      legend = "right")
low_consistencyB_prev_jm
high_consistencyB_prev_tdcm <- consistency_plot_prev(data = performance_data,
                                                         stat = "brier",
                                                         dim = "1 longitudinal",
                                                         prob = "20% missing",
                                                         dgm_arg = "DGM = TDCM",
                                                         sample = "n = 1000",
                                                         ylabel = "Brier score range\n across landmark times",
                                                         ylimits = c(0, 0.1),
                                                         legend = "right")
high_consistencyB_prev_tdcm
low_consistencyB_prev_tdcm <- consistency_plot_prev(data = performance_data,
                                                        stat = "brier",
                                                        dim = "1 longitudinal",
                                                        prob = "80% missing",
                                                        dgm_arg = "DGM = TDCM",
                                                        sample = "n = 200",
                                                        ylabel = "Brier score range\n across landmark times",
                                                        ylimits = c(0, 0.1),
                                                        legend = "right")
low_consistencyB_prev_tdcm

#####################
### Fitting times ###
#####################

high_fitting_prev_jm <- fitting_time_prev(data = performance_data,
                                              stat = "fit_time", 
                                              ylimits = c(0,5), 
                                              ylabel= "Model fitting time (minutes)", 
                                              dim = "1 longitudinal", 
                                              prob = "20% missing", 
                                              dgm_arg = "DGM = JM", 
                                              sample = "n = 1000",
                                              legend = "none")
high_fitting_prev_jm
low_fitting_prev_jm <- fitting_time_prev(data = performance_data,
                                             stat = "fit_time", 
                                             ylimits = c(0,2), 
                                             ylabel= "Model fitting time (minutes)", 
                                             dim = "1 longitudinal", 
                                             prob = "80% missing",
                                             dgm_arg = "DGM = JM",
                                             sample = "n = 200",
                                             legend = "right")
low_fitting_prev_jm
high_fitting_prev_tdcm <- fitting_time_prev(data = performance_data,
                                                stat = "fit_time", 
                                                ylimits = c(0,5), 
                                                ylabel= "Model fitting time (minutes)", 
                                                dim = "1 longitudinal", 
                                                prob = "20% missing", 
                                                dgm_arg = "DGM = TDCM", 
                                                sample = "n = 1000",
                                                legend = "none")
high_fitting_prev_tdcm
low_fitting_prev_tdcm <- fitting_time_prev(data = performance_data,
                                               stat = "fit_time", 
                                               ylimits = c(0,20), 
                                               ylabel= "Model fitting time (minutes)", 
                                               dim = "1 longitudinal", 
                                               prob = "80% missing",
                                               dgm_arg = "DGM = TDCM",
                                               sample = "n = 200",
                                               legend = "right")
low_fitting_prev_tdcm


########################
### Prediction times ###
########################

high_pred_time_prev_jm <- performance_plot_prev(data = performance_data,
                                                    stat = "pred_time", 
                                                    ylimits = c(0,10), 
                                                    ylabel= "Prediction time (minutes)", 
                                                    dim = "1 longitudinal", 
                                                    prob = "20% missing", 
                                                    dgm_arg = "DGM = JM", 
                                                    sample = "n = 1000",
                                                    legend = "right")
high_pred_time_prev_jm
low_pred_time_prev_jm <- performance_plot_prev(data = performance_data,
                                                   stat = "pred_time", 
                                                   ylimits = c(0,10), 
                                                   ylabel= "Prediction time (minutes)", 
                                                   dim = "1 longitudinal", 
                                                   prob = "80% missing", 
                                                   dgm_arg = "DGM = JM", 
                                                   sample = "n = 200",
                                                   legend = "right")
low_pred_time_prev_jm 
high_pred_time_prev_tdcm <- performance_plot_prev(data = performance_data,
                                                      stat = "pred_time", 
                                                      ylimits = c(0,10), 
                                                      ylabel= "Prediction time (minutes)", 
                                                      dim = "1 longitudinal", 
                                                      prob = "20% missing", 
                                                      dgm_arg = "DGM = TDCM", 
                                                      sample = "n = 1000",
                                                      legend = "right")
high_pred_time_prev_tdcm
low_pred_time_prev_tdcm <- performance_plot_prev(data = performance_data,
                                                     stat = "pred_time", 
                                                     ylimits = c(0,30), 
                                                     ylabel= "Prediction time (minutes)", 
                                                     dim = "1 longitudinal", 
                                                     prob = "80% missing", 
                                                     dgm_arg = "DGM = TDCM", 
                                                     sample = "n = 200",
                                                     legend = "right")
low_pred_time_prev_tdcm 


#######################################
### Follow-up missingness influence ###
#######################################
#*************************************#

#######################
####### C-index #######
#######################


high_c_prob_jm <- performance_plot_prob(data = performance_data,
                                        stat = "c_index", 
                                        ylimits = c(0.5,0.9), 
                                        ylabel= "Harrell's C-index (95% quantile range)", 
                                        dim = "1 longitudinal", 
                                        prev = "prev ~ 40%", 
                                        dgm_arg = "DGM = JM", 
                                        sample = "n = 1000",
                                        legend = "right")
high_c_prob_jm

low_c_prob_jm <- performance_plot_prob(data = performance_data, 
                                       stat = "c_index", 
                                       ylimits = c(0.5,0.9), 
                                       ylabel= "Harrell's C-index (95% quantile range)", 
                                       dim = "1 longitudinal", 
                                       prev = "prev ~ 10%",
                                       dgm_arg = "DGM = JM",
                                       sample = "n = 200",
                                       legend = "right")
low_c_prob_jm

high_c_prob_tdcm <- performance_plot_prob(data = performance_data,
                                          stat = "c_index", 
                                          ylimits = c(0.5,0.9), 
                                          ylabel= "Harrell's C-index (95% quantile range)", 
                                          dim = "1 longitudinal", 
                                          prev = "prev ~ 40%",
                                          dgm_arg = "DGM = TDCM", 
                                          sample = "n = 1000",
                                          legend = "right")
high_c_prob_tdcm

low_c_prob_tdcm <- performance_plot_prob(data = performance_data, 
                                         stat = "c_index", 
                                         ylimits = c(0.4,0.9), 
                                         ylabel= "Harrell's C-index (95% quantile range)", 
                                         dim = "1 longitudinal", 
                                         prev = "prev ~ 10%",
                                         dgm_arg = "DGM = TDCM",
                                         sample = "n = 200",
                                         legend = "right")
low_c_prob_tdcm

#######################
####### Brier   #######
#######################

high_brier_prob_jm <- performance_plot_prob(data = performance_data,
                                            stat = "brier", 
                                            ylimits = c(0, 0.175), 
                                            ylabel= "Survival Brier score (95% quantile range)", 
                                            dim = "1 longitudinal", 
                                            prev = "prev ~ 40%", 
                                            dgm_arg = "DGM = JM",
                                            sample = "n = 1000",
                                            legend = "right")
high_brier_prob_jm

low_brier_prob_jm <- performance_plot_prob(data = performance_data,
                                           stat = "brier", 
                                           ylimits = c(0, 0.05), 
                                           ylabel= "Survival Brier score (95% quantile range)", 
                                           dim = "1 longitudinal", 
                                           prev = "prev ~ 10%", 
                                           dgm_arg = "DGM = JM",
                                           sample = "n = 200",
                                           legend = "right")
low_brier_prob_jm

high_brier_prob_tdcm <- performance_plot_prob(data = performance_data,
                                              stat = "brier", 
                                              ylimits = c(0, 0.175), 
                                              ylabel= "Survival Brier score (95% quantile range)", 
                                              dim = "1 longitudinal", 
                                              prev = "prev ~ 40%", 
                                              dgm_arg = "DGM = TDCM",
                                              sample = "n = 1000",
                                              legend = "right")
high_brier_prob_tdcm

low_brier_prob_tdcm <- performance_plot_prob(data = performance_data,
                                             stat = "brier", 
                                             ylimits = c(0, 0.05), 
                                             ylabel= "Survival Brier score (95% quantile range)", 
                                             dim = "1 longitudinal", 
                                             prev = "prev ~ 10%",
                                             dgm_arg = "DGM = TDCM",
                                             sample = "n = 200",
                                             legend = "right")
low_brier_prob_tdcm


#######################
####### IPA     #######
#######################

high_ipa_prob_jm <- performance_plot_prob(data = brier_ipa,
                                          stat = "ipa", 
                                          ylimits = c(-25, 75), 
                                          ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                          dim = "1 longitudinal", 
                                          prev = "prev ~ 40%", 
                                          dgm_arg = "DGM = JM", 
                                          sample = "n = 1000",
                                          legend = "right")
high_ipa_prob_jm
low_ipa_prob_jm <- performance_plot_prob(data = brier_ipa,
                                         stat = "ipa", 
                                         ylimits = c(-25, 75), 
                                         ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                         dim = "1 longitudinal", 
                                         prev = "prev ~ 10%", 
                                         dgm_arg = "DGM = JM",  
                                         sample = "n = 200",
                                         legend = "right")
low_ipa_prob_jm
high_ipa_prob_tdcm <- performance_plot_prob(data = brier_ipa,
                                            stat = "ipa", 
                                            ylimits = c(-25, 75), 
                                            ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                            dim = "1 longitudinal", 
                                            prev = "prev ~ 40%",
                                            dgm_arg = "DGM = TDCM", 
                                            sample = "n = 1000",
                                            legend = "right")
high_ipa_prob_tdcm
low_ipa_prob_tdcm <- performance_plot_prob(data = brier_ipa,
                                           stat = "ipa", 
                                           ylimits = c(-25, 75), 
                                           ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                           dim = "1 longitudinal", 
                                           prev = "prev ~ 10%", 
                                           dgm_arg = "DGM = TDCM",  
                                           sample = "n = 200",
                                           legend = "right")
low_ipa_prob_tdcm


##########################
### Optimism - C-index ###
##########################

high_optC_prob_jm <- performance_plot_prob(data = performance_data,
                                           stat = "c_index_adjusted", 
                                           ylimits = c(-0.2, 0.2), 
                                           ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                           dim = "1 longitudinal", 
                                           prev = "prev ~ 40%", 
                                           dgm_arg = "DGM = JM", 
                                           sample = "n = 1000",
                                           legend = "right")
high_optC_prob_jm
low_optC_prob_jm <- performance_plot_prob(data = performance_data,
                                          stat = "c_index_adjusted", 
                                          ylimits = c(-0.25, 0.5), 
                                          ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                          dim = "1 longitudinal", 
                                          prev = "prev ~ 10%", 
                                          dgm_arg = "DGM = JM",  
                                          sample = "n = 200",
                                          legend = "right")
low_optC_prob_jm
high_optC_prob_tdcm <- performance_plot_prob(data = performance_data,
                                             stat = "c_index_adjusted", 
                                             ylimits = c(-0.2, 0.2), 
                                             ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                             dim = "1 longitudinal", 
                                             prev = "prev ~ 40%", 
                                             dgm_arg = "DGM = TDCM", 
                                             sample = "n = 1000",
                                             legend = "right")
high_optC_prob_tdcm
low_optC_prob_tdcm <- performance_plot_prob(data = performance_data,
                                            stat = "c_index_adjusted", 
                                            ylimits = c(-0.25, 0.5), 
                                            ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                            dim = "1 longitudinal", 
                                            prev = "prev ~ 10%", 
                                            dgm_arg = "DGM = TDCM",  
                                            sample = "n = 200",
                                            legend = "right")
low_optC_prob_tdcm 

##########################
### Optimism - Brier   ###
##########################

high_optB_prob_jm <- performance_plot_prob(data = performance_data,
                                           stat = "brier_adjusted", 
                                           ylimits = c(-0.07, 0.07), 
                                           ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                           dim = "1 longitudinal", 
                                           prev = "prev ~ 40%", 
                                           dgm_arg = "DGM = JM",  
                                           sample = "n = 1000",
                                           legend = "right")
high_optB_prob_jm
low_optB_prob_jm <- performance_plot_prob(data = performance_data,
                                          stat = "brier_adjusted", 
                                          ylimits = c(-0.07, 0.07), 
                                          ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                          dim = "1 longitudinal", 
                                          prev = "prev ~ 10%", 
                                          dgm_arg = "DGM = JM", 
                                          sample = "n = 200",
                                          legend = "right")

low_optB_prob_jm
high_optB_prob_tdcm <- performance_plot_prob(data = performance_data,
                                             stat = "brier_adjusted", 
                                             ylimits = c(-0.07, 0.07), 
                                             ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                             dim = "1 longitudinal", 
                                             prev = "prev ~ 40%",
                                             dgm_arg = "DGM = TDCM",  
                                             sample = "n = 1000",
                                             legend = "right")
high_optB_prob_tdcm
low_optB_prob_tdcm <- performance_plot_prob(data = performance_data,
                                            stat = "brier_adjusted", 
                                            ylimits = c(-0.07, 0.07), 
                                            ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                            dim = "1 longitudinal", 
                                            prev = "prev ~ 10%", 
                                            dgm_arg = "DGM = TDCM", 
                                            sample = "n = 200",
                                            legend = "right")

low_optB_prob_tdcm


#####################
#### Stability   ####
#####################

high_stability_prob_jm <- stability_plot_prob(stability_data, 
                                              ylimits = c(0, 0.9), 
                                              ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                              dim = "1 longitudinal", 
                                              prev = "prev ~ 40%", 
                                              dgm_arg = "DGM = JM", 
                                              sample = "n = 1000", 
                                              legend = "right")
high_stability_prob_jm
low_stability_prob_jm <- stability_plot_prob(stability_data, 
                                             ylimits = c(0, 0.9), 
                                             ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                             dim = "1 longitudinal", 
                                             prev = "prev ~ 10%", 
                                             dgm_arg = "DGM = JM", 
                                             sample = "n = 200", 
                                             legend = "right")
low_stability_prob_jm
high_stability_prob_tdcm <- stability_plot_prob(stability_data, 
                                                ylimits = c(0, 0.9), 
                                                ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                dim = "1 longitudinal", 
                                                prev = "prev ~ 40%", 
                                                dgm_arg = "DGM = TDCM", 
                                                sample = "n = 1000", 
                                                legend = "right")
high_stability_prob_tdcm 
low_stability_prob_tdcm <- stability_plot_prob(stability_data, 
                                               ylimits = c(0, 0.9), 
                                               ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                               dim = "1 longitudinal", 
                                               prev = "prev ~ 10%", 
                                               dgm_arg = "DGM = TDCM", 
                                               sample = "n = 200", 
                                               legend = "right")
low_stability_prob_tdcm



###################
### Consistency ###
###################

high_consistencyC_prob_jm <- consistency_plot_prob(data = performance_data,
                                                   stat = "c_index",
                                                   dim = "1 longitudinal",
                                                   prev = "prev ~ 40%",
                                                   dgm_arg = "DGM = JM",
                                                   sample = "n = 1000",
                                                   ylabel = "Harrell's C-index range\n across landmark times",
                                                   ylimits = c(0, 0.1),
                                                   legend = "right")
high_consistencyC_prob_jm
low_consistencyC_prob_jm <- consistency_plot_prob(data = performance_data,
                                                  stat = "c_index",
                                                  dim = "1 longitudinal",
                                                  prev = "prev ~ 10%",
                                                  dgm_arg = "DGM = JM",
                                                  sample = "n = 200",
                                                  ylabel = "Harrell's C-index range\n across landmark times",
                                                  ylimits = c(0, 0.07),
                                                  legend = "right")
low_consistencyC_prob_jm
high_consistencyC_prob_tdcm <- consistency_plot_prob(data = performance_data,
                                                     stat = "c_index",
                                                     dim = "1 longitudinal",
                                                     prev = "prev ~ 40%",
                                                     dgm_arg = "DGM = TDCM",
                                                     sample = "n = 1000",
                                                     ylabel = "Harrell's C-index range\n across landmark times",
                                                     ylimits = c(0, 0.07),
                                                     legend = "right")
high_consistencyC_prob_tdcm
low_consistencyC_prob_tdcm <- consistency_plot_prob(data = performance_data,
                                                    stat = "c_index",
                                                    dim = "1 longitudinal",
                                                    prev = "prev ~ 10%",
                                                    dgm_arg = "DGM = TDCM",
                                                    sample = "n = 200",
                                                    ylabel = "Harrell's C-index range\n across landmark times",
                                                    ylimits = c(0, 0.07),
                                                    legend = "right")
low_consistencyC_prob_tdcm



high_consistencyB_prob_jm <- consistency_plot_prob(data = performance_data,
                                                   stat = "brier",
                                                   dim = "1 longitudinal",
                                                   prev = "prev ~ 40%",
                                                   dgm_arg = "DGM = JM",
                                                   sample = "n = 1000",
                                                   ylabel = "Brier score range\n across landmark times",
                                                   ylimits = c(0, 0.1),
                                                   legend = "right")
high_consistencyB_prob_jm
low_consistencyB_prob_jm <- consistency_plot_prob(data = performance_data,
                                                  stat = "brier",
                                                  dim = "1 longitudinal",
                                                  prev = "prev ~ 10%",
                                                  dgm_arg = "DGM = JM",
                                                  sample = "n = 200",
                                                  ylabel = "Brier score range\n across landmark times",
                                                  ylimits = c(0, 0.1),
                                                  legend = "right")
low_consistencyB_prob_jm
high_consistencyB_prob_tdcm <- consistency_plot_prob(data = performance_data,
                                                     stat = "brier",
                                                     dim = "1 longitudinal",
                                                     prev = "prev ~ 40%",
                                                     dgm_arg = "DGM = TDCM",
                                                     sample = "n = 1000",
                                                     ylabel = "Brier score range\n across landmark times",
                                                     ylimits = c(0, 0.1),
                                                     legend = "right")
high_consistencyB_prob_tdcm
low_consistencyB_prob_tdcm <- consistency_plot_prob(data = performance_data,
                                                    stat = "brier",
                                                    dim = "1 longitudinal",
                                                    prev = "prev ~ 10%",
                                                    dgm_arg = "DGM = TDCM",
                                                    sample = "n = 200",
                                                    ylabel = "Brier score range\n across landmark times",
                                                    ylimits = c(0, 0.1),
                                                    legend = "right")
low_consistencyB_prob_tdcm

#####################
### Fitting times ###
#####################

high_fitting_prob_jm <- fitting_time_prob(data = performance_data,
                                          stat = "fit_time", 
                                          ylimits = c(0,5), 
                                          ylabel= "Model fitting time (minutes)", 
                                          dim = "1 longitudinal", 
                                          prev = "prev ~ 40%", 
                                          dgm_arg = "DGM = JM", 
                                          sample = "n = 1000",
                                          legend = "right")
high_fitting_prob_jm
low_fitting_prob_jm <- fitting_time_prob(data = performance_data,
                                         stat = "fit_time", 
                                         ylimits = c(0,5), 
                                         ylabel= "Model fitting time (minutes)", 
                                         dim = "1 longitudinal", 
                                         prev = "prev ~ 10%",
                                         dgm_arg = "DGM = JM",
                                         sample = "n = 200",
                                         legend = "right")
low_fitting_prob_jm
high_fitting_prob_tdcm <- fitting_time_prob(data = performance_data,
                                            stat = "fit_time", 
                                            ylimits = c(0,5), 
                                            ylabel= "Model fitting time (minutes)", 
                                            dim = "1 longitudinal", 
                                            prev = "prev ~ 40%", 
                                            dgm_arg = "DGM = TDCM", 
                                            sample = "n = 1000",
                                            legend = "right")
high_fitting_prob_tdcm
low_fitting_prob_tdcm <- fitting_time_prob(data = performance_data,
                                           stat = "fit_time", 
                                           ylimits = c(0,20), 
                                           ylabel= "Model fitting time (minutes)", 
                                           dim = "1 longitudinal", 
                                           prev = "prev ~ 10%",
                                           dgm_arg = "DGM = TDCM",
                                           sample = "n = 200",
                                           legend = "right")
low_fitting_prob_tdcm


########################
### Prediction times ###
########################

high_pred_time_prob_jm <- performance_plot_prob(data = performance_data,
                                                stat = "pred_time", 
                                                ylimits = c(0,10), 
                                                ylabel= "Prediction time (minutes)", 
                                                dim = "1 longitudinal", 
                                                prev = "prev ~ 40%", 
                                                dgm_arg = "DGM = JM", 
                                                sample = "n = 1000",
                                                legend = "right")
high_pred_time_prob_jm
low_pred_time_prob_jm <- performance_plot_prob(data = performance_data,
                                               stat = "pred_time", 
                                               ylimits = c(0,10), 
                                               ylabel= "Prediction time (minutes)", 
                                               dim = "1 longitudinal", 
                                               prev = "prev ~ 10%", 
                                               dgm_arg = "DGM = JM", 
                                               sample = "n = 200",
                                               legend = "right")
low_pred_time_prob_jm 
high_pred_time_prob_tdcm <- performance_plot_prob(data = performance_data,
                                                  stat = "pred_time", 
                                                  ylimits = c(0,10), 
                                                  ylabel= "Prediction time (minutes)", 
                                                  dim = "1 longitudinal", 
                                                  prev = "prev ~ 40%", 
                                                  dgm_arg = "DGM = TDCM", 
                                                  sample = "n = 1000",
                                                  legend = "right")
high_pred_time_prob_tdcm
low_pred_time_prob_tdcm <- performance_plot_prob(data = performance_data,
                                                 stat = "pred_time", 
                                                 ylimits = c(0,30), 
                                                 ylabel= "Prediction time (minutes)", 
                                                 dim = "1 longitudinal", 
                                                 prev = "prev ~ 10%", 
                                                 dgm_arg = "DGM = TDCM", 
                                                 sample = "n = 200",
                                                 legend = "right")
low_pred_time_prob_tdcm 

#################################
# Three longitudinal predictors #
#################################

##########################
### Comparing extremes ###
##########################
#************************#

#######################
####### C-index #######
#######################


high_extreme_c <- performance_plot_dgm(data = performance_data,
                                       stat = "c_index", 
                                       ylimits = c(0.4,0.9), 
                                       ylabel= " ", 
                                       dim = "3 longitudinal", 
                                       prob = "20% missing", 
                                       sample = "n = 1000", 
                                       prev = "prev ~ 40%",
                                       legend = "right")
low_extreme_c <- performance_plot_dgm(data = performance_data, 
                                      stat = "c_index", 
                                      ylimits = c(0.4,0.9), 
                                      ylabel= "Harrell's C-index (95% quantile range)", 
                                      dim = "3 longitudinal", 
                                      prob = "80% missing", 
                                      sample = "n = 200", 
                                      prev = "prev ~ 10%",
                                      legend = "none")

low_extreme_c | high_extreme_c

#######################
####### Brier   #######
#######################
high_extreme_brier <- performance_plot_dgm(data = performance_data,
                                           stat = "brier", 
                                           ylimits = c(0, 0.175), 
                                           ylabel= " ", 
                                           dim = "3 longitudinal", 
                                           prob = "20% missing", 
                                           sample = "n = 1000", 
                                           prev = "prev ~ 40%",
                                           legend = "right")
low_extreme_brier <- performance_plot_dgm(data = performance_data,
                                          stat = "brier", 
                                          ylimits = c(0, 0.175), 
                                          ylabel= "Survival Brier score (95% quantile range)", 
                                          dim = "3 longitudinal", 
                                          prob = "80% missing", 
                                          sample = "n = 200", 
                                          prev = "prev ~ 10%",
                                          legend = "none")

low_extreme_brier | high_extreme_brier

#######################
####### IPA     #######
#######################

high_extreme_ipa <- performance_plot_dgm(data = brier_ipa,
                                         stat = "ipa", 
                                         ylimits = c(-25, 100), 
                                         ylabel= " ", 
                                         dim = "3 longitudinal", 
                                         prob = "20% missing", 
                                         sample = "n = 1000", 
                                         prev = "prev ~ 40%",
                                         legend = "right")
low_extreme_ipa <- performance_plot_dgm(data = brier_ipa,
                                        stat = "ipa", 
                                        ylimits = c(-25, 100), 
                                        ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                        dim = "3 longitudinal", 
                                        prob = "80% missing", 
                                        sample = "n = 200", 
                                        prev = "prev ~ 10%",
                                        legend = "none")

low_extreme_ipa | high_extreme_ipa

##########################
### Optimism - C-index ###
##########################

high_extreme_optC <- performance_plot_dgm(data = performance_data,
                                          stat = "c_index_adjusted", 
                                          ylimits = c(-0.25, 0.5), 
                                          ylabel= " ", 
                                          dim = "3 longitudinal", 
                                          prob = "20% missing", 
                                          sample = "n = 1000", 
                                          prev = "prev ~ 40%",
                                          legend = "right")

low_extreme_optC <- performance_plot_dgm(data = performance_data,
                                         stat = "c_index_adjusted", 
                                         ylimits = c(-0.25, 0.5), 
                                         ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                         dim = "3 longitudinal", 
                                         prob = "80% missing", 
                                         sample = "n = 200", 
                                         prev = "prev ~ 10%",
                                         legend = "none")

low_extreme_optC | high_extreme_optC

##########################
### Optimism - Brier   ###
##########################

high_extreme_optB <- performance_plot_dgm(data = performance_data,
                                          stat = "brier_adjusted", 
                                          ylimits = c(-0.03, 0.03), 
                                          ylabel= " ", 
                                          dim = "3 longitudinal", 
                                          prob = "20% missing", 
                                          sample = "n = 1000", 
                                          prev = "prev ~ 40%",
                                          legend = "right")

low_extreme_optB <- performance_plot_dgm(data = performance_data,
                                         stat = "brier_adjusted", 
                                         ylimits = c(-0.03, 0.03), 
                                         ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                         dim = "3 longitudinal", 
                                         prob = "80% missing", 
                                         sample = "n = 200", 
                                         prev = "prev ~ 10%",
                                         legend = "none")

low_extreme_optB | high_extreme_optB

#####################
#### Stability   ####
#####################

trial <- stability_data %>% filter(!is.na(landmark_time),
                                   dimension == "3 longitudinal",
                                   probability_missing == "80% missing",
                                   n == "n = 200",
                                   event_prevalence == "prev ~ 10%",
                                   model == 2) %>% group_by(dgm, model, 
                                                                          landmark_time) %>% 
  summarise(value_up = quantile(range, 0.975, na.rm = TRUE), 
            value_down = quantile(range, 0.025, na.rm =TRUE),
            value = mean(range, na.rm=TRUE))

high_extreme_stability <- stability_plot_dgm(stability_data, 
                                             ylimits = c(0, 1), 
                                             ylabel = " ", 
                                             dim = "3 longitudinal", 
                                             prob = "20% missing", 
                                             sample = "n = 1000", 
                                             prev = "prev ~ 40%", 
                                             legend = "right")

low_extreme_stability <- stability_plot_dgm(stability_data, 
                                            ylimits = c(0, 1), 
                                            ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                            dim = "3 longitudinal", 
                                            prob = "20% missing", 
                                            sample = "n = 200", 
                                            prev = "prev ~ 10%", 
                                            legend = "none")

low_extreme_stability | high_extreme_stability


###################
### Consistency ###
###################

high_extreme_consistencyC <- consistency_plot_dgm(data = performance_data,
                                                  stat = "c_index",
                                                  dim = "3 longitudinal",
                                                  prob = "20% missing",
                                                  sample = "n = 1000",
                                                  prev = "prev ~ 40%",
                                                  ylabel = " ",
                                                  ylimits = c(0, 0.07),
                                                  legend = "none")
low_extreme_consistencyC <- consistency_plot_dgm(data = performance_data,
                                                 stat = "c_index",
                                                 dim = "3 longitudinal",
                                                 prob = "80% missing",
                                                 sample = "n = 200",
                                                 prev = "prev ~ 10%",
                                                 ylabel = "Harrell's C-index range\n across landmark times",
                                                 ylimits = c(0, 0.07),
                                                 legend = "none")

low_extreme_consistencyC | high_extreme_consistencyC



high_extreme_consistencyB <- consistency_plot_dgm(data = performance_data,
                                                  stat = "brier",
                                                  dim = "3 longitudinal",
                                                  prob = "20% missing",
                                                  sample = "n = 1000",
                                                  prev = "prev ~ 40%",
                                                  ylabel = " ",
                                                  ylimits = c(0, 0.1),
                                                  legend = "none")
low_extreme_consistencyB <- consistency_plot_dgm(data = performance_data,
                                                 stat = "brier",
                                                 dim = "3 longitudinal",
                                                 prob = "80% missing",
                                                 sample = "n = 200",
                                                 prev = "prev ~ 10%",
                                                 ylabel = "Harrell's C-index range\n across landmark times",
                                                 ylimits = c(0, 0.1),
                                                 legend = "none")

low_extreme_consistencyB | high_extreme_consistencyB


########################
### Prediction times ###
########################

high_extreme_pred_time <- performance_plot_dgm(data = performance_data,
                                               stat = "pred_time", 
                                               ylimits = c(0,50), 
                                               ylabel= " ", 
                                               dim = "3 longitudinal", 
                                               prob = "20% missing", 
                                               sample = "n = 1000", 
                                               prev = "prev ~ 40%",
                                               legend = "right")

low_extreme_pred_time <- performance_plot_dgm(data = performance_data,
                                              stat = "pred_time", 
                                              ylimits = c(0,50), 
                                              ylabel= "Prediction time (minutes)", 
                                              dim = "3 longitudinal", 
                                              prob = "80% missing", 
                                              sample = "n = 200", 
                                              prev = "prev ~ 10%",
                                              legend = "none")
low_extreme_pred_time | high_extreme_pred_time


#############################
### Sample size influence ###
#############################
#***************************#

#######################
####### C-index #######
#######################


high_c_sample_jm <- performance_plot_sample(data = performance_data,
                                            stat = "c_index", 
                                            ylimits = c(0.5,0.9), 
                                            ylabel= "Harrell's C-index (95% quantile range)", 
                                            dim = "3 longitudinal", 
                                            prob = "20% missing", 
                                            dgm_arg = "DGM = JM", 
                                            prev = "prev ~ 40%",
                                            legend = "right")
high_c_sample_jm

low_c_sample_jm <- performance_plot_sample(data = performance_data, 
                                           stat = "c_index", 
                                           ylimits = c(0.4,0.9), 
                                           ylabel= "Harrell's C-index (95% quantile range)", 
                                           dim = "3 longitudinal", 
                                           prob = "80% missing", 
                                           dgm_arg = "DGM = JM",
                                           prev = "prev ~ 10%",
                                           legend = "right")
low_c_sample_jm

high_c_sample_tdcm <- performance_plot_sample(data = performance_data,
                                              stat = "c_index", 
                                              ylimits = c(0.5,0.9), 
                                              ylabel= "Harrell's C-index (95% quantile range)", 
                                              dim = "3 longitudinal", 
                                              prob = "20% missing", 
                                              dgm_arg = "DGM = TDCM", 
                                              prev = "prev ~ 40%",
                                              legend = "right")
high_c_sample_tdcm

low_c_sample_tdcm <- performance_plot_sample(data = performance_data, 
                                             stat = "c_index", 
                                             ylimits = c(0.4,0.9), 
                                             ylabel= "Harrell's C-index (95% quantile range)", 
                                             dim = "3 longitudinal", 
                                             prob = "80% missing", 
                                             dgm_arg = "DGM = TDCM",
                                             prev = "prev ~ 10%",
                                             legend = "right")
low_c_sample_tdcm

#######################
####### Brier   #######
#######################

high_brier_sample_jm <- performance_plot_sample(data = performance_data,
                                                stat = "brier", 
                                                ylimits = c(0, 0.175), 
                                                ylabel= "Survival Brier score (95% quantile range)", 
                                                dim = "3 longitudinal", 
                                                prob = "20% missing", 
                                                dgm_arg = "DGM = JM",
                                                prev = "prev ~ 40%",
                                                legend = "right")
high_brier_sample_jm

low_brier_sample_jm <- performance_plot_sample(data = performance_data,
                                               stat = "brier", 
                                               ylimits = c(0, 0.05), 
                                               ylabel= "Survival Brier score (95% quantile range)", 
                                               dim = "3 longitudinal", 
                                               prob = "80% missing", 
                                               dgm_arg = "DGM = JM",
                                               prev = "prev ~ 10%",
                                               legend = "right")
low_brier_sample_jm

high_brier_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                  stat = "brier", 
                                                  ylimits = c(0, 0.175), 
                                                  ylabel= "Survival Brier score (95% quantile range)", 
                                                  dim = "3 longitudinal", 
                                                  prob = "20% missing", 
                                                  dgm_arg = "DGM = TDCM",
                                                  prev = "prev ~ 40%",
                                                  legend = "right")
high_brier_sample_tdcm

low_brier_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                 stat = "brier", 
                                                 ylimits = c(0, 0.05), 
                                                 ylabel= "Survival Brier score (95% quantile range)", 
                                                 dim = "3 longitudinal", 
                                                 prob = "80% missing", 
                                                 dgm_arg = "DGM = TDCM",
                                                 prev = "prev ~ 10%",
                                                 legend = "right")
low_brier_sample_tdcm


#######################
####### IPA     #######
#######################

high_ipa_sample_jm <- performance_plot_sample(data = brier_ipa,
                                              stat = "ipa", 
                                              ylimits = c(-25, 75), 
                                              ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                              dim = "3 longitudinal", 
                                              prob = "20% missing", 
                                              dgm_arg = "DGM = JM", 
                                              prev = "prev ~ 40%",
                                              legend = "right")
high_ipa_sample_jm
low_ipa_sample_jm <- performance_plot_sample(data = brier_ipa,
                                             stat = "ipa", 
                                             ylimits = c(-25, 100), 
                                             ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                             dim = "3 longitudinal", 
                                             prob = "80% missing", 
                                             dgm_arg = "DGM = JM",  
                                             prev = "prev ~ 10%",
                                             legend = "right")
low_ipa_sample_jm
high_ipa_sample_tdcm <- performance_plot_sample(data = brier_ipa,
                                                stat = "ipa", 
                                                ylimits = c(-25, 75), 
                                                ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                                dim = "3 longitudinal", 
                                                prob = "20% missing", 
                                                dgm_arg = "DGM = TDCM", 
                                                prev = "prev ~ 40%",
                                                legend = "right")
high_ipa_sample_tdcm
low_ipa_sample_tdcm <- performance_plot_sample(data = brier_ipa,
                                               stat = "ipa", 
                                               ylimits = c(-25, 75), 
                                               ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                               dim = "3 longitudinal", 
                                               prob = "80% missing", 
                                               dgm_arg = "DGM = TDCM",  
                                               prev = "prev ~ 10%",
                                               legend = "right")
low_ipa_sample_tdcm


##########################
### Optimism - C-index ###
##########################

high_optC_sample_jm <- performance_plot_sample(data = performance_data,
                                               stat = "c_index_adjusted", 
                                               ylimits = c(-0.2, 0.2), 
                                               ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                               dim = "3 longitudinal", 
                                               prob = "20% missing", 
                                               dgm_arg = "DGM = JM", 
                                               prev = "prev ~ 40%",
                                               legend = "right")
high_optC_sample_jm
low_optC_sample_jm <- performance_plot_sample(data = performance_data,
                                              stat = "c_index_adjusted", 
                                              ylimits = c(-0.25, 0.5), 
                                              ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                              dim = "3 longitudinal", 
                                              prob = "80% missing", 
                                              dgm_arg = "DGM = JM",  
                                              prev = "prev ~ 10%",
                                              legend = "right")
low_optC_sample_jm
high_optC_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                 stat = "c_index_adjusted", 
                                                 ylimits = c(-0.2, 0.2), 
                                                 ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                                 dim = "3 longitudinal", 
                                                 prob = "20% missing", 
                                                 dgm_arg = "DGM = TDCM", 
                                                 prev = "prev ~ 40%",
                                                 legend = "right")
high_optC_sample_tdcm
low_optC_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                stat = "c_index_adjusted", 
                                                ylimits = c(-0.25, 0.5), 
                                                ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                                dim = "3 longitudinal", 
                                                prob = "80% missing", 
                                                dgm_arg = "DGM = TDCM",  
                                                prev = "prev ~ 10%",
                                                legend = "right")
low_optC_sample_tdcm 

##########################
### Optimism - Brier   ###
##########################

high_optB_sample_jm <- performance_plot_sample(data = performance_data,
                                               stat = "brier_adjusted", 
                                               ylimits = c(-0.07, 0.07), 
                                               ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                               dim = "3 longitudinal", 
                                               prob = "20% missing", 
                                               dgm_arg = "DGM = JM",  
                                               prev = "prev ~ 40%",
                                               legend = "right")
high_optB_sample_jm
low_optB_sample_jm <- performance_plot_sample(data = performance_data,
                                              stat = "brier_adjusted", 
                                              ylimits = c(-0.07, 0.07), 
                                              ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                              dim = "3 longitudinal", 
                                              prob = "80% missing", 
                                              dgm_arg = "DGM = JM", 
                                              prev = "prev ~ 10%",
                                              legend = "right")

low_optB_sample_jm
high_optB_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                 stat = "brier_adjusted", 
                                                 ylimits = c(-0.07, 0.07), 
                                                 ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                                 dim = "3 longitudinal", 
                                                 prob = "20% missing", 
                                                 dgm_arg = "DGM = TDCM",  
                                                 prev = "prev ~ 40%",
                                                 legend = "right")
high_optB_sample_tdcm
low_optB_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                stat = "brier_adjusted", 
                                                ylimits = c(-0.07, 0.07), 
                                                ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                                dim = "3 longitudinal", 
                                                prob = "80% missing", 
                                                dgm_arg = "DGM = TDCM", 
                                                prev = "prev ~ 10%",
                                                legend = "right")

low_optB_sample_tdcm

#####################
#### Stability   ####
#####################

high_stability_sample_jm <- stability_plot_sample(stability_data, 
                                                  ylimits = c(0, 0.9), 
                                                  ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                  dim = "3 longitudinal", 
                                                  prob = "20% missing", 
                                                  dgm_arg = "DGM = JM", 
                                                  prev = "prev ~ 40%", 
                                                  legend = "right")
high_stability_sample_jm
low_stability_sample_jm <- stability_plot_sample(stability_data, 
                                                 ylimits = c(0, 1), 
                                                 ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                 dim = "3 longitudinal", 
                                                 prob = "20% missing", 
                                                 dgm_arg = "DGM = JM", 
                                                 prev = "prev ~ 10%", 
                                                 legend = "right")
low_stability_sample_jm
high_stability_sample_tdcm <- stability_plot_sample(stability_data, 
                                                    ylimits = c(0, 0.9), 
                                                    ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                    dim = "3 longitudinal", 
                                                    prob = "20% missing", 
                                                    dgm_arg = "DGM = TDCM", 
                                                    prev = "prev ~ 40%", 
                                                    legend = "right")
high_stability_sample_tdcm 
low_stability_sample_tdcm <- stability_plot_sample(stability_data, 
                                                   ylimits = c(0, 0.9), 
                                                   ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                   dim = "3 longitudinal", 
                                                   prob = "80% missing", 
                                                   dgm_arg = "DGM = TDCM", 
                                                   prev = "prev ~ 10%", 
                                                   legend = "right")
low_stability_sample_tdcm



###################
### Consistency ###
###################

high_consistencyC_sample_jm <- consistency_plot_sample(data = performance_data,
                                                       stat = "c_index",
                                                       dim = "3 longitudinal",
                                                       prob = "20% missing",
                                                       dgm_arg = "DGM = JM",
                                                       prev = "prev ~ 40%",
                                                       ylabel = "Harrell's C-index range\n across landmark times",
                                                       ylimits = c(0, 0.1),
                                                       legend = "right")
high_consistencyC_sample_jm
low_consistencyC_sample_jm <- consistency_plot_sample(data = performance_data,
                                                      stat = "c_index",
                                                      dim = "3 longitudinal",
                                                      prob = "80% missing",
                                                      dgm_arg = "DGM = JM",
                                                      prev = "prev ~ 10%",
                                                      ylabel = "Harrell's C-index range\n across landmark times",
                                                      ylimits = c(0, 0.07),
                                                      legend = "right")
low_consistencyC_sample_jm
high_consistencyC_sample_tdcm <- consistency_plot_sample(data = performance_data,
                                                         stat = "c_index",
                                                         dim = "3 longitudinal",
                                                         prob = "20% missing",
                                                         dgm_arg = "DGM = TDCM",
                                                         prev = "prev ~ 40%",
                                                         ylabel = "Harrell's C-index range\n across landmark times",
                                                         ylimits = c(0, 0.07),
                                                         legend = "right")
high_consistencyC_sample_tdcm
low_consistencyC_sample_tdcm <- consistency_plot_sample(data = performance_data,
                                                        stat = "c_index",
                                                        dim = "3 longitudinal",
                                                        prob = "80% missing",
                                                        dgm_arg = "DGM = TDCM",
                                                        prev = "prev ~ 10%",
                                                        ylabel = "Harrell's C-index range\n across landmark times",
                                                        ylimits = c(0, 0.07),
                                                        legend = "right")
low_consistencyC_sample_tdcm



high_consistencyB_sample_jm <- consistency_plot_sample(data = performance_data,
                                                       stat = "brier",
                                                       dim = "3 longitudinal",
                                                       prob = "20% missing",
                                                       dgm_arg = "DGM = JM",
                                                       prev = "prev ~ 40%",
                                                       ylabel = "Brier score range\n across landmark times",
                                                       ylimits = c(0, 0.1),
                                                       legend = "right")
high_consistencyB_sample_jm
low_consistencyB_sample_jm <- consistency_plot_sample(data = performance_data,
                                                      stat = "brier",
                                                      dim = "3 longitudinal",
                                                      prob = "80% missing",
                                                      dgm_arg = "DGM = JM",
                                                      prev = "prev ~ 10%",
                                                      ylabel = "Brier score range\n across landmark times",
                                                      ylimits = c(0, 0.1),
                                                      legend = "right")
low_consistencyB_sample_jm
high_consistencyB_sample_tdcm <- consistency_plot_sample(data = performance_data,
                                                         stat = "brier",
                                                         dim = "3 longitudinal",
                                                         prob = "20% missing",
                                                         dgm_arg = "DGM = TDCM",
                                                         prev = "prev ~ 40%",
                                                         ylabel = "Brier score range\n across landmark times",
                                                         ylimits = c(0, 0.1),
                                                         legend = "right")
high_consistencyB_sample_tdcm
low_consistencyB_sample_tdcm <- consistency_plot_sample(data = performance_data,
                                                        stat = "brier",
                                                        dim = "3 longitudinal",
                                                        prob = "80% missing",
                                                        dgm_arg = "DGM = TDCM",
                                                        prev = "prev ~ 10%",
                                                        ylabel = "Brier score range\n across landmark times",
                                                        ylimits = c(0, 0.1),
                                                        legend = "right")
low_consistencyB_sample_tdcm

#####################
### Fitting times ###
#####################

high_fitting_sample_jm <- fitting_time_sample(data = performance_data,
                                              stat = "fit_time", 
                                              ylimits = c(0,10), 
                                              ylabel= "Model fitting time (minutes)", 
                                              dim = "3 longitudinal", 
                                              prob = "20% missing", 
                                              dgm_arg = "DGM = JM", 
                                              prev = "prev ~ 40%",
                                              legend = "none")
high_fitting_sample_jm
low_fitting_sample_jm <- fitting_time_sample(data = performance_data,
                                             stat = "fit_time", 
                                             ylimits = c(0,10), 
                                             ylabel= "Model fitting time (minutes)", 
                                             dim = "3 longitudinal", 
                                             prob = "80% missing",
                                             dgm_arg = "DGM = JM",
                                             prev = "prev ~ 10%",
                                             legend = "right")
low_fitting_sample_jm
high_fitting_sample_tdcm <- fitting_time_sample(data = performance_data,
                                                stat = "fit_time", 
                                                ylimits = c(0,10), 
                                                ylabel= "Model fitting time (minutes)", 
                                                dim = "3 longitudinal", 
                                                prob = "20% missing", 
                                                dgm_arg = "DGM = TDCM", 
                                                prev = "prev ~ 40%",
                                                legend = "none")
high_fitting_sample_tdcm
low_fitting_sample_tdcm <- fitting_time_sample(data = performance_data,
                                               stat = "fit_time", 
                                               ylimits = c(0,10), 
                                               ylabel= "Model fitting time (minutes)", 
                                               dim = "3 longitudinal", 
                                               prob = "80% missing",
                                               dgm_arg = "DGM = TDCM",
                                               prev = "prev ~ 10%",
                                               legend = "right")
low_fitting_sample_tdcm


########################
### Prediction times ###
########################

high_pred_time_sample_jm <- performance_plot_sample(data = performance_data,
                                                    stat = "pred_time", 
                                                    ylimits = c(0,50), 
                                                    ylabel= "Prediction time (minutes)", 
                                                    dim = "3 longitudinal", 
                                                    prob = "20% missing", 
                                                    dgm_arg = "DGM = JM", 
                                                    prev = "prev ~ 40%",
                                                    legend = "right")
high_pred_time_sample_jm
low_pred_time_sample_jm <- performance_plot_sample(data = performance_data,
                                                   stat = "pred_time", 
                                                   ylimits = c(0,50), 
                                                   ylabel= "Prediction time (minutes)", 
                                                   dim = "3 longitudinal", 
                                                   prob = "80% missing", 
                                                   dgm_arg = "DGM = JM", 
                                                   prev = "prev ~ 10%",
                                                   legend = "right")
low_pred_time_sample_jm 
high_pred_time_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                      stat = "pred_time", 
                                                      ylimits = c(0,50), 
                                                      ylabel= "Prediction time (minutes)", 
                                                      dim = "3 longitudinal", 
                                                      prob = "20% missing", 
                                                      dgm_arg = "DGM = TDCM", 
                                                      prev = "prev ~ 40%",
                                                      legend = "right")
high_pred_time_sample_tdcm
low_pred_time_sample_tdcm <- performance_plot_sample(data = performance_data,
                                                     stat = "pred_time", 
                                                     ylimits = c(0,50), 
                                                     ylabel= "Prediction time (minutes)", 
                                                     dim = "3 longitudinal", 
                                                     prob = "80% missing", 
                                                     dgm_arg = "DGM = TDCM", 
                                                     prev = "prev ~ 10%",
                                                     legend = "right")
low_pred_time_sample_tdcm 

##################################
### Event prevalence influence ###
##################################
#********************************#

#######################
####### C-index #######
#######################


high_c_prev_jm <- performance_plot_prev(data = performance_data,
                                        stat = "c_index", 
                                        ylimits = c(0.5,0.9), 
                                        ylabel= "Harrell's C-index (95% quantile range)", 
                                        dim = "3 longitudinal", 
                                        prob = "20% missing", 
                                        dgm_arg = "DGM = JM", 
                                        sample = "n = 1000",
                                        legend = "right")
high_c_prev_jm

low_c_prev_jm <- performance_plot_prev(data = performance_data, 
                                       stat = "c_index", 
                                       ylimits = c(0.4,0.9), 
                                       ylabel= "Harrell's C-index (95% quantile range)", 
                                       dim = "3 longitudinal", 
                                       prob = "80% missing", 
                                       dgm_arg = "DGM = JM",
                                       sample = "n = 200",
                                       legend = "right")
low_c_prev_jm

high_c_prev_tdcm <- performance_plot_prev(data = performance_data,
                                          stat = "c_index", 
                                          ylimits = c(0.5,0.9), 
                                          ylabel= "Harrell's C-index (95% quantile range)", 
                                          dim = "3 longitudinal", 
                                          prob = "20% missing", 
                                          dgm_arg = "DGM = TDCM", 
                                          sample = "n = 1000",
                                          legend = "right")
high_c_prev_tdcm

low_c_prev_tdcm <- performance_plot_prev(data = performance_data, 
                                         stat = "c_index", 
                                         ylimits = c(0.4,0.9), 
                                         ylabel= "Harrell's C-index (95% quantile range)", 
                                         dim = "3 longitudinal", 
                                         prob = "80% missing", 
                                         dgm_arg = "DGM = TDCM",
                                         sample = "n = 200",
                                         legend = "right")
low_c_prev_tdcm

#######################
####### Brier   #######
#######################

high_brier_prev_jm <- performance_plot_prev(data = performance_data,
                                            stat = "brier", 
                                            ylimits = c(0, 0.175), 
                                            ylabel= "Survival Brier score (95% quantile range)", 
                                            dim = "3 longitudinal", 
                                            prob = "20% missing", 
                                            dgm_arg = "DGM = JM",
                                            sample = "n = 1000",
                                            legend = "right")
high_brier_prev_jm

low_brier_prev_jm <- performance_plot_prev(data = performance_data,
                                           stat = "brier", 
                                           ylimits = c(0, 0.175), 
                                           ylabel= "Survival Brier score (95% quantile range)", 
                                           dim = "3 longitudinal", 
                                           prob = "80% missing", 
                                           dgm_arg = "DGM = JM",
                                           sample = "n = 200",
                                           legend = "right")
low_brier_prev_jm

high_brier_prev_tdcm <- performance_plot_prev(data = performance_data,
                                              stat = "brier", 
                                              ylimits = c(0, 0.175), 
                                              ylabel= "Survival Brier score (95% quantile range)", 
                                              dim = "3 longitudinal", 
                                              prob = "20% missing", 
                                              dgm_arg = "DGM = TDCM",
                                              sample = "n = 1000",
                                              legend = "right")
high_brier_prev_tdcm

low_brier_prev_tdcm <- performance_plot_prev(data = performance_data,
                                             stat = "brier", 
                                             ylimits = c(0, 0.175), 
                                             ylabel= "Survival Brier score (95% quantile range)", 
                                             dim = "3 longitudinal", 
                                             prob = "80% missing", 
                                             dgm_arg = "DGM = TDCM",
                                             sample = "n = 200",
                                             legend = "right")
low_brier_prev_tdcm


#######################
####### IPA     #######
#######################

high_ipa_prev_jm <- performance_plot_prev(data = brier_ipa,
                                          stat = "ipa", 
                                          ylimits = c(-25, 75), 
                                          ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                          dim = "3 longitudinal", 
                                          prob = "20% missing", 
                                          dgm_arg = "DGM = JM", 
                                          sample = "n = 1000",
                                          legend = "right")
high_ipa_prev_jm
low_ipa_prev_jm <- performance_plot_prev(data = brier_ipa,
                                         stat = "ipa", 
                                         ylimits = c(-25, 100), 
                                         ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                         dim = "3 longitudinal", 
                                         prob = "80% missing", 
                                         dgm_arg = "DGM = JM",  
                                         sample = "n = 200",
                                         legend = "right")
low_ipa_prev_jm
high_ipa_prev_tdcm <- performance_plot_prev(data = brier_ipa,
                                            stat = "ipa", 
                                            ylimits = c(-25, 75), 
                                            ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                            dim = "3 longitudinal", 
                                            prob = "20% missing", 
                                            dgm_arg = "DGM = TDCM", 
                                            sample = "n = 1000",
                                            legend = "right")
high_ipa_prev_tdcm
low_ipa_prev_tdcm <- performance_plot_prev(data = brier_ipa,
                                           stat = "ipa", 
                                           ylimits = c(-25, 75), 
                                           ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                           dim = "3 longitudinal", 
                                           prob = "80% missing", 
                                           dgm_arg = "DGM = TDCM",  
                                           sample = "n = 200",
                                           legend = "right")
low_ipa_prev_tdcm


##########################
### Optimism - C-index ###
##########################

high_optC_prev_jm <- performance_plot_prev(data = performance_data,
                                           stat = "c_index_adjusted", 
                                           ylimits = c(-0.2, 0.2), 
                                           ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                           dim = "3 longitudinal", 
                                           prob = "20% missing", 
                                           dgm_arg = "DGM = JM", 
                                           sample = "n = 1000",
                                           legend = "right")
high_optC_prev_jm
low_optC_prev_jm <- performance_plot_prev(data = performance_data,
                                          stat = "c_index_adjusted", 
                                          ylimits = c(-0.25, 0.5), 
                                          ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                          dim = "3 longitudinal", 
                                          prob = "80% missing", 
                                          dgm_arg = "DGM = JM",  
                                          sample = "n = 200",
                                          legend = "right")
low_optC_prev_jm
high_optC_prev_tdcm <- performance_plot_prev(data = performance_data,
                                             stat = "c_index_adjusted", 
                                             ylimits = c(-0.2, 0.2), 
                                             ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                             dim = "3 longitudinal", 
                                             prob = "20% missing", 
                                             dgm_arg = "DGM = TDCM", 
                                             sample = "n = 1000",
                                             legend = "right")
high_optC_prev_tdcm
low_optC_prev_tdcm <- performance_plot_prev(data = performance_data,
                                            stat = "c_index_adjusted", 
                                            ylimits = c(-0.25, 0.5), 
                                            ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                            dim = "3 longitudinal", 
                                            prob = "80% missing", 
                                            dgm_arg = "DGM = TDCM",  
                                            sample = "n = 200",
                                            legend = "right")
low_optC_prev_tdcm 

##########################
### Optimism - Brier   ###
##########################

high_optB_prev_jm <- performance_plot_prev(data = performance_data,
                                           stat = "brier_adjusted", 
                                           ylimits = c(-0.07, 0.07), 
                                           ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                           dim = "3 longitudinal", 
                                           prob = "20% missing", 
                                           dgm_arg = "DGM = JM",  
                                           sample = "n = 1000",
                                           legend = "right")
high_optB_prev_jm
low_optB_prev_jm <- performance_plot_prev(data = performance_data,
                                          stat = "brier_adjusted", 
                                          ylimits = c(-0.07, 0.07), 
                                          ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                          dim = "3 longitudinal", 
                                          prob = "80% missing", 
                                          dgm_arg = "DGM = JM", 
                                          sample = "n = 200",
                                          legend = "right")

low_optB_prev_jm
high_optB_prev_tdcm <- performance_plot_prev(data = performance_data,
                                             stat = "brier_adjusted", 
                                             ylimits = c(-0.07, 0.07), 
                                             ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                             dim = "3 longitudinal", 
                                             prob = "20% missing", 
                                             dgm_arg = "DGM = TDCM",  
                                             sample = "n = 1000",
                                             legend = "right")
high_optB_prev_tdcm
low_optB_prev_tdcm <- performance_plot_prev(data = performance_data,
                                            stat = "brier_adjusted", 
                                            ylimits = c(-0.07, 0.07), 
                                            ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                            dim = "3 longitudinal", 
                                            prob = "80% missing", 
                                            dgm_arg = "DGM = TDCM", 
                                            sample = "n = 200",
                                            legend = "right")

low_optB_prev_tdcm


#####################
#### Stability   ####
#####################

high_stability_prev_jm <- stability_plot_prev(stability_data, 
                                              ylimits = c(0, 0.9), 
                                              ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                              dim = "3 longitudinal", 
                                              prob = "20% missing", 
                                              dgm_arg = "DGM = JM", 
                                              sample = "n = 1000", 
                                              legend = "right")
high_stability_prev_jm
low_stability_prev_jm <- stability_plot_prev(stability_data, 
                                             ylimits = c(0, 1), 
                                             ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                             dim = "3 longitudinal", 
                                             prob = "20% missing", 
                                             dgm_arg = "DGM = JM", 
                                             sample = "n = 200", 
                                             legend = "right")
low_stability_prev_jm
high_stability_prev_tdcm <- stability_plot_prev(stability_data, 
                                                ylimits = c(0, 0.9), 
                                                ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                dim = "3 longitudinal", 
                                                prob = "20% missing", 
                                                dgm_arg = "DGM = TDCM", 
                                                sample = "n = 1000", 
                                                legend = "right")
high_stability_prev_tdcm 
low_stability_prev_tdcm <- stability_plot_prev(stability_data, 
                                               ylimits = c(0, 0.9), 
                                               ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                               dim = "3 longitudinal", 
                                               prob = "80% missing", 
                                               dgm_arg = "DGM = TDCM", 
                                               sample = "n = 200", 
                                               legend = "right")
low_stability_prev_tdcm



###################
### Consistency ###
###################

high_consistencyC_prev_jm <- consistency_plot_prev(data = performance_data,
                                                   stat = "c_index",
                                                   dim = "3 longitudinal",
                                                   prob = "20% missing",
                                                   dgm_arg = "DGM = JM",
                                                   sample = "n = 1000",
                                                   ylabel = "Harrell's C-index range\n across landmark times",
                                                   ylimits = c(0, 0.1),
                                                   legend = "right")
high_consistencyC_prev_jm
low_consistencyC_prev_jm <- consistency_plot_prev(data = performance_data,
                                                  stat = "c_index",
                                                  dim = "3 longitudinal",
                                                  prob = "80% missing",
                                                  dgm_arg = "DGM = JM",
                                                  sample = "n = 200",
                                                  ylabel = "Harrell's C-index range\n across landmark times",
                                                  ylimits = c(0, 0.07),
                                                  legend = "right")
low_consistencyC_prev_jm
high_consistencyC_prev_tdcm <- consistency_plot_prev(data = performance_data,
                                                     stat = "c_index",
                                                     dim = "3 longitudinal",
                                                     prob = "20% missing",
                                                     dgm_arg = "DGM = TDCM",
                                                     sample = "n = 1000",
                                                     ylabel = "Harrell's C-index range\n across landmark times",
                                                     ylimits = c(0, 0.07),
                                                     legend = "right")
high_consistencyC_prev_tdcm
low_consistencyC_prev_tdcm <- consistency_plot_prev(data = performance_data,
                                                    stat = "c_index",
                                                    dim = "3 longitudinal",
                                                    prob = "80% missing",
                                                    dgm_arg = "DGM = TDCM",
                                                    sample = "n = 200",
                                                    ylabel = "Harrell's C-index range\n across landmark times",
                                                    ylimits = c(0, 0.07),
                                                    legend = "right")
low_consistencyC_prev_tdcm



high_consistencyB_prev_jm <- consistency_plot_prev(data = performance_data,
                                                   stat = "brier",
                                                   dim = "3 longitudinal",
                                                   prob = "20% missing",
                                                   dgm_arg = "DGM = JM",
                                                   sample = "n = 1000",
                                                   ylabel = "Brier score range\n across landmark times",
                                                   ylimits = c(0, 0.1),
                                                   legend = "right")
high_consistencyB_prev_jm
low_consistencyB_prev_jm <- consistency_plot_prev(data = performance_data,
                                                  stat = "brier",
                                                  dim = "3 longitudinal",
                                                  prob = "80% missing",
                                                  dgm_arg = "DGM = JM",
                                                  sample = "n = 200",
                                                  ylabel = "Brier score range\n across landmark times",
                                                  ylimits = c(0, 0.1),
                                                  legend = "right")
low_consistencyB_prev_jm
high_consistencyB_prev_tdcm <- consistency_plot_prev(data = performance_data,
                                                     stat = "brier",
                                                     dim = "1 longitudinal",
                                                     prob = "20% missing",
                                                     dgm_arg = "DGM = TDCM",
                                                     sample = "n = 1000",
                                                     ylabel = "Brier score range\n across landmark times",
                                                     ylimits = c(0, 0.1),
                                                     legend = "right")
high_consistencyB_prev_tdcm
low_consistencyB_prev_tdcm <- consistency_plot_prev(data = performance_data,
                                                    stat = "brier",
                                                    dim = "3 longitudinal",
                                                    prob = "80% missing",
                                                    dgm_arg = "DGM = TDCM",
                                                    sample = "n = 200",
                                                    ylabel = "Brier score range\n across landmark times",
                                                    ylimits = c(0, 0.1),
                                                    legend = "right")
low_consistencyB_prev_tdcm

#####################
### Fitting times ###
#####################

high_fitting_prev_jm <- fitting_time_prev(data = performance_data,
                                          stat = "fit_time", 
                                          ylimits = c(0,10), 
                                          ylabel= "Model fitting time (minutes)", 
                                          dim = "3 longitudinal", 
                                          prob = "20% missing", 
                                          dgm_arg = "DGM = JM", 
                                          sample = "n = 1000",
                                          legend = "none")
high_fitting_prev_jm
low_fitting_prev_jm <- fitting_time_prev(data = performance_data,
                                         stat = "fit_time", 
                                         ylimits = c(0,5), 
                                         ylabel= "Model fitting time (minutes)", 
                                         dim = "3 longitudinal", 
                                         prob = "80% missing",
                                         dgm_arg = "DGM = JM",
                                         sample = "n = 200",
                                         legend = "right")
low_fitting_prev_jm
high_fitting_prev_tdcm <- fitting_time_prev(data = performance_data,
                                            stat = "fit_time", 
                                            ylimits = c(0,5), 
                                            ylabel= "Model fitting time (minutes)", 
                                            dim = "3 longitudinal", 
                                            prob = "20% missing", 
                                            dgm_arg = "DGM = TDCM", 
                                            sample = "n = 1000",
                                            legend = "none")
high_fitting_prev_tdcm
low_fitting_prev_tdcm <- fitting_time_prev(data = performance_data,
                                           stat = "fit_time", 
                                           ylimits = c(0,5), 
                                           ylabel= "Model fitting time (minutes)", 
                                           dim = "3 longitudinal", 
                                           prob = "80% missing",
                                           dgm_arg = "DGM = TDCM",
                                           sample = "n = 200",
                                           legend = "right")
low_fitting_prev_tdcm


########################
### Prediction times ###
########################

high_pred_time_prev_jm <- performance_plot_prev(data = performance_data,
                                                stat = "pred_time", 
                                                ylimits = c(0,50), 
                                                ylabel= "Prediction time (minutes)", 
                                                dim = "3 longitudinal", 
                                                prob = "20% missing", 
                                                dgm_arg = "DGM = JM", 
                                                sample = "n = 1000",
                                                legend = "right")
high_pred_time_prev_jm
low_pred_time_prev_jm <- performance_plot_prev(data = performance_data,
                                               stat = "pred_time", 
                                               ylimits = c(0,50), 
                                               ylabel= "Prediction time (minutes)", 
                                               dim = "3 longitudinal", 
                                               prob = "80% missing", 
                                               dgm_arg = "DGM = JM", 
                                               sample = "n = 200",
                                               legend = "right")
low_pred_time_prev_jm 
high_pred_time_prev_tdcm <- performance_plot_prev(data = performance_data,
                                                  stat = "pred_time", 
                                                  ylimits = c(0,50), 
                                                  ylabel= "Prediction time (minutes)", 
                                                  dim = "3 longitudinal", 
                                                  prob = "20% missing", 
                                                  dgm_arg = "DGM = TDCM", 
                                                  sample = "n = 1000",
                                                  legend = "right")
high_pred_time_prev_tdcm
low_pred_time_prev_tdcm <- performance_plot_prev(data = performance_data,
                                                 stat = "pred_time", 
                                                 ylimits = c(0,50), 
                                                 ylabel= "Prediction time (minutes)", 
                                                 dim = "3 longitudinal", 
                                                 prob = "80% missing", 
                                                 dgm_arg = "DGM = TDCM", 
                                                 sample = "n = 200",
                                                 legend = "right")
low_pred_time_prev_tdcm 


#######################################
### Follow-up missingness influence ###
#######################################
#*************************************#

#######################
####### C-index #######
#######################


high_c_prob_jm <- performance_plot_prob(data = performance_data,
                                        stat = "c_index", 
                                        ylimits = c(0.5,0.9), 
                                        ylabel= "Harrell's C-index (95% quantile range)", 
                                        dim = "3 longitudinal", 
                                        prev = "prev ~ 40%", 
                                        dgm_arg = "DGM = JM", 
                                        sample = "n = 1000",
                                        legend = "right")
high_c_prob_jm

low_c_prob_jm <- performance_plot_prob(data = performance_data, 
                                       stat = "c_index", 
                                       ylimits = c(0.4,0.9), 
                                       ylabel= "Harrell's C-index (95% quantile range)", 
                                       dim = "3 longitudinal", 
                                       prev = "prev ~ 10%",
                                       dgm_arg = "DGM = JM",
                                       sample = "n = 200",
                                       legend = "right")
low_c_prob_jm

high_c_prob_tdcm <- performance_plot_prob(data = performance_data,
                                          stat = "c_index", 
                                          ylimits = c(0.5,0.9), 
                                          ylabel= "Harrell's C-index (95% quantile range)", 
                                          dim = "3 longitudinal", 
                                          prev = "prev ~ 40%",
                                          dgm_arg = "DGM = TDCM", 
                                          sample = "n = 1000",
                                          legend = "right")
high_c_prob_tdcm

low_c_prob_tdcm <- performance_plot_prob(data = performance_data, 
                                         stat = "c_index", 
                                         ylimits = c(0.4,0.9), 
                                         ylabel= "Harrell's C-index (95% quantile range)", 
                                         dim = "3 longitudinal", 
                                         prev = "prev ~ 10%",
                                         dgm_arg = "DGM = TDCM",
                                         sample = "n = 200",
                                         legend = "right")
low_c_prob_tdcm

#######################
####### Brier   #######
#######################

high_brier_prob_jm <- performance_plot_prob(data = performance_data,
                                            stat = "brier", 
                                            ylimits = c(0, 0.175), 
                                            ylabel= "Survival Brier score (95% quantile range)", 
                                            dim = "3 longitudinal", 
                                            prev = "prev ~ 40%", 
                                            dgm_arg = "DGM = JM",
                                            sample = "n = 1000",
                                            legend = "right")
high_brier_prob_jm

low_brier_prob_jm <- performance_plot_prob(data = performance_data,
                                           stat = "brier", 
                                           ylimits = c(0, 0.05), 
                                           ylabel= "Survival Brier score (95% quantile range)", 
                                           dim = "3 longitudinal", 
                                           prev = "prev ~ 10%", 
                                           dgm_arg = "DGM = JM",
                                           sample = "n = 200",
                                           legend = "right")
low_brier_prob_jm

high_brier_prob_tdcm <- performance_plot_prob(data = performance_data,
                                              stat = "brier", 
                                              ylimits = c(0, 0.175), 
                                              ylabel= "Survival Brier score (95% quantile range)", 
                                              dim = "3 longitudinal", 
                                              prev = "prev ~ 40%", 
                                              dgm_arg = "DGM = TDCM",
                                              sample = "n = 1000",
                                              legend = "right")
high_brier_prob_tdcm

low_brier_prob_tdcm <- performance_plot_prob(data = performance_data,
                                             stat = "brier", 
                                             ylimits = c(0, 0.05), 
                                             ylabel= "Survival Brier score (95% quantile range)", 
                                             dim = "3 longitudinal", 
                                             prev = "prev ~ 10%",
                                             dgm_arg = "DGM = TDCM",
                                             sample = "n = 200",
                                             legend = "right")
low_brier_prob_tdcm


#######################
####### IPA     #######
#######################

high_ipa_prob_jm <- performance_plot_prob(data = brier_ipa,
                                          stat = "ipa", 
                                          ylimits = c(-25, 100), 
                                          ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                          dim = "3 longitudinal", 
                                          prev = "prev ~ 40%", 
                                          dgm_arg = "DGM = JM", 
                                          sample = "n = 1000",
                                          legend = "right")
high_ipa_prob_jm
low_ipa_prob_jm <- performance_plot_prob(data = brier_ipa,
                                         stat = "ipa", 
                                         ylimits = c(-25, 100), 
                                         ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                         dim = "3 longitudinal", 
                                         prev = "prev ~ 10%", 
                                         dgm_arg = "DGM = JM",  
                                         sample = "n = 200",
                                         legend = "right")
low_ipa_prob_jm
high_ipa_prob_tdcm <- performance_plot_prob(data = brier_ipa,
                                            stat = "ipa", 
                                            ylimits = c(-25, 75), 
                                            ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                            dim = "3 longitudinal", 
                                            prev = "prev ~ 40%",
                                            dgm_arg = "DGM = TDCM", 
                                            sample = "n = 1000",
                                            legend = "right")
high_ipa_prob_tdcm
low_ipa_prob_tdcm <- performance_plot_prob(data = brier_ipa,
                                           stat = "ipa", 
                                           ylimits = c(-30, 75), 
                                           ylabel= "Index of Predictive Accuracy\n(%, 95% quantile range)", 
                                           dim = "3 longitudinal", 
                                           prev = "prev ~ 10%", 
                                           dgm_arg = "DGM = TDCM",  
                                           sample = "n = 200",
                                           legend = "right")
low_ipa_prob_tdcm


##########################
### Optimism - C-index ###
##########################

high_optC_prob_jm <- performance_plot_prob(data = performance_data,
                                           stat = "c_index_adjusted", 
                                           ylimits = c(-0.2, 0.2), 
                                           ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                           dim = "3 longitudinal", 
                                           prev = "prev ~ 40%", 
                                           dgm_arg = "DGM = JM", 
                                           sample = "n = 1000",
                                           legend = "right")
high_optC_prob_jm
low_optC_prob_jm <- performance_plot_prob(data = performance_data,
                                          stat = "c_index_adjusted", 
                                          ylimits = c(-0.25, 0.5), 
                                          ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                          dim = "3 longitudinal", 
                                          prev = "prev ~ 10%", 
                                          dgm_arg = "DGM = JM",  
                                          sample = "n = 200",
                                          legend = "right")
low_optC_prob_jm
high_optC_prob_tdcm <- performance_plot_prob(data = performance_data,
                                             stat = "c_index_adjusted", 
                                             ylimits = c(-0.2, 0.2), 
                                             ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                             dim = "3 longitudinal", 
                                             prev = "prev ~ 40%", 
                                             dgm_arg = "DGM = TDCM", 
                                             sample = "n = 1000",
                                             legend = "right")
high_optC_prob_tdcm
low_optC_prob_tdcm <- performance_plot_prob(data = performance_data,
                                            stat = "c_index_adjusted", 
                                            ylimits = c(-0.25, 0.5), 
                                            ylabel= "Optimism estimate for C-index\n (95% quantile range)", 
                                            dim = "3 longitudinal", 
                                            prev = "prev ~ 10%", 
                                            dgm_arg = "DGM = TDCM",  
                                            sample = "n = 200",
                                            legend = "right")
low_optC_prob_tdcm 

##########################
### Optimism - Brier   ###
##########################

high_optB_prob_jm <- performance_plot_prob(data = performance_data,
                                           stat = "brier_adjusted", 
                                           ylimits = c(-0.07, 0.07), 
                                           ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                           dim = "3 longitudinal", 
                                           prev = "prev ~ 40%", 
                                           dgm_arg = "DGM = JM",  
                                           sample = "n = 1000",
                                           legend = "right")
high_optB_prob_jm
low_optB_prob_jm <- performance_plot_prob(data = performance_data,
                                          stat = "brier_adjusted", 
                                          ylimits = c(-0.07, 0.07), 
                                          ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                          dim = "3 longitudinal", 
                                          prev = "prev ~ 10%", 
                                          dgm_arg = "DGM = JM", 
                                          sample = "n = 200",
                                          legend = "right")

low_optB_prob_jm
high_optB_prob_tdcm <- performance_plot_prob(data = performance_data,
                                             stat = "brier_adjusted", 
                                             ylimits = c(-0.07, 0.07), 
                                             ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                             dim = "3 longitudinal", 
                                             prev = "prev ~ 40%",
                                             dgm_arg = "DGM = TDCM",  
                                             sample = "n = 1000",
                                             legend = "right")
high_optB_prob_tdcm
low_optB_prob_tdcm <- performance_plot_prob(data = performance_data,
                                            stat = "brier_adjusted", 
                                            ylimits = c(-0.07, 0.07), 
                                            ylabel= "Optimism estimate for Brier score\n (95% quantile range)", 
                                            dim = "3 longitudinal", 
                                            prev = "prev ~ 10%", 
                                            dgm_arg = "DGM = TDCM", 
                                            sample = "n = 200",
                                            legend = "right")

low_optB_prob_tdcm


#####################
#### Stability   ####
#####################

high_stability_prob_jm <- stability_plot_prob(stability_data, 
                                              ylimits = c(0, 0.9), 
                                              ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                              dim = "3 longitudinal", 
                                              prev = "prev ~ 40%", 
                                              dgm_arg = "DGM = JM", 
                                              sample = "n = 1000", 
                                              legend = "right")
high_stability_prob_jm
low_stability_prob_jm <- stability_plot_prob(stability_data, 
                                             ylimits = c(0, 1), 
                                             ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                             dim = "3 longitudinal", 
                                             prev = "prev ~ 40%", 
                                             dgm_arg = "DGM = JM", 
                                             sample = "n = 200", 
                                             legend = "right")
low_stability_prob_jm
high_stability_prob_tdcm <- stability_plot_prob(stability_data, 
                                                ylimits = c(0, 0.9), 
                                                ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                                dim = "3 longitudinal", 
                                                prev = "prev ~ 40%", 
                                                dgm_arg = "DGM = TDCM", 
                                                sample = "n = 1000", 
                                                legend = "right")
high_stability_prob_tdcm 
low_stability_prob_tdcm <- stability_plot_prob(stability_data, 
                                               ylimits = c(0, 0.9), 
                                               ylabel = "Prediction range across simulated models\n(95% quantile range across subjects)", 
                                               dim = "3 longitudinal", 
                                               prev = "prev ~ 10%", 
                                               dgm_arg = "DGM = TDCM", 
                                               sample = "n = 200", 
                                               legend = "right")
low_stability_prob_tdcm



###################
### Consistency ###
###################

high_consistencyC_prob_jm <- consistency_plot_prob(data = performance_data,
                                                   stat = "c_index",
                                                   dim = "3 longitudinal",
                                                   prev = "prev ~ 40%",
                                                   dgm_arg = "DGM = JM",
                                                   sample = "n = 1000",
                                                   ylabel = "Harrell's C-index range\n across landmark times",
                                                   ylimits = c(0, 0.1),
                                                   legend = "right")
high_consistencyC_prob_jm
low_consistencyC_prob_jm <- consistency_plot_prob(data = performance_data,
                                                  stat = "c_index",
                                                  dim = "3 longitudinal",
                                                  prev = "prev ~ 10%",
                                                  dgm_arg = "DGM = JM",
                                                  sample = "n = 200",
                                                  ylabel = "Harrell's C-index range\n across landmark times",
                                                  ylimits = c(0, 0.07),
                                                  legend = "right")
low_consistencyC_prob_jm
high_consistencyC_prob_tdcm <- consistency_plot_prob(data = performance_data,
                                                     stat = "c_index",
                                                     dim = "3 longitudinal",
                                                     prev = "prev ~ 40%",
                                                     dgm_arg = "DGM = TDCM",
                                                     sample = "n = 1000",
                                                     ylabel = "Harrell's C-index range\n across landmark times",
                                                     ylimits = c(0, 0.07),
                                                     legend = "right")
high_consistencyC_prob_tdcm
low_consistencyC_prob_tdcm <- consistency_plot_prob(data = performance_data,
                                                    stat = "c_index",
                                                    dim = "3 longitudinal",
                                                    prev = "prev ~ 10%",
                                                    dgm_arg = "DGM = TDCM",
                                                    sample = "n = 200",
                                                    ylabel = "Harrell's C-index range\n across landmark times",
                                                    ylimits = c(0, 0.07),
                                                    legend = "right")
low_consistencyC_prob_tdcm



high_consistencyB_prob_jm <- consistency_plot_prob(data = performance_data,
                                                   stat = "brier",
                                                   dim = "3 longitudinal",
                                                   prev = "prev ~ 40%",
                                                   dgm_arg = "DGM = JM",
                                                   sample = "n = 1000",
                                                   ylabel = "Brier score range\n across landmark times",
                                                   ylimits = c(0, 0.1),
                                                   legend = "right")
high_consistencyB_prob_jm
low_consistencyB_prob_jm <- consistency_plot_prob(data = performance_data,
                                                  stat = "brier",
                                                  dim = "3 longitudinal",
                                                  prev = "prev ~ 10%",
                                                  dgm_arg = "DGM = JM",
                                                  sample = "n = 200",
                                                  ylabel = "Brier score range\n across landmark times",
                                                  ylimits = c(0, 0.1),
                                                  legend = "right")
low_consistencyB_prob_jm
high_consistencyB_prob_tdcm <- consistency_plot_prob(data = performance_data,
                                                     stat = "brier",
                                                     dim = "3 longitudinal",
                                                     prev = "prev ~ 40%",
                                                     dgm_arg = "DGM = TDCM",
                                                     sample = "n = 1000",
                                                     ylabel = "Brier score range\n across landmark times",
                                                     ylimits = c(0, 0.1),
                                                     legend = "right")
high_consistencyB_prob_tdcm
low_consistencyB_prob_tdcm <- consistency_plot_prob(data = performance_data,
                                                    stat = "brier",
                                                    dim = "3 longitudinal",
                                                    prev = "prev ~ 10%",
                                                    dgm_arg = "DGM = TDCM",
                                                    sample = "n = 200",
                                                    ylabel = "Brier score range\n across landmark times",
                                                    ylimits = c(0, 0.1),
                                                    legend = "right")
low_consistencyB_prob_tdcm

#####################
### Fitting times ###
#####################

high_fitting_prob_jm <- fitting_time_prob(data = performance_data,
                                          stat = "fit_time", 
                                          ylimits = c(0,10), 
                                          ylabel= "Model fitting time (minutes)", 
                                          dim = "3 longitudinal", 
                                          prev = "prev ~ 40%", 
                                          dgm_arg = "DGM = JM", 
                                          sample = "n = 1000",
                                          legend = "right")
high_fitting_prob_jm
low_fitting_prob_jm <- fitting_time_prob(data = performance_data,
                                         stat = "fit_time", 
                                         ylimits = c(0,5), 
                                         ylabel= "Model fitting time (minutes)", 
                                         dim = "3 longitudinal", 
                                         prev = "prev ~ 10%",
                                         dgm_arg = "DGM = JM",
                                         sample = "n = 200",
                                         legend = "right")
low_fitting_prob_jm
high_fitting_prob_tdcm <- fitting_time_prob(data = performance_data,
                                            stat = "fit_time", 
                                            ylimits = c(0,5), 
                                            ylabel= "Model fitting time (minutes)", 
                                            dim = "3 longitudinal", 
                                            prev = "prev ~ 40%", 
                                            dgm_arg = "DGM = TDCM", 
                                            sample = "n = 1000",
                                            legend = "right")
high_fitting_prob_tdcm
low_fitting_prob_tdcm <- fitting_time_prob(data = performance_data,
                                           stat = "fit_time", 
                                           ylimits = c(0,5), 
                                           ylabel= "Model fitting time (minutes)", 
                                           dim = "3 longitudinal", 
                                           prev = "prev ~ 10%",
                                           dgm_arg = "DGM = TDCM",
                                           sample = "n = 200",
                                           legend = "right")
low_fitting_prob_tdcm


########################
### Prediction times ###
########################

high_pred_time_prob_jm <- performance_plot_prob(data = performance_data,
                                                stat = "pred_time", 
                                                ylimits = c(0,50), 
                                                ylabel= "Prediction time (minutes)", 
                                                dim = "3 longitudinal", 
                                                prev = "prev ~ 40%", 
                                                dgm_arg = "DGM = JM", 
                                                sample = "n = 1000",
                                                legend = "right")
high_pred_time_prob_jm
low_pred_time_prob_jm <- performance_plot_prob(data = performance_data,
                                               stat = "pred_time", 
                                               ylimits = c(0,50), 
                                               ylabel= "Prediction time (minutes)", 
                                               dim = "3 longitudinal", 
                                               prev = "prev ~ 10%", 
                                               dgm_arg = "DGM = JM", 
                                               sample = "n = 200",
                                               legend = "right")
low_pred_time_prob_jm 
high_pred_time_prob_tdcm <- performance_plot_prob(data = performance_data,
                                                  stat = "pred_time", 
                                                  ylimits = c(0,50), 
                                                  ylabel= "Prediction time (minutes)", 
                                                  dim = "3 longitudinal", 
                                                  prev = "prev ~ 40%", 
                                                  dgm_arg = "DGM = TDCM", 
                                                  sample = "n = 1000",
                                                  legend = "right")
high_pred_time_prob_tdcm
low_pred_time_prob_tdcm <- performance_plot_prob(data = performance_data,
                                                 stat = "pred_time", 
                                                 ylimits = c(0,50), 
                                                 ylabel= "Prediction time (minutes)", 
                                                 dim = "3 longitudinal", 
                                                 prev = "prev ~ 10%", 
                                                 dgm_arg = "DGM = TDCM", 
                                                 sample = "n = 200",
                                                 legend = "right")
low_pred_time_prob_tdcm 