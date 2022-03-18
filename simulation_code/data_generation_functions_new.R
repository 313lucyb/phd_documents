library(splus2R)
library(JMbayes2)
library(dplyr)
library(survival)
library(truncnorm)
##### DGF function : Joint model, 1 longitudinal and 5 fixed

jm_simulate_one<-function(seed, n, K, t_max, betas, random_cov_matrix, 
                          noise, upp_Cens, shape_wb, alpha, gammas, probability_missing,
                          prop_female){
  
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max
  set.seed(seed)
  time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max)))))
  print(length(time))
  DF <- data.frame(id = rep(seq_len(n), each = K),
                   time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                   sex = rep(rbinom(n, 1, prop_female), each = K))
  
  DF$age <- rep(rtruncnorm(n, a=0, b=Inf, mean = 56.2, sd = 12.3), each = K)
  DF$disease_duration <- round(rep(rtruncnorm(n, a=0, b=Inf, mean = 11.96, sd = 9.7), each = K ), 0)
  DF$bmi <- rep(rtruncnorm(n, a=0, b=Inf, mean = 27.4, sd = 6.6), each = K)
  DF$patient_disease <- round(rep(rtruncnorm(n, a= 0, b = 3, 2, 0.69), each = K), 0)
  
  # design matrices for the fixed and random effects
  X <- model.matrix(~ sex * time, data = DF)
  Z <- model.matrix(~ time, data = DF)
  
  
  # we simulate random effects
  b <- rmvnorm(n, mean = c(0,0), cov = random_cov_matrix)
  # linear predictor
  eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
  # longitudinal predictors assumed to be normally distributed
  mu_y = eta_y
  DF$clinician_disease <- mu_y + rnorm(n*K, 0, noise)
  
  
  W <- model.matrix(~ sex + age + disease_duration + bmi + patient_disease, data = DF[!duplicated(DF$id), ])
  # linear predictor for the survival model
  eta_t <- as.vector(W %*% gammas)
  # to simulate event times we use inverse transform sampling
  # (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
  # to find t, such that S(t) = u, where S(.) is the survival function, and u a 
  # number from the Unif(0, 1) distribution. The function below calculates 
  # log(u) - log(S(t)), and for a given u, we want to find t for which it equals
  # zero. We do that below using the uniroot() function
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f <- as.vector(X_at_s %*% betas +
                       rowSums(Z_at_s * b[rep(i, nrow(Z_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f * alpha)
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  # we simulate the event times
  u <- runif(n)
  trueTimes <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  # we use fixed Type I right censoring denoting the end of the trial.
  Ctimes <- upp_Cens
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  
  # we keep the longitudinal measurements before the event times
  DF$Time <- Time[DF$id]
  DF$event <- event[DF$id]
  DF <- DF[DF$time <= DF$Time, ]
  
  # Missing indicators generated to monitor sparsity
  DF$missing_1 <- rbinom(length(DF$time), 1, probability_missing)
  DF <- DF %>% mutate(clinician_disease = ifelse((missing_1 == 1 & time != 0), NA, clinician_disease))
  DF <- DF %>% dplyr::select(-missing_1)
  return(DF)
}


##### DGF function: Joint model, 2 longitudinal and 4 fixed

jm_simulate_two<-function(seed, n, K, t_max, betas, betas2, random_cov_matrix, 
                          noise1, noise2, upp_Cens, shape_wb, alpha, alpha2, gammas, 
                          probability_missing, prop_female){
  set.seed(seed)
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max
  DF <- data.frame(id = rep(seq_len(n), each = K),
                   time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                   sex = rep(rbinom(n, 1, prop_female), each = K))
  
  DF$age <- rep(rtruncnorm(n, a=0, b=Inf, mean = 56.2, sd = 12.3), each = K)
  DF$disease_duration <- round(rep(rtruncnorm(n, a = 0 , b = Inf , mean = 11.96, sd = 9.7), each = K ), 0)
  DF$bmi <- rep(rtruncnorm(n, a=0, b=Inf, mean = 27.4, sd = 6.6), each = K)
  # design matrices for the fixed and random effects
  X <- model.matrix(~ sex * time, data = DF)
  Z <- model.matrix(~ time, data = DF)
  
  
  # we simulate random effects for all longitudinal outcomes
  b <- rmvnorm(n, mean = c(0,0,0,0), cov = random_cov_matrix)
  b1 = b[,(1:2)]
  b2 = b[,(3:4)]
  
  # linear predictor (longitudinal outcome 1 - DAS28)
  eta_y1 <- as.vector(X %*% betas + rowSums(Z * b1[DF$id, ]))
  # mean of the Normal distribution
  mu_y1 = eta_y1
  DF$clinician_disease <- mu_y1 + rnorm(n*K, 0, noise1)
  
  
  # linear predictor (longitudinal outcome 2 - HAQ)
  eta_y2 <- as.vector(X %*% betas2 + rowSums(Z * b2[DF$id, ]))
  # mean of the Normal distribution
  mu_y2 = eta_y2
  DF$patient_disease <- mu_y2 + rnorm(n*K, 0, noise2)
  
  W <- model.matrix(~ sex + age + disease_duration + bmi, data = DF[!duplicated(DF$id), ])
  # linear predictor for the survival model
  eta_t <- as.vector(W %*% gammas)
  # to simulate event times we use inverse transform sampling
  # (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
  # to find t, such that S(t) = u, where S(.) is the survival function, and u a 
  # number from the Unif(0, 1) distribution. The function below calculates 
  # log(u) - log(S(t)), and for a given u, we want to find t for which it equals
  # zero. We do that below using the uniroot() function
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f <- as.vector(X_at_s %*% betas +
                       rowSums(Z_at_s * b1[rep(i, nrow(Z_at_s)), ]))
      
      # the linear predictor from the mixed model evaluated at time s
      f2 <- as.vector(X_at_s %*% betas2 +
                        rowSums(Z_at_s * b2[rep(i, nrow(Z_at_s)), ]))
      
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f*alpha + f2*alpha2)
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  # we simulate the event times
  u <- runif(n)
  trueTimes <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  # we use fixed Type I right censoring denoting the end of the trial.
  Ctimes <- upp_Cens
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  
  # we keep the longitudinal measurements before the event times
  DF$Time <- Time[DF$id]
  DF$event <- event[DF$id]
  DF <- DF[DF$time <= DF$Time, ]
  
  # Missing indicators generated to monitor sparsity of longitudinal outcomes
  DF$missing_1 <- rbinom(length(DF$time), 1, probability_missing)
  DF$missing_2 <- rbinom(length(DF$time), 1, probability_missing)
  DF <- DF %>% mutate(clinician_disease = ifelse(missing_1 == 1 & time != 0, NA, clinician_disease),
                      patient_disease = ifelse(missing_2 == 1 & time != 0, NA, patient_disease))
  DF <- DF %>% dplyr::select(-missing_1)
  DF <- DF %>% dplyr::select(-missing_2)
  return(DF)
}


##### DGF function: Joint model, 3 longitudinal and 2 fixed

jm_simulate_three<-function(seed, n, K, t_max, betas, betas2, betas3, random_cov_matrix, 
                          noise1, noise2, noise3, upp_Cens, shape_wb, alpha, alpha2, alpha3, gammas,
                          probability_missing, prop_female){
  
  set.seed(seed)
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max
  DF <- data.frame(id = rep(seq_len(n), each = K),
                   time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                   sex = rep(rbinom(n, 1, prop_female), each = K))
  
  DF$age <- rep(rtruncnorm(n, a=0, b=Inf, mean = 56.2, sd = 12.3), each = K)
  DF$disease_duration <- round(rep(rtruncnorm(n, a=0, b = Inf, mean = 11.96, sd = 9.7), each = K ), 0)
  
  # design matrices for the fixed and random effects
  X <- model.matrix(~ sex * time, data = DF)
  Z <- model.matrix(~ time, data = DF)
  
  
  # we simulate random effects
  b <- rmvnorm(n, mean = c(0,0,0,0,0,0), cov = random_cov_matrix)
  b1 = b[,(1:2)]
  b2 = b[,(3:4)]
  b3 = b[,(5:6)]
  
  # linear predictor (DAS28)
  eta_y1 <- as.vector(X %*% betas + rowSums(Z * b1[DF$id, ]))
  mu_y1 = eta_y1
  DF$clinician_disease <- mu_y1 + rnorm(n*K, 0, noise1)
  
  
  # linear predictor (HAQ)
  eta_y2 <- as.vector(X %*% betas2 + rowSums(Z * b2[DF$id, ]))
  mu_y2 = eta_y2
  DF$patient_disease <- mu_y2 + rnorm(n*K, 0, noise2)
  
  # linear predictor (BMI)
  eta_y3 <- as.vector(X %*% betas3 + rowSums(Z * b3[DF$id, ]))
  mu_y3 = eta_y3
  DF$bmi <- mu_y3 + rnorm(n*K, 0, noise3)
  
  
  W <- model.matrix(~ sex + age + disease_duration, data = DF[!duplicated(DF$id), ])
  # linear predictor for the survival model
  eta_t <- as.vector(W %*% gammas)
  # to simulate event times we use inverse transform sampling
  # (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
  # to find t, such that S(t) = u, where S(.) is the survival function, and u a 
  # number from the Unif(0, 1) distribution. The function below calculates 
  # log(u) - log(S(t)), and for a given u, we want to find t for which it equals
  # zero. We do that below using the uniroot() function
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z_at_s <- cbind(1, s)
      # the linear predictor from mixed model 1 (DAS28) evaluated at time s
      f <- as.vector(X_at_s %*% betas +
                       rowSums(Z_at_s * b1[rep(i, nrow(Z_at_s)), ]))
      
      # the linear predictor from mixed model 2 (HAQ) evaluated at time s
      f2 <- as.vector(X_at_s %*% betas2 +
                        rowSums(Z_at_s * b2[rep(i, nrow(Z_at_s)), ]))
      
      # the linear predictor from mixed model 3 (BMI) evaluated at time s
      f3 <- as.vector(X_at_s %*% betas3 +
                        rowSums(Z_at_s * b3[rep(i, nrow(Z_at_s)), ]))
      
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f*alpha + f2*alpha2 + f3*alpha3)
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  # we simulate the event times
  u <- runif(n)
  trueTimes <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  # we use fixed Type I right censoring denoting the end of the trial.
  Ctimes <- upp_Cens
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  
  # we keep the longitudinal measurements before the event times
  DF$Time <- Time[DF$id]
  DF$event <- event[DF$id]
  DF <- DF[DF$time <= DF$Time, ]
  
  # Missing indicators generated to monitor sparsity of longitudinal outcomes
  DF$missing_1 <- rbinom(length(DF$time), 1, probability_missing)
  DF$missing_2 <- rbinom(length(DF$time), 1, probability_missing)
  DF$missing_3 <- rbinom(length(DF$time), 1, probability_missing)
  DF <- DF %>% mutate(clinician_disease = ifelse(missing_1 == 1 & time != 0, NA, clinician_disease),
                      patient_disease = ifelse(missing_2 == 1 & time != 0, NA, patient_disease),
                      bmi = ifelse(missing_3 == 1 & time != 0, NA, bmi))
  DF <- DF %>% dplyr::select(-missing_1)
  DF <- DF %>% dplyr::select(-missing_2)
  DF <- DF %>% dplyr::select(-missing_3)
  return(DF)
}

##### DGF function: TDCM

simulate_tdcm <- function(n, p, max_time, measurement_count, 
                          lambda, rho, n_fixed, n_vary, beta, 
                          fixed_mean, fixed_var, vary_mean,
                          vary_var, fixed_left, fixed_right, vary_left, vary_right, fixed_names, vary_names, re_mean, 
                          re_var, p_female, cens_time, seed){
  
  # n = sample size
  # p = probability of covariate change [single value across subjects]
  # max_time = latest measurement time
  # measurement_count = number of measurement times
  # lambda = scale parameter for Weibull hazard distribution
  # rho = shape parameter for Weibull hazard distribution
  # n_fixed = number of time-constant covariates (not including gender)
  # n_vary = number of time-varying covariates
  # beta = association parameter vector
  # fixed_means = vector of mean parameters for fixed covariates
  # fixed_vars = vector of variance parameters for fixed covariates
  # fixed/vary_left/right = Vector of upper and lower bounds for truncated Normal distributions
  # re_mean = mean of within-subject standard deviations for each of the time-varying covariates
  # re_var = standard deviation of within-subject standard deviations for each of the time-varying covariates
  # p_female = proportion female in the dataset
  # cens_time = end of study/censoring time
  
  set.seed(seed)
  id = 1:n
  data = data.frame()
  
  # Pop'n distribution of within-subject variability for longitudinal covariates
  
  re<-matrix(NA, nrow=n, ncol=3)
  re[,1] <- rnorm(n, re_mean[1], re_var[1])  
  re[,2] <- rnorm(n, re_mean[2], re_var[2])
  re[,3] <- rnorm(n, re_mean[3], re_var[3])
  
  # Subject specific matrices 
  
  for (i in 1:n){
    subj_covariate_changes = rbinom(measurement_count, 1, p)
    #times<-seq(0,max_time, length.out=measurement_count)
    times <- c(0, sort(runif(measurement_count - 1, 0, max_time)))
    subject_data = matrix(0, nrow=measurement_count, ncol = n_fixed + n_vary + 5)
    subject_data[,1]<-rep(id[i], times = measurement_count)
    subject_data[,2]= times
    subject_data[,3] <- rep(rbinom(1,1, p_female), times = measurement_count)
    for (g in 1:(n_fixed)){
      subject_data[, (3 + g)] <- rep(rtruncnorm(1, a = fixed_left[g], b = fixed_right[g], mean = fixed_mean[g], sd =fixed_var[g]), times = measurement_count)
    }
    for (k in 1:n_vary){
      subject_data[1, (3 + n_fixed + k)] = rtruncnorm(1, a = vary_left[k], b = vary_right[k], mean = vary_mean[k], sd = vary_var[k])
      
      # longitudinal predictors only updated if change indicator equal to one
      
      for (j in (2:measurement_count)){
        if (subj_covariate_changes[j] == 1){
          subject_data[j, (3 + n_fixed + k)] = subject_data[(j-1), (3 + n_fixed + k)] + rnorm(1, 0, sqrt(re[i,k]^2))
        } else {
          subject_data[j, (3 + n_fixed + k)] = subject_data[(j-1), (3 + n_fixed + k)]
        }
      }
    }
    subject_data[, (3 + n_fixed + n_vary + 1)] <- subj_covariate_changes
    u = runif(measurement_count)
    # Subject design matrix for all observation times
    X_mat = subject_data[,(3:(3 + n_fixed + n_vary))]
    # Survival times generated using covariates available at each measurement time (LOCF if not updated)
    for (t in (1:measurement_count)){
      subject_data[t, (5 + n_fixed + n_vary)] = (-log(u[t])/ (lambda*exp(X_mat[t,] %*% beta)))^(1 / rho)
    }
    
    subject_data <- as.data.frame(subject_data)
    names(subject_data) <- c("id", "start", "sex", fixed_names, vary_names, "change", "T")
    
    # Only updated survival times, due to covariate changes considered for comparison with next measurement time
    # Survival times added to measurement time, as baseline assumed zero above
    
    subject_data <- subject_data %>% 
      mutate(stop=lead(start), stop = ifelse(is.na(stop), T, stop),
             T = ifelse(change == 0 & start != 0, lag(T), T),
             eventT = ifelse(change == 1 | start == 0, start + T, NA)) 
    
    # Append all subject-specific matrices together
    
    data = rbind(data, subject_data)
  }        
  data<-as.data.frame(data)
  
  # Stop observation when survival time does not exceed subsequent observation time
  # Remove observations after the true survival time
  
  data<-data %>%
    dplyr::group_by(id) %>%
    tidyr::fill(eventT, .direction = "down") %>%
    mutate(cut = ifelse(stop > eventT, 1, 0),
           del = cumsum(cut)) %>%
    filter(del == 0 | cut == del) %>%
    mutate(event = cut, stop = ifelse(event == 1, eventT, stop)) %>%
    dplyr::select(id, start, stop, event, sex, all_of(fixed_names), all_of(vary_names))
  
  # Put data into same format as joint model 
  
  data<-as.data.frame(data)
  new_data <- data %>% mutate(timevent = ifelse(event == 1, stop, 0)) %>%
    group_by(id) %>% mutate(Event = max(event),  time = start) %>% 
    mutate(timevent = ifelse(Event == 0, max(time), timevent), EventTime = max(timevent)) %>% 
    dplyr::select(id, time, sex, all_of(fixed_names), all_of(vary_names), EventTime, Event)
  names(new_data)<-c("id", "time", "sex", all_of(fixed_names), all_of(vary_names), "Time", "event")
  
  # Impose censoring at cens_time
  new_data <- new_data %>% mutate(event = ifelse(Time <= cens_time & event == 1, 1, 0), 
                                  Time = ifelse(Time > cens_time, cens_time, Time)) %>% filter(time <= Time)
  return(new_data)
}

simulate_data <-function(dgm, dimension, probability_missing, n, j, weibull_shape, weibull_scale, i){
  if (dgm == "jm"){
    if (dimension == 1){
      data <- jm_simulate_one(
        seed = i,
        n=n,
        K=j,
        t_max = 10,
        betas=c(5.3614, 0.2177, -0.8825, 0.0789),
        random_cov_matrix = matrix(c(0.439, 0.0942, 0.0942, 0.027),2,2),
        upp_Cens=4,
        shape_wb = weibull_shape,
        alpha = 0.2315,
        gammas=c("intercept" = log(weibull_scale), "sex" = -0.3693, "age" = 0.0275, 
                 "disease_duration" = 0.0096,
                 "patient_disease" = 0.235, "bmi" = 0.0081),
        noise = 1.2668,
        probability_missing = probability_missing,
        prop_female = 0.75
      )
    }
    else if (dimension == 2){
      data <- jm_simulate_two(seed=i, 
                              n=n, 
                              K=j, 
                              t_max=10, 
                              betas=c(5.3543, 0.2265, -0.8751, 0.0766), 
                              betas2=c(1.7586, 0.1840, -0.1573, 0.0286),
                              random_cov_matrix = matrix(c(0.472, 0.0828, 0.271, 0.0219, 
                                                           0.0828, 0.0297, 0.0226, 0.0142, 
                                                           0.271, 0.0226, 0.383, 0.0134, 
                                                           0.0219, 0.0142, 0.0134, 0.0179),4,4), 
                              upp_Cens=4, 
                              shape_wb=0.75, 
                              alpha=0.1669,
                              alpha2 =0.3527,
                              gammas=c("intercept" = log(weibull_scale), "sex" = -0.3821, "age" = 0.0262, 
                                       "disease_duration" = 0.00788, "bmi" = 0.0065),
                              noise1 = 1.2613,
                              noise2 = 0.3401,
                              probability_missing = probability_missing,
                              prop_female = 0.75
      )
    }
    else if (dimension == 3){
      data<-jm_simulate_three(
        seed=i, 
        n=n, 
        K=j, 
        t_max=10, 
        betas=c(5.3541, 0.2235, -0.8757, 0.0839), 
        betas2=c(1.7582, 0.1848, -0.156, 0.0286),
        betas3=c(27.5809, -0.1630, 0.2459, -0.0567),
        random_cov_matrix = matrix(c(0.469, 0.0872, 0.267, 0.0249, 0.316, -0.00355,
                                     0.0872, 0.0282, 0.0261, 0.0129, 0.148, -0.0238,
                                     0.267, 0.0261, 0.383, 0.0135, 0.473, -0.00507,
                                     0.0249, 0.0129, 0.0135, 0.0182, 0.0775, -0.00372,
                                     0.316, 0.148, 0.473, 0.0775, 41.0175, -0.704,
                                     -0.00355, -0.0238, -0.00507, -0.00372, -0.704, 0.602),6,6), 
        upp_Cens=4, 
        shape_wb=weibull_shape, 
        alpha=0.1914,
        alpha2 =0.3364,
        alpha3 = 0.0071,
        gammas=c("intercept" = log(weibull_scale), "sex" = -0.3878, "age" = 0.0263, 
                 "disease_duration" = 0.0079),
        noise1 = 1.262,
        noise2 = 0.34,
        noise3 = 1.5867,
        probability_missing = probability_missing,
        prop_female = 0.75
      )
    }
  } else if (dgm == "tdcm"){
    if (dimension == 1){
      data <- simulate_tdcm(n=n, 
                            p= 1- probability_missing, 
                            max_time = 10,
                            cens_time = 4,
                            measurement_count = j, 
                            lambda = weibull_scale, 
                            rho = weibull_shape, 
                            n_fixed = 5 - dimension, 
                            n_vary = dimension, 
                            beta = c(-0.312, 0.027, 0.0102, 0.0078, 0.259, 0.0458), 
                            fixed_mean= c(56.2, 11.95, 27.38, 2),
                            fixed_var = c(12.27, 9.67, 6.6, 0.69),
                            vary_mean = c(6.4),
                            vary_var = c(1.05),
                            fixed_left = c(0, 0, 0, 0, 0),
                            fixed_right = c(Inf, Inf, Inf, Inf, 3),
                            vary_left = c(0),
                            vary_right = c(10),
                            fixed_names = c("age", "disease_duration", "bmi", "patient_disease"),
                            vary_names = c("clinician_disease"), 
                            re_mean = c(0.58,0.18,0.93),
                            re_var = c(1.07, 0.22, 0.7),
                            p_female = 0.75, 
                            seed = i)
  } else if (dimension == 2){
    data <- simulate_tdcm(n=n, 
                          p= 1- probability_missing, 
                          max_time = 10, 
                          cens_time = 4,
                          measurement_count = j, 
                          lambda = weibull_scale, 
                          rho = weibull_shape, 
                          n_fixed = 5 - dimension, 
                          n_vary = dimension, 
                          beta = c(-0.320, 0.026, 0.0097, 0.007, 0.277, 0.0301), 
                          fixed_mean= c(56.2, 11.95, 27.38),
                          fixed_var = c(12.27, 9.67, 6.6),
                          vary_mean = c(2, 6.4),
                          vary_var = c(0.69, 1.05),
                          fixed_left = c(0, 0, 0, 0),
                          fixed_right = c(Inf, Inf, Inf, Inf),
                          vary_left = c(0, 0),
                          vary_right = c(3, 10),
                          fixed_names = c("age", "disease_duration", "bmi"),
                          vary_names = c("patient_disease", "clinician_disease"), 
                          re_mean = c(0.58,0.18,0.93),
                          re_var = c(1.07, 0.22, 0.7),
                          p_female = 0.75, 
                          seed = i)
  } else if (dimension == 3){
    data <- simulate_tdcm(n=n, 
                          p= 1 - probability_missing, 
                          max_time = 10, 
                          cens_time = 4,
                          measurement_count = j, 
                          lambda = weibull_scale, 
                          rho = weibull_shape, 
                          n_fixed = 5 - dimension, 
                          n_vary = dimension, 
                          beta = c(-0.320, 0.026, 0.0097, 0.0077, 0.277, 0.0305), 
                          fixed_mean= c(56.2, 11.95),
                          fixed_var = c(12.27, 9.67),
                          vary_mean = c(27.38, 2, 6.4),
                          vary_var = c(6.6, 0.69, 1.05),
                          fixed_left = c(0, 0, 0),
                          fixed_right = c(Inf, Inf, Inf),
                          vary_left = c(0, 0, 0),
                          vary_right = c(Inf, 3, 10),
                          fixed_names = c("age", "disease_duration"),
                          vary_names = c("bmi", "patient_disease", "clinician_disease"), 
                          re_mean = c(0.58,0.18,0.93),
                          re_var = c(1.07, 0.22, 0.7),
                          p_female = 0.75, 
                          seed = i)
  }
}
return(data)
}


##### simulated data to baseline survival model

to_baseline <- function(data){
  data2 <- data
  data2 <- data2 %>% filter(time == 0)
  return(data2)
}
  
##### simulated data to TDCM

to_tdcm <- function(data){
  data2 <- data
  data2 <- data2 %>% group_by(id) %>% mutate(start = time, 
                                             stop = lead(start), 
                                             stop = ifelse(is.na(stop), start+0.5, stop),
                                             T = ifelse(Time <= stop, Time, stop),
                                             E = ifelse(Time <= stop, 1, 0),
                                             stop = ifelse(T <= stop, T, stop),
                                             event = E, bmi = last(bmi),
                                             clinician_disease=last(clinician_disease),
                                             patient_disease=last(patient_disease)) %>% filter(start < Time)
  return(data2)
}

##### simulated data to LM (aggregated)

to_lm_agg <- function(data, w, max_time){
  data2 <- data
  events <- data2 %>% filter(event == 1)
  event_times <- unique(events$Time)
  lm_points <- c(0, unique(round(sort(event_times, decreasing=FALSE), digits = 1)))
  lm_points <- unique(lm_points)
  lm_points <- lm_points[lm_points < max_time - w]
  new_data = list()
  for (i in 1:length(lm_points)) {
    landmark_data <- data2 %>% filter(Time > lm_points[i], time <= lm_points[i]) %>%
      group_by(id) %>% dplyr::summarise(median_bmi = median(bmi, na.rm=TRUE),
                                 min_bmi = min(bmi, na.rm=TRUE),
                                 max_bmi = max(bmi, na.rm=TRUE),
                                 sd_bmi = sd(as.numeric(bmi), na.rm=TRUE),
                                 median_clinician_disease = median(clinician_disease, na.rm=TRUE),
                                 min_clinician_disease = min(clinician_disease, na.rm=TRUE),
                                 max_clinician_disease = max(clinician_disease, na.rm=TRUE),
                                 sd_clinician_disease = sd(clinician_disease, na.rm=TRUE),
                                 median_patient_disease = median(patient_disease, na.rm=TRUE),
                                 min_patient_disease = min(patient_disease, na.rm=TRUE),
                                 max_patient_disease = max(patient_disease, na.rm=TRUE),
                                 sd_patient_disease = sd(patient_disease, na.rm=TRUE),
                                 Time = max(Time, na.rm = TRUE),
                                 event = max(event, na.rm = TRUE),
                                 disease_duration = max(disease_duration, na.rm = TRUE),
                                 sex = max(sex, na.rm = TRUE),
                                 age = max(age, na.rm = TRUE),
                                 bmi = last(bmi),
                                 clinician_disease=last(clinician_disease),
                                 patient_disease=last(patient_disease)) %>%
      mutate(sd_clinician_disease = ifelse(is.na(sd_clinician_disease), 0, sd_clinician_disease),
             sd_patient_disease = ifelse(is.na(sd_patient_disease), 0, sd_patient_disease),
             sd_bmi = ifelse(is.na(sd_bmi), 0, sd_bmi))
    
    landmark_data <- landmark_data %>% group_by(id) %>% mutate(event = ifelse(Time < lm_points[i] + w & event == 1, 1, 0),
                                                               Time = ifelse(Time > lm_points[i] + w & event == 0, lm_points[i] + w, Time),
                                                               LM = lm_points[i])
    new_data[[i]] <- landmark_data
  }
  return(new_data)
}

##### simulated data to TS-LM (ME model)

to_lm_long<-function(data, w, max_time){
  data2 <- data
  events <- data2 %>% filter(event == 1)
  event_times <- unique(events$Time)
  lm_points <- c(0, unique(round(sort(event_times, decreasing=FALSE), digits = 1)))
  lm_points <- unique(lm_points)
  lm_points <- lm_points[lm_points < max_time - w]
  new_data = list()
  for (i in 1:length(lm_points)) {
    landmark_data <- data2 %>% filter(Time > lm_points[i], time <= lm_points[i])
    landmark_data <- landmark_data %>% group_by(id) %>% mutate(event = ifelse(Time < lm_points[i] + w & event == 1, 1, 0),
                                                               Time = ifelse(Time > lm_points[i] + w & event == 0, lm_points[i] + w, Time),
                                                               LM = lm_points[i])
    new_data[[i]] <- landmark_data
    }
  return(new_data)
}

### temporal assessment datasets

prediction_datasets <- function(data, prediction_window, landmark_times, type){
  output = list()
  for (i in 1:length(landmark_times)){
    landmark_time = landmark_times[i]
    if (type == "short"){
      temp_data <- risk_set(data, prediction_window, landmark_time)
    }
    if (type == "long"){
      temp_data <- risk_set_long(data, prediction_window, landmark_time)
    }
    output[[i]] <- temp_data
    rm(temp_data)
  }
  return(output)
}

risk_set<-function(data, prediction_window, landmark_time){
  set <- data %>% dplyr::filter(Time > landmark_time, time <= landmark_time) %>%
    dplyr::group_by(id) %>% dplyr::summarise(median_bmi = median(bmi, na.rm=TRUE),
                                                   min_bmi = min(bmi, na.rm=TRUE),
                                                   max_bmi = max(bmi, na.rm=TRUE),
                                                   sd_bmi = sd(bmi, na.rm=TRUE),
                                                   median_clinician_disease = median(clinician_disease, na.rm=TRUE),
                                                   min_clinician_disease = min(clinician_disease, na.rm=TRUE),
                                                   max_clinician_disease = max(clinician_disease, na.rm=TRUE),
                                                   sd_clinician_disease = sd(clinician_disease, na.rm=TRUE),
                                                   median_patient_disease = median(patient_disease, na.rm=TRUE),
                                                   min_patient_disease = min(patient_disease, na.rm=TRUE),
                                                   max_patient_disease = max(patient_disease, na.rm=TRUE),
                                                   sd_patient_disease = sd(patient_disease, na.rm=TRUE),
                                                   Time = max(Time, na.rm = TRUE),
                                                   event = max(event, na.rm = TRUE),
                                                   disease_duration = max(disease_duration, na.rm = TRUE),
                                                   sex = max(sex, na.rm = TRUE),
                                                   age = max(age, na.rm = TRUE),
                                                   bmi = last(na.omit(bmi)),
                                                   clinician_disease=last(na.omit(clinician_disease)),
                                                   patient_disease=last(na.omit(patient_disease))) %>%
    mutate(sd_clinician_disease = ifelse(is.na(sd_clinician_disease), 0, sd_clinician_disease),
           sd_patient_disease = ifelse(is.na(sd_patient_disease), 0, sd_patient_disease),
           sd_bmi = ifelse(is.na(sd_bmi), 0, sd_bmi))
  
  set <- set %>% mutate(event = ifelse((Time <= landmark_time + prediction_window) & (event == 1), 1, 0),
           Time = ifelse(Time > landmark_time + prediction_window, landmark_time + prediction_window, Time),
           LM = landmark_time)
 return(set)
}

risk_set_long<-function(data, prediction_window, landmark_time){
  set <- data %>%
    filter(time <= landmark_time) %>%
    dplyr::group_by(id) %>%
  set <- set %>% filter(Time > landmark_time) %>% 
    mutate(event = ifelse((Time <= landmark_time + prediction_window) & (event == 1), 1, 0),
           Time = ifelse(Time > landmark_time + prediction_window, landmark_time + prediction_window, Time),
           LM = landmark_time)
  return(set)
}

