#!/usr/bin/env Rscript
setwd("/home/sidiwang/longitudinal")

args = commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
print(args)
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
} else {
  
  if (args[[1]] == 250) {
    load(paste0("simulated_datasets/external_sample_size_50_scenario_", args[[2]], ".RData"))
    external_data = datasets
  } else {
    load(paste0("simulated_datasets/external_sample_size_10_scenario_", args[[2]], ".RData"))
    external_data = datasets
  }
  
  if (args[[1]] == 25) {
    load(paste0("simulated_datasets/sample_size_25_exchangeability_", args[[3]], "_scenario_", args[[4]], ".RData"))
    current_data = datasets
  } else if (args[[1]] == 50){
    # load(paste0("external_sample_size_20_scenario_", args[[2]]))
    #  external_data = datasets
    load(paste0("simulated_datasets/sample_size_50_exchangeability_", args[[3]], "_scenario_", args[[4]], ".RData"))
    current_data = datasets
  } else {
    load(paste0("simulated_datasets/sample_size_250_exchangeability_", args[[3]], "_scenario_", args[[4]], ".RData"))
    current_data = datasets
  }
}

SIMULATION_N = 1000 #length(external_data)

# do simulation in parallel
require("foreach")
require("doParallel")
parallel::detectCores()
n.cores <- parallel::detectCores()
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

# check cluster definition (optional)
print(my.cluster)
# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# check if it is registered (optional)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

comb <- function(x, ...) {
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}

result <- foreach(
  i = 1:SIMULATION_N,
  .combine = comb,
  .inorder = TRUE,
  .multicombine = TRUE,
  .errorhandling='remove',
  .verbose = TRUE
) %dopar% {
  print(i)
  
  if (args[4] == "1") {
    pl <- 0
    ph <- 0
    theta_p_true = 0
    theta_1l_true = 0
    theta_1h_true = 0
  } else if (args[4] == "2") {
    pl <- 0
    ph <- 6
    theta_p_true = 0
    theta_1l_true = 0
    theta_1h_true = 6
  } else if (args[4] == "3") {
    pl <- 2
    ph <- 6
    theta_p_true = 0
    theta_1l_true = 2
    theta_1h_true = 6
  } else if (args[4] == "4") {
    pl <- 4
    ph <- 8
    theta_p_true = 0
    theta_1l_true = 4
    theta_1h_true = 8
  }
  
  library(readr)
  library(nlme)
  library(rjags)
  library(coda)
  library(stringr)
  library(tidyr)
  library(dplyr)
  
  external_subjects = external_data[[i]]
  current_subjects = current_data[[i]]
  external_subjects$id = paste0("e", external_subjects$id)
  current_subjects$id = paste0("c", current_subjects$id)
  current_subjects$stage = as.factor(current_subjects$stage)
  current_subjects$treatment = as.factor(current_subjects$treatment)
  
  external_wide <- external_subjects %>% 
    pivot_wider(
      id_cols = c(id, age, baseline),
      names_from = time,
      values_from = obs
    )
  
  if (args[[1]] == 25) {
    external_wide = external_wide[1:5,]
  }
  
  current_subject_stage1 = subset(current_subjects, stage == 1)
  
  current_wide <- current_subject_stage1 %>%
    pivot_wider(
      id_cols = c(id, age, baseline, treatment),
      names_from = time,
      values_from = obs
    )
  
  external_wide$study_indicator = 0
  current_wide$study_indicator = 1
  
  
  current_wide_placebo = current_wide #subset(current_wide, treatment == 1)
  
  
  IPTW_data = rbind(external_wide[, c("id", "age", "baseline", "study_indicator")], current_wide_placebo[, c("id", "age", "baseline", "study_indicator")])
  IPTW_data$study_indicator = as.factor(IPTW_data$study_indicator)
  
  ## Calculate Propensity Score
  IPTW_model = glm(study_indicator ~ age + baseline, data = IPTW_data, family = binomial("logit"))
  summary(IPTW_model)
  IPTW_data$propensity_score = 1
  IPTW_data$propensity_score[which(rowSums(is.na(IPTW_data)) == 0)] = predict(IPTW_model, type = "response", newdata = NULL)
  min_propensity_score = min(IPTW_data$propensity_score[which(IPTW_data$study_indicator == 1)])
  max_propensity_score = max(IPTW_data$propensity_score[which(IPTW_data$study_indicator == 1)])
  IPTW_data$propensity_score = ifelse(IPTW_data$study_indicator == 1, 1, IPTW_data$propensity_score)
  IPTW_data$weights = ifelse(IPTW_data$study_indicator == 1, 1, (IPTW_data$propensity_score) / (1 - IPTW_data$propensity_score))
  IPTW_data$weights[which(IPTW_data$study_indicator == 0)] = IPTW_data$weights[which(IPTW_data$study_indicator == 0)]/(sum(IPTW_data$weights[which(IPTW_data$study_indicator == 0)])/length(which(IPTW_data$study_indicator == 0)))
  IPTW_data$weights[which(IPTW_data$weights > 1)] = 1
  IPTW_data2 = IPTW_data
  to_delete = which(IPTW_data$propensity_score < min_propensity_score | IPTW_data$propensity_score > max_propensity_score & IPTW_data$propensity_score < 1)
  if (length(to_delete) != 0) {
    IPTW_data = IPTW_data[-which(IPTW_data$propensity_score < min_propensity_score),]
  }
  
  
  ### Model fitting
  # 1. Frequentist MMRM with data from 1st stage and continuous time (without External control)
  
  Frequentist <- gls(obs ~ baseline + treatment * time, na.action = na.omit, data = current_subject_stage1,
                     correlation = nlme::corSymm(form = ~ time | id),
                     weights = nlme::varIdent(form = ~ 1 | time))
  summary(Frequentist)
  
  # 2. BJSM with data from both stages, cross-sectional at week 48 (without External Control)
  current_BJSM = current_subjects
  current_BJSM$time = ifelse(current_BJSM$stage == 2, current_BJSM$time + 4, current_BJSM$time)
  current_BJSM_wide <- current_BJSM %>% 
    pivot_wider(
      id_cols = c(id, age, baseline),
      names_from = time,
      values_from = c(obs, treatment)
    )
  
  current_BJSM_wide = current_BJSM_wide[, c("id", "age", "baseline", "obs_4", "obs_8", "treatment_1", "treatment_5")]
  colnames(current_BJSM_wide)[6:7] = c("trt1", "trt2")
  
  current_BJSM_wide$bw48 = as.numeric(current_BJSM_wide$obs_4) - as.numeric(current_BJSM_wide$baseline)
  current_BJSM_wide$bw96 = as.numeric(current_BJSM_wide$obs_8) - as.numeric(current_BJSM_wide$obs_4)
  
  n_MCMC_chain = 2
  n.adapt = 10000
  MCMC_sample = 10000
  n.thin = 2
  
  # 6MWD
  posterior_sample_BJSM <- NULL
  attempt <- 1
  while (is.null(posterior_sample_BJSM) && attempt <= 5) {
    attempt <- attempt + 1
    
    try({
      jags_BJSM <- rjags::jags.model(
        file = "JointStageBayes_mixture.bug",
        data = list(
          overall_sample_size = nrow(current_BJSM_wide),
          stage2size = which(!is.na(current_BJSM_wide$bw96)),
          data_effect_stageI = current_BJSM_wide$bw48,
          data_effect_stageII = current_BJSM_wide$bw96,
          treatment_stageI = current_BJSM_wide$trt1,
          treatment_stageII = current_BJSM_wide$trt2,
          mu_guess = c(50, 50, 50),
          var_prior = 1/c(10000, 10000, 10000),
          mu_mixture = c(0, 0, 0),
          var_value = 0.0001
        ),
        n.chains = n_MCMC_chain, n.adapt = n.adapt
      )
      posterior_sample_BJSM <- rjags::coda.samples(
        jags_BJSM,
        c("mu", "alpha", "beta"),
        MCMC_sample, thin = n.thin
      )
    })
  }
  
  
  
  
  #### 3. Bayesian longitudinal piecewise model with data from 2 stages, continuous time.
  ## we need to create two groups of data to feed into the MCMC model
  # firstly we need the stage 1 Lilly trial + cureduchenne +  external control
  
  external_select = subset(external_wide, id %in% IPTW_data$id[which(str_detect(IPTW_data$id, "e"))])
  external_select$trt = as.factor(1)
  external_select$study_indicator = 0
  external_select$propensity_score = IPTW_data$weights[which(str_detect(IPTW_data$id, "e"))]
  current_wide$propensity_score = 1
  current_wide = current_wide[, c("id", "age", "baseline", "1", "2", "3", "4", "study_indicator", "treatment", "propensity_score")]
  colnames(current_wide)[9] = "trt"
  dataset1 = current_wide
  
  
  
  dataset1_long = reshape2::melt(dataset1, id.vars = c("id", "age", "trt", "study_indicator", "propensity_score", "baseline"),
                                 measure.vars = c("1", "2", "3", "4"),
                                 variable.name = "time",
                                 value.name = "outcome")
  
  dataset1_long = as.data.frame(dataset1_long %>% arrange(id, time))
  
  
  # secondly, we need the stage 2 Lilly trial + stage 1 last observations as the baseline for stage 2 Lilly trial
  current_subject_stage2 = subset(current_subjects, stage == 2)
  
  current_stage2 <- current_subject_stage2 %>%
    pivot_wider(
      id_cols = c(id, age, baseline, treatment),
      names_from = time,
      values_from = obs
    )
  
  current_stage2 = subset(current_stage2, id %in% current_wide$id)
  colnames(current_stage2)[4] = c("trt")
  dataset2 = current_stage2
  dataset2 = subset(dataset2, !is.na(trt))
  dataset2_long = reshape2::melt(dataset2, id.vars = c("id", "age", "trt", "baseline"),
                                 measure.vars = c("1", "2", "3", "4"),
                                 variable.name = "time",
                                 value.name = "outcome")
  
  dataset2_long = as.data.frame(dataset2_long %>% arrange(id, time))
  
  
  ### Calculate the MAC 
  
  external_long = reshape2::melt(external_select, id.vars = c("id", "age", "trt", "study_indicator", "propensity_score", "baseline"),
                                 measure.vars = c("1", "2", "3", "4"),
                                 variable.name = "time",
                                 value.name = "outcome")
  
  external_long = as.data.frame(external_long %>% arrange(id, time))
  
  
  dataset1_long = merge(dataset1_long, unique(dataset2_long[, c("id", "trt")]), by = "id", all.x = TRUE)
  colnames(dataset1_long) = c("id", "age", "trt1", "study_indicator", "propensity_score", "baseline", "time", "outcome", "trt2")
  dataset2_long = merge(dataset2_long, unique(dataset1_long[, c("id", "trt1")]), by = "id")
  dataset2_long$study_indicator = 1
  dataset2_long$propensity_score = 1
  dataset2_long = dataset2_long[, c("id", "age", "trt1", "study_indicator", "propensity_score", "baseline", "time", "outcome", "trt")]
  colnames(dataset2_long)[ncol(dataset2_long)] = "trt2" 
  dataset1_long$stage2 = 0
  dataset2_long$stage2 = 1
  dataset1_long$time1 = dataset1_long$time
  dataset1_long$time2 = 0
  dataset1_long$time = as.numeric(dataset1_long$time1) + dataset1_long$time2
  dataset2_long$time1 = 4
  dataset2_long$time2 = dataset2_long$time
  dataset2_long$time = dataset2_long$time1 + as.numeric(dataset2_long$time2)
  dataset = rbind(dataset1_long, dataset2_long)
  dataset$trt1 = as.numeric(dataset$trt1)
  dataset$trt2 = as.character(dataset$trt2)
  dataset$trt2 = ifelse(is.na(dataset$trt2), 1, dataset$trt2)
  dataset$trt2 = as.numeric(dataset$trt2)
  dataset$time1 = as.numeric(dataset$time1)
  dataset$time2 = as.numeric(dataset$time2)
  
  
  # now we can fit the jags model
  # NSAA
  n.adapt = ifelse(args[[1]] == 25, 200000, ifelse(args[[1]] == 50, 100000, 20000))
  MCMC_sample = ifelse(args[[1]] == 25, 2000000, ifelse(args[[1]] == 50, 1000000, 500000))
  n.burnin = ifelse(args[[1]] == 25, 200000, ifelse(args[[1]] == 50, 100000, 20000))
  n_MCMC_chain = 2
  n.thin = ifelse(args[[1]] == 25, 200, ifelse(args[[1]] == 50, 100, 50))
  
  posterior_sample_BLPM <- NULL
  attempt <- 1
  while (is.null(posterior_sample_BLPM) && attempt <= 5) {
    attempt <- attempt + 1
    
    try({
      BLPM <- rjags::jags.model(
        file = "BLPM_test.bug",
        data = list(
          
          # external control
          N_subj_ext = length(unique(external_long$id)),
          N_obs_ext = nrow(external_long),
          age_ext = external_long$age,
          y_ext = external_long$outcome,
          baseline_ext = external_long$baseline,
          time_ext = external_long$time,
          subj_ext = as.numeric(as.factor(external_long$id)),
          W_ext = (external_long$propensity_score),
          
          # current trial
          time = dataset$time,
          time1 = dataset$time1,
          time2 = dataset$time2,
          age = dataset$age,
          trt1 = dataset$trt1,
          trt2 = dataset$trt2,
          baseline = dataset$baseline,
          y = dataset$outcome,
          N_subj = length(unique(dataset$id)),
          N_obs = nrow(dataset),
          subj = as.numeric(as.factor(dataset$id)),
          stage2 = dataset$stage2,
          p.exch = 0.5, 
          Prior.tau_0 = 0.125,
          Prior.tau_1 = 0.125,
          Prior.mu_0 = c(0, 5),
          Prior.mu_1 = c(0, 5),
          Prior.nex_0 = 0,
          Prior.nex_1 = 0,
          mu_nex_0 = 0,
          mu_nex_1 = 0,
          Prior.tau_treatment = 0.125
        ),
        inits = function(){
          list(beta = c(0, 0),
               beta_age = -1,
               beta_baseline = 1,
               beta_treatment = c(NA, pl/4, ph/4),
               mu_new = c(NA, pl/4, ph/4)
          )
        },
        n.chains = n_MCMC_chain, n.adapt = n.adapt
      )
      update(BLPM, n.iter = n.burnin)
      posterior_sample_BLPM <- rjags::coda.samples(
        BLPM,
        c("beta_ext", "beta_age_ext", "beta_baseline_ext", 
          "sigma_b0_ext", "sigma_b1_ext", "sigma_b0_log_ext", "sigma_b1_log_ext", "sigma_y_ext",
          "beta", "beta_baseline", "beta_age", "beta_treatment", 
          "beta_treatment2", "Z_current", "mu_new",
          "Z", "sigma_b0", "sigma_b1", "y_external_tau", 
          "sigma_y", "mu_nex_0", "mu_nex_1", "Prior.nex_0", "Prior.nex_1",
          "tau_0", "tau_1", "tau_new", "tau_new2", "PH", "PL", "mu_0", "mu_1",
          "tau_0_beta", "tau_1_beta"),
        n.iter = MCMC_sample, thin = n.thin
      )
    })
  }
  
  summary(posterior_sample_BLPM)[[1]]
  
  
  ## 4. Bayesian Longitudinal Piecewise Model with robust treatment effect estimation
  
  posterior_sample_BLPM_robust <- NULL
  attempt <- 1
  while (is.null(posterior_sample_BLPM_robust) && attempt <= 5) {
    attempt <- attempt + 1
    
    try({
      BLPM_robust <- rjags::jags.model(
        file = "BLPM_robust_test.bug",
        data = list(
          
          
          # external control
          N_subj_ext = length(unique(external_long$id)),
          N_obs_ext = nrow(external_long),
          age_ext = external_long$age,
          y_ext = external_long$outcome,
          baseline_ext = external_long$baseline,
          time_ext = external_long$time,
          subj_ext = as.numeric(as.factor(external_long$id)),
          W_ext = (external_long$propensity_score),
          
          # current trial
          time = dataset$time,
          time1 = dataset$time1,
          time2 = dataset$time2,
          age = dataset$age,
          trt1 = dataset$trt1,
          trt2 = dataset$trt2,
          baseline = dataset$baseline,
          y = dataset$outcome,
          N_subj = length(unique(dataset$id)),
          N_obs = nrow(dataset),
          subj = as.numeric(as.factor(dataset$id)),
          stage2 = dataset$stage2,
          p.exch = 0.5,
          p.treatment = 0.5, 
          #N_ExternalControl = 1,
          Prior.tau_0 = 0.125,
          Prior.tau_1 = 0.125,
          Prior.mu_0 = c(0, 5),
          Prior.mu_1 = c(0, 5),
          Prior.nex_0 = 0,
          Prior.nex_1 = 0,
          mu_nex_0 = 0,
          mu_nex_1 = 0,
          mu_new_nex = c(0,
                         mean(dataset$outcome[which(dataset$trt1 == 2 & dataset$stage2 == 0 & dataset$time == 4)])/4,
                         mean(dataset$outcome[which(dataset$trt1 == 3 & dataset$stage2 == 0 & dataset$time == 4)])/4),
          mu_new2_nex = matrix(c(rep(0, 3), 
                                 rep(mean(dataset$outcome[which(dataset$trt1 == 2 & dataset$stage2 == 0 & dataset$time == 4)])/4, 3), 
                                 rep(mean(dataset$outcome[which(dataset$trt1 == 3 & dataset$stage2 == 0 & dataset$time == 4)])/4, 3)), nrow = 3),
          Prior.tau_treatment = 0.125
        ),
        inits = function(){
          list(beta = c(0, 0),
               beta_age = -1,
               beta_baseline = 1,
               mu_new = c(NA, pl/4, ph/4)
          )
        },
        n.chains = n_MCMC_chain, n.adapt = n.adapt
      )
      update(BLPM_robust, n.burnin)
      posterior_sample_BLPM_robust <- rjags::coda.samples(
        BLPM_robust,
        c("beta_ext", "beta_age_ext", "beta_baseline_ext", 
          "sigma_b0_ext", "sigma_b1_ext", "sigma_b0_log_ext", "sigma_b1_log_ext", 
          "sigma_y_ext", "beta", "beta_baseline", "beta_age", "beta_treatment",
          "beta_treatment2", "Z_current", "mu_new", "mu_0", "mu_1",
          "Z", "sigma_b0", "sigma_b1", "Z_treatment", 
          "Z_treatment2", "y_external_tau", "tau_new", "tau_new2", "PH", "PL"),
        MCMC_sample, thin = n.thin
      )
    })
  }
  
  summary(posterior_sample_BLPM_robust)[[1]]
  
  
  
  
  if (args[4] == "1") {
    pl <- 0
    ph <- 0
    theta_p_true = 0
    theta_1l_true = 0
    theta_1h_true = 0
  } else if (args[4] == "2") {
    pl <- 0
    ph <- 6
    theta_p_true = 0
    theta_1l_true = 0
    theta_1h_true = 6
  } else if (args[4] == "3") {
    pl <- 2
    ph <- 6
    theta_p_true = 0
    theta_1l_true = 2
    theta_1h_true = 6
  } else if (args[4] == "4") {
    pl <- 4
    ph <- 8
    theta_p_true = 0
    theta_1l_true = 4
    theta_1h_true = 8
  }
  
  print(i)
  library(readr)
  library(nlme)
  library(rjags)
  library(coda)
  library(stringr)
  library(tidyr)
  library(dplyr)
  
  # frequentist
  coefs = coef(Frequentist)
  cov_mat = vcov(Frequentist)
  visitnum_value = 4
  
  sum_coefs_12 = coefs["treatment2"] + coefs["treatment2:time"] * visitnum_value
  SE_coef1_12 = sqrt(cov_mat["treatment2", "treatment2"])
  SE_coef2_12 = sqrt(cov_mat["treatment2:time", "treatment2:time"])
  cov_coef1_coef2_12 = cov_mat["treatment2", "treatment2:time"]
  SE_sum_coefs_12 = sqrt(SE_coef1_12^2 + (visitnum_value^2) * SE_coef2_12^2 + 2 * visitnum_value * cov_coef1_coef2_12)
  
  sum_coefs_13 = coefs["treatment3"] + coefs["treatment3:time"] * visitnum_value
  SE_coef1_13 = sqrt(cov_mat["treatment3", "treatment3"])
  SE_coef2_13 = sqrt(cov_mat["treatment3:time", "treatment3:time"])
  cov_coef1_coef2_13 = cov_mat["treatment3", "treatment3:time"]
  SE_sum_coefs_13 = sqrt(SE_coef1_13^2 + (visitnum_value^2) * SE_coef2_13^2 + 2 * visitnum_value * cov_coef1_coef2_13)
  
  frequentist_outcome_i = c(sum_coefs_12, SE_sum_coefs_12, sum_coefs_13, SE_sum_coefs_13, pl, ph)
  names(frequentist_outcome_i) = c("pl_est", "pl_sd", "ph_est", "ph_sd", "true_pl", "true_ph")
  
  # BJSM
  
  BJSM <- as.data.frame(posterior_sample_BJSM[[1]])
  BJSM$PL <- BJSM$`mu[2]` - BJSM$`mu[1]`
  BJSM$PH <- BJSM$`mu[3]` - BJSM$`mu[1]`
  
  BJSM_mean_estimate <- colMeans(BJSM)
  BJSM_tmp_mean <- c(BJSM_mean_estimate["mu[1]"], BJSM_mean_estimate["mu[2]"], BJSM_mean_estimate["mu[3]"], BJSM_mean_estimate["PL"], BJSM_mean_estimate["PH"], pl, ph)
  names(BJSM_tmp_mean)[6:7] = c("true_pl", "true_ph")
  BJSM_tmp_response_rate_posterior_mean <- colMeans(BJSM[, c("mu[1]", "mu[2]", "mu[3]", "PL", "PH")])
  BJSM_tmp_hdi <- HDInterval::hdi(BJSM, 0.95)
  
  BJSM_tmp_response_rate_hdi_coverage_rate_theta_p <- as.numeric(BJSM_tmp_hdi["lower", "mu[1]"] <= theta_p_true & BJSM_tmp_hdi["upper", "mu[1]"] >= theta_p_true)
  BJSM_tmp_response_rate_hdi_coverage_rate_theta_l <- as.numeric(BJSM_tmp_hdi["lower", "mu[2]"] <= theta_1l_true & BJSM_tmp_hdi["upper", "mu[2]"] >= theta_1l_true)
  BJSM_tmp_response_rate_hdi_coverage_rate_theta_h <- as.numeric(BJSM_tmp_hdi["lower", "mu[3]"] <= theta_1h_true & BJSM_tmp_hdi["upper", "mu[3]"] >= theta_1h_true)
  BJSM_tmp_response_rate_hdi_coverage_rate_pl <- as.numeric(BJSM_tmp_hdi["lower", "PL"] <= 0 & BJSM_tmp_hdi["upper", "PL"] >= 0)
  BJSM_tmp_response_rate_hdi_coverage_rate_ph <- as.numeric(BJSM_tmp_hdi["lower", "PH"] <= 0 & BJSM_tmp_hdi["upper", "PH"] >= 0)
  BJSM_tmp_response_rate_hdi_coverage_rate_pl_2 <- as.numeric(BJSM_tmp_hdi["lower", "PL"] <= pl & BJSM_tmp_hdi["upper", "PL"] >= pl)
  BJSM_tmp_response_rate_hdi_coverage_rate_ph_2 <- as.numeric(BJSM_tmp_hdi["lower", "PH"] <= ph & BJSM_tmp_hdi["upper", "PH"] >= ph)
  BJSM_tmp_response_rate_hdi_coverage_rate <- cbind(
    BJSM_tmp_response_rate_hdi_coverage_rate_theta_p, BJSM_tmp_response_rate_hdi_coverage_rate_theta_l,
    BJSM_tmp_response_rate_hdi_coverage_rate_theta_h, BJSM_tmp_response_rate_hdi_coverage_rate_pl,
    BJSM_tmp_response_rate_hdi_coverage_rate_ph, BJSM_tmp_response_rate_hdi_coverage_rate_pl_2,
    BJSM_tmp_response_rate_hdi_coverage_rate_ph_2
  )
  
  BJSM_tmp_response_rate_hdi_length_theta_p <- abs(BJSM_tmp_hdi["lower", "mu[1]"] - BJSM_tmp_hdi["upper", "mu[1]"])
  BJSM_tmp_response_rate_hdi_length_theta_l <- abs(BJSM_tmp_hdi["lower", "mu[2]"] - BJSM_tmp_hdi["upper", "mu[2]"])
  BJSM_tmp_response_rate_hdi_length_theta_h <- abs(BJSM_tmp_hdi["lower", "mu[3]"] - BJSM_tmp_hdi["upper", "mu[3]"])
  BJSM_tmp_response_rate_hdi_length_pl <- abs(BJSM_tmp_hdi["lower", "PL"] - BJSM_tmp_hdi["upper", "PL"])
  BJSM_tmp_response_rate_hdi_length_ph <- abs(BJSM_tmp_hdi["lower", "PH"] - BJSM_tmp_hdi["upper", "PH"])
  BJSM_tmp_response_rate_hdi_length <- cbind(
    BJSM_tmp_response_rate_hdi_length_theta_p, BJSM_tmp_response_rate_hdi_length_theta_l,
    BJSM_tmp_response_rate_hdi_length_theta_h, BJSM_tmp_response_rate_hdi_length_pl,
    BJSM_tmp_response_rate_hdi_length_ph
  )
  
  # BLPM
  
  BLPM = summary(posterior_sample_BLPM)[[1]]
  BLPM_tmp_mean <- c(BLPM["PL", 1], BLPM["PH", 1], pl, ph)
  names(BLPM_tmp_mean) = c("PL", "PH", "true_pl", "true_ph")
  BLPM_tmp_response_rate_posterior_mean <- c(BLPM["PL", 1], BLPM["PH", 1])
  names(BLPM_tmp_response_rate_posterior_mean) = c("PL", "PH")
  BLPM_tmp_hdi <- HDInterval::hdi(posterior_sample_BLPM, 0.95)
  
  BLPM_tmp_response_rate_hdi_coverage_rate_pl <- as.numeric(BLPM_tmp_hdi["lower", "PL"] <= 0 & BLPM_tmp_hdi["upper", "PL"] >= 0)
  BLPM_tmp_response_rate_hdi_coverage_rate_ph <- as.numeric(BLPM_tmp_hdi["lower", "PH"] <= 0 & BLPM_tmp_hdi["upper", "PH"] >= 0)
  BLPM_tmp_response_rate_hdi_coverage_rate_pl_2 <- as.numeric(BLPM_tmp_hdi["lower", "PL"] <= pl & BLPM_tmp_hdi["upper", "PL"] >= pl)
  BLPM_tmp_response_rate_hdi_coverage_rate_ph_2 <- as.numeric(BLPM_tmp_hdi["lower", "PH"] <= ph & BLPM_tmp_hdi["upper", "PH"] >= ph)
  BLPM_tmp_response_rate_hdi_coverage_rate <- cbind(
    BLPM_tmp_response_rate_hdi_coverage_rate_pl,
    BLPM_tmp_response_rate_hdi_coverage_rate_ph, BLPM_tmp_response_rate_hdi_coverage_rate_pl_2,
    BLPM_tmp_response_rate_hdi_coverage_rate_ph_2
  )
  
  BLPM_tmp_response_rate_hdi_length_pl <- abs(BLPM_tmp_hdi["lower", "PL"] - BLPM_tmp_hdi["upper", "PL"])
  BLPM_tmp_response_rate_hdi_length_ph <- abs(BLPM_tmp_hdi["lower", "PH"] - BLPM_tmp_hdi["upper", "PH"])
  BLPM_tmp_response_rate_hdi_length <- cbind(
    BLPM_tmp_response_rate_hdi_length_pl,
    BLPM_tmp_response_rate_hdi_length_ph
  )
  
  
  # BLPM robust
  
  BLPM_robust = summary(posterior_sample_BLPM_robust)[[1]]
  BLPM_robust_tmp_mean <- c(BLPM_robust["PL", 1], BLPM_robust["PH", 1], pl, ph)
  names(BLPM_robust_tmp_mean) = c("PL", "PH", "true_pl", "true_ph")
  BLPM_robust_tmp_response_rate_posterior_mean <- c(BLPM_robust["PL", 1], BLPM_robust["PH", 1])
  names(BLPM_robust_tmp_response_rate_posterior_mean) = c("PL", "PH")
  BLPM_robust_tmp_hdi <- HDInterval::hdi(posterior_sample_BLPM_robust, 0.95)
  
  BLPM_robust_tmp_response_rate_hdi_coverage_rate_pl <- as.numeric(BLPM_robust_tmp_hdi["lower", "PL"] <= 0 & BLPM_robust_tmp_hdi["upper", "PL"] >= 0)
  BLPM_robust_tmp_response_rate_hdi_coverage_rate_ph <- as.numeric(BLPM_robust_tmp_hdi["lower", "PH"] <= 0 & BLPM_robust_tmp_hdi["upper", "PH"] >= 0)
  BLPM_robust_tmp_response_rate_hdi_coverage_rate_pl_2 <- as.numeric(BLPM_robust_tmp_hdi["lower", "PL"] <= pl & BLPM_robust_tmp_hdi["upper", "PL"] >= pl)
  BLPM_robust_tmp_response_rate_hdi_coverage_rate_ph_2 <- as.numeric(BLPM_robust_tmp_hdi["lower", "PH"] <= ph & BLPM_robust_tmp_hdi["upper", "PH"] >= ph)
  BLPM_robust_tmp_response_rate_hdi_coverage_rate <- cbind(
    BLPM_robust_tmp_response_rate_hdi_coverage_rate_pl,
    BLPM_robust_tmp_response_rate_hdi_coverage_rate_ph, BLPM_robust_tmp_response_rate_hdi_coverage_rate_pl_2,
    BLPM_robust_tmp_response_rate_hdi_coverage_rate_ph_2
  )
  
  BLPM_robust_tmp_response_rate_hdi_length_pl <- abs(BLPM_robust_tmp_hdi["lower", "PL"] - BLPM_robust_tmp_hdi["upper", "PL"])
  BLPM_robust_tmp_response_rate_hdi_length_ph <- abs(BLPM_robust_tmp_hdi["lower", "PH"] - BLPM_robust_tmp_hdi["upper", "PH"])
  BLPM_robust_tmp_response_rate_hdi_length <- cbind(
    BLPM_robust_tmp_response_rate_hdi_length_pl,
    BLPM_robust_tmp_response_rate_hdi_length_ph
  )
  
  
  return(list(
    "frequentist_outcome" = frequentist_outcome_i, #frequentist_outcome,
    "BJSM_final_mean" = BJSM_tmp_mean, #BJSM_final_mean,
    "BJSM_response_rate_hdi_coverage_rate" = BJSM_tmp_response_rate_hdi_coverage_rate, #BJSM_response_rate_hdi_coverage_rate,
    "BJSM_response_rate_hdi_length" = BJSM_tmp_response_rate_hdi_length, #BJSM_response_rate_hdi_length,
    "BJSM_response_rate_posterior_mean" = BJSM_tmp_response_rate_posterior_mean, #BJSM_response_rate_posterior_mean,
    "BLPM_final_mean" = BLPM_tmp_mean, #BLPM_final_mean,
    "BLPM_response_rate_hdi_coverage_rate" = BLPM_tmp_response_rate_hdi_coverage_rate, #BLPM_response_rate_hdi_coverage_rate,
    "BLPM_response_rate_hdi_length" = BLPM_tmp_response_rate_hdi_length, # BLPM_response_rate_hdi_length,
    "BLPM_response_rate_posterior_mean" = BLPM_tmp_response_rate_posterior_mean, #BLPM_response_rate_posterior_mean,
    "BLPM_robust_final_mean" = BLPM_robust_tmp_mean, #BLPM_robust_final_mean,
    "BLPM_robust_response_rate_hdi_coverage_rate" = BLPM_robust_tmp_response_rate_hdi_coverage_rate, #BLPM_robust_response_rate_hdi_coverage_rate,
    "BLPM_robust_response_rate_hdi_length" = BLPM_robust_tmp_response_rate_hdi_length, #BLPM_robust_response_rate_hdi_length,
    "BLPM_robust_response_rate_posterior_mean" = BLPM_robust_tmp_response_rate_posterior_mean #BLPM_robust_response_rate_posterior_mean
  ))
}

# combine result from parallel for loops
frequentist_outcome <- result$frequentist_outcome
BJSM_final_mean <- result$BJSM_final_mean
BJSM_response_rate_hdi_coverage_rate <- result$BJSM_response_rate_hdi_coverage_rate
BJSM_response_rate_hdi_length <- result$BJSM_response_rate_hdi_length
BJSM_response_rate_posterior_mean <- result$BJSM_response_rate_posterior_mean
BLPM_final_mean <- result$BLPM_final_mean
BLPM_response_rate_hdi_coverage_rate <- result$BLPM_response_rate_hdi_coverage_rate
BLPM_response_rate_hdi_length <- result$BLPM_response_rate_hdi_length
BLPM_response_rate_posterior_mean <- result$BLPM_response_rate_posterior_mean
BLPM_robust_final_mean <- result$BLPM_robust_final_mean
BLPM_robust_response_rate_hdi_coverage_rate <- result$BLPM_robust_response_rate_hdi_coverage_rate
BLPM_robust_response_rate_hdi_length <- result$BLPM_robust_response_rate_hdi_length
BLPM_robust_response_rate_posterior_mean <- result$BLPM_robust_response_rate_posterior_mean

stopCluster(my.cluster)


frequentist_bias_pl = mean(frequentist_outcome[, 1] - frequentist_outcome[, 5])
frequentist_bias_ph = mean(frequentist_outcome[, 3] - frequentist_outcome[, 6])

frequentist_rmse_pl = sqrt(mean((frequentist_outcome[, 1] - frequentist_outcome[, 5])^2))
frequentist_rmse_ph = sqrt(mean((frequentist_outcome[, 3] - frequentist_outcome[, 6])^2))

frequentist_coverage_pl = mean(frequentist_outcome[, 1] - 1.96 * frequentist_outcome[, 2] < frequentist_outcome[, 5] & frequentist_outcome[, 1] + 1.96 * frequentist_outcome[, 2] > frequentist_outcome[, 5])
frequentist_coverage_ph = mean(frequentist_outcome[, 3] - 1.96 * frequentist_outcome[, 4] < frequentist_outcome[, 6] & frequentist_outcome[, 3] + 1.96 * frequentist_outcome[, 4] > frequentist_outcome[, 6])

frequentist_width_pl = mean(frequentist_outcome[, 2] * 1.96 * 2)
frequentist_width_ph = mean(frequentist_outcome[, 4] * 1.96 * 2)

BJSM_bias_pl = mean(BJSM_final_mean[, 4] - BJSM_final_mean[, 6])
BJSM_bias_ph = mean(BJSM_final_mean[, 5] - BJSM_final_mean[, 7])

BJSM_rmse_pl = sqrt(mean((BJSM_final_mean[, 4] - BJSM_final_mean[, 6])^2))
BJSM_rmse_ph = sqrt(mean((BJSM_final_mean[, 5] - BJSM_final_mean[, 7])^2))

BJSM_coverage_pl = colMeans(BJSM_response_rate_hdi_coverage_rate)[6]
BJSM_coverage_ph = colMeans(BJSM_response_rate_hdi_coverage_rate)[7]

BJSM_width_pl = colMeans(BJSM_response_rate_hdi_length)[4]
BJSM_width_ph = colMeans(BJSM_response_rate_hdi_length)[5]

BLPM_bias_pl = mean(BLPM_final_mean[, 1] - BLPM_final_mean[, 3])
BLPM_bias_ph = mean(BLPM_final_mean[, 2] - BLPM_final_mean[, 4])

BLPM_rmse_pl = sqrt(mean((BLPM_final_mean[, 1] - BLPM_final_mean[, 3])^2))
BLPM_rmse_ph = sqrt(mean((BLPM_final_mean[, 2] - BLPM_final_mean[, 4])^2))

BLPM_coverage_pl = colMeans(BLPM_response_rate_hdi_coverage_rate)[3]
BLPM_coverage_ph = colMeans(BLPM_response_rate_hdi_coverage_rate)[4]

BLPM_width_pl = colMeans(BLPM_response_rate_hdi_length)[1]
BLPM_width_ph = colMeans(BLPM_response_rate_hdi_length)[2]

BLPM_robust_bias_pl = mean(BLPM_robust_final_mean[, 1] - BLPM_robust_final_mean[, 3])
BLPM_robust_bias_ph = mean(BLPM_robust_final_mean[, 2] - BLPM_robust_final_mean[, 4])

BLPM_robust_rmse_pl = sqrt(mean((BLPM_robust_final_mean[, 1] - BLPM_robust_final_mean[, 3])^2))
BLPM_robust_rmse_ph = sqrt(mean((BLPM_robust_final_mean[, 2] - BLPM_robust_final_mean[, 4])^2))

BLPM_robust_coverage_pl = colMeans(BLPM_robust_response_rate_hdi_coverage_rate)[3]
BLPM_robust_coverage_ph = colMeans(BLPM_robust_response_rate_hdi_coverage_rate)[4]

BLPM_robust_width_pl = colMeans(BLPM_robust_response_rate_hdi_length)[1]
BLPM_robust_width_ph = colMeans(BLPM_robust_response_rate_hdi_length)[2]


frequentist_bias_pl
frequentist_bias_ph
BJSM_bias_pl
BJSM_bias_ph
BLPM_bias_pl
BLPM_bias_ph
BLPM_robust_bias_pl
BLPM_robust_bias_ph

frequentist_rmse_pl
frequentist_rmse_ph
BJSM_rmse_pl
BJSM_rmse_ph
BLPM_rmse_pl
BLPM_rmse_ph
BLPM_robust_rmse_pl
BLPM_robust_rmse_ph

frequentist_coverage_pl
frequentist_coverage_ph
BJSM_coverage_pl
BJSM_coverage_ph
BLPM_coverage_pl
BLPM_coverage_ph
BLPM_robust_coverage_pl
BLPM_robust_coverage_ph

frequentist_width_pl
frequentist_width_ph
BJSM_width_pl
BJSM_width_ph
BLPM_width_pl
BLPM_width_ph
BLPM_robust_width_pl
BLPM_robust_width_ph

# save outcome
save.image(file = paste0("/home/sidiwang/longitudinal/result/test_result_", args[[1]], "_", args[[2]], "_", args[[3]], "_", args[[4]], "_outcome.RData"))
