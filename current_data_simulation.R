
library(simstudy)
library(dplyr)



# Function to generate observations for each time point
generate_observation <- function(data, time, stage, scenario, obs, exchangeability){
  if (stage == 2 & exchangeability == 0){
    effect <- case_when(
      data$treatment == 1 ~ 0,                                     
      data$treatment == 2 ~ ifelse(scenario %in% c(4,5,6), (4 + 1)/4, ifelse(scenario %in% c(3), (2 + 1)/4, 1/4)),
      data$treatment == 3 ~ ifelse(scenario %in% c(2,3), (6 + 1)/4, ifelse(scenario %in% c(4,5,6), (8 + 1)/4, 1/4))
    )
  } else {
    effect <- case_when(
      data$treatment == 1 ~ 0,                                     
      data$treatment == 2 ~ ifelse(scenario %in% c(4,5,6), 4/4, ifelse(scenario %in% c(3), 2/4, 0)),
      data$treatment == 3 ~ ifelse(scenario %in% c(2,3), 6/4, ifelse(scenario %in% c(4,5,6), 8/4, 0))
    )
  }
  
  if (time == 1 & stage == 1){
    obs_value <- data$baseline + data$age + data$time_effect + rnorm(nrow(data), mean = effect, sd = 0.125)
  } else {
    obs_value <- obs + data$time_effect + rnorm(nrow(data), mean = effect, sd = 0.125)
  }
  return(data.frame(id = data$id, 
                    time = time, 
                    stage = stage, 
                    treatment = data$treatment, 
                    baseline = data$baseline, 
                    age = data$age, 
                    obs = obs_value))
}

for (N in c(25, 50, 250)){
  print(N)
  for (exchangeability in c(1, 0)){
    print(exchangeability)
  
    
    for (scenario in c(1:4)) { # Change as per scenario 
      
      print(scenario)
      
      datasets = list()
      
      for (iii in c(1:3000)) {
        
        # Generate longitudinal data
        long_data <- list()
        
        # Define baseline data
        def <- defData(varname = "age", formula = 0, variance = 0.25, dist = "normal")
        def <- defData(def, varname = "baseline", formula = 0, variance = 0.25, dist = "normal")
        def <- defData(def, varname = "treatment", formula = "0.2;0.4;0.4", dist = "categorical")
        def <- defData(def, varname = "time_effect", formula = 0, variance = 0.25, dist = "normal" )
        dt <- genData(N, def)
        
        for(stage in 1) {
          for(time in 1:4) {
            long_data[[length(long_data) + 1]] <- generate_observation(dt, time, stage, scenario, long_data[[length(long_data)]]$obs, exchangeability)
          }
        }
        
        # Add logic to categorize responders/non-responders at end of Stage 1
        threshold <- 4
        dt$response_stage1 <- ifelse(long_data[[4]]$obs >= threshold, 1, 0) # 1 for responders, 0 for non-responders
        
        # Logic for treatment assignment in Stage 2 based on response and initial treatment
        dt$treatment <- case_when(
          dt$treatment == 1 ~ sample(c(2,3), size = nrow(dt), replace = TRUE),            # Placebo to low/high dose
          dt$treatment == 2 & dt$response_stage1 == 1 ~ 2,                                # Low dose responders stay
          dt$treatment == 2 & dt$response_stage1 == 0 ~ 3,                                # Low dose non-responders to high dose
          dt$treatment == 3 & dt$response_stage1 == 1 ~ sample(c(2,3), size = nrow(dt), replace = TRUE), # High dose responders rerandomized
          dt$treatment == 3 & dt$response_stage1 == 0 ~ NA_real_                          # High dose non-responders leave
        )
        
        # Generate data for Stage 2 for the same scenario
        for(stage in 2) {
          for(time in 1:4) {
            long_data[[length(long_data) + 1]] <- generate_observation(dt, time, stage, scenario, long_data[[length(long_data)]]$obs, exchangeability)
          }
        }
        
        
        # Combine data into a single long format DataFrame
        longitudinal_data <- bind_rows(long_data)
        datasets[[iii]] = longitudinal_data
      }
      save(datasets, file = paste0("sample_size_", N, "_exchangeability_", exchangeability, "_scenario_", scenario, ".RData"))
    }
  }
}
