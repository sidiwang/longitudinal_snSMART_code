
library(simstudy)
library(dplyr)



# Function to generate observations for each time point
generate_observation <- function(data, time, scenario, obs){

  if (time == 1){
    obs_value <- data$baseline + data$age + data$time_effect + rnorm(nrow(data), mean = dt$scenario, sd = 0.125)
  } else {
    obs_value <- obs + data$time_effect + rnorm(nrow(data), mean = dt$scenario, sd = 0.125)
  }
  return(data.frame(id = data$id, 
                    time = time, 
                    baseline = data$baseline, 
                    age = data$age, 
                    obs = obs_value))
}

for (N in c(5, 10, 20, 50)){
  print(N)
    
    
    for (scenario in c(1:3)) { # Change as per scenario 
      
      print(scenario)
      
      datasets = list()
      
      for (iii in c(1:3000)) {
        
        # Generate longitudinal data
        long_data <- list()
        
        # Define baseline data
        def <- defData(varname = "age", formula = 0, variance = 0.25, dist = "normal")
        def <- defData(def, varname = "baseline", formula = 0, variance = 0.25, dist = "normal")
        def <- defData(def, varname = "time_effect", formula = 0, variance = 0.25, dist = "normal" )
        dt <- genData(N, def)
        dt$scenario = rbinom(N, 1, ifelse(scenario == 1, 0, ifelse(scenario == 2, 0.5, 1))) * 0.25
        
        for(stage in 1) {
          for(time in 1:4) {
            long_data[[length(long_data) + 1]] <- generate_observation(dt, time, scenario, long_data[[length(long_data)]]$obs)
          }
        }
        

        # Combine data into a single long format DataFrame
        longitudinal_data <- bind_rows(long_data)
        datasets[[iii]] = longitudinal_data
      }
      save(datasets, file = paste0("test_external_", "sample_size_", N, "_scenario_", scenario, ".RData"))
    }
  }

