source("D:/Abheepsa_and_rashi/Stochastic Simulations v3/plot_theme.R")
source("D:/Abheepsa_and_rashi/Stochastic Simulations v3/plotFuncs.R")
source("D:/Abheepsa_and_rashi/Stochastic Simulations v3/simFuncs.R")

#################################################################################################
## Run simulations for all parameters in a paramter files and create csv file (with all initial conditions) for each parameter set.
## Then make plots for each paramters set with different trajectories in one plot.
# work_dir <- "D:/Abheepsa_and_rashi/Stochastic Simulations v3/input_files"
# setwd(work_dir)
# simulateFull("ESMR.topo", scaledNoise = F, deltaT = 0.1,noiseLevel=0 ,tMax = 200, nInit = 100, parallel = T, workers = 10)
# trajectories_full_folder("./ESMR/")

#################################################################################################
## Run simulations for all parameters in a paramter files and create csv file (with all initial conditions) for each parameter set.
## Then make plots for each paramters set with different trajectories in one plot.
work_dir <- "D:/Abheepsa_and_rashi/Stochastic Simulations v3/input_files"
setwd(work_dir)
#simulateFull("ESMR.topo", scaledNoise = T, deltaT = 0.1,noiseLevel=0.75 ,tMax = 200, nInit = 100, parallel = T, workers = 10)

landscape_full_folder("./100_0.75/",filenames_list=c(118,124,100,114,217,255)) # 

#################################################################################################
## Code to get multiple dynamic trajectories for different initial conditions for a list of parameter sets:
# 
# parameter_id_array <- c(217,255,100,114,118,124)
# 
# for( parameter_id in parameter_id_array){
# print(parameter_id)
# work_dir <- "D:/Abheepsa_and_rashi/Stochastic Simulations v3/input_files"
# setwd(work_dir)
# 
# topoFile <- "ESMR.topo"
# 
# simSwitching(topoFile, ParID = parameter_id,noiseLevel = 0.75,scaledNoise=T,nInit = 10,workers = 10)
# 
# dfTraj <- read_csv(paste0("./ESMR/",parameter_id, "_simDat.csv"))
# 
# dir.create(file.path("./../", paste0("scaled_noise_0.75/Parameter_id=",parameter_id)),recursive = TRUE)
# folder<- paste0("./../scaled_noise_0.75/Parameter_id=",parameter_id,"/")
# 
# ids <- unique(dfTraj$Id)
# sapply(ids, function(x){
#     trajectoriesSwitching(dfTraj, x, tMin = 0, tMax = NULL,folder) # use tMin and tMax for zoomed plots (make sure tMax is not NULL)
# })
# 
# }

#################################################################################################
## Code to get multiple dynamic trajectories for different initial conditions for a list of parameter sets:

# parameter_id_array <- c(118,124)
# 
# for( parameter_id in parameter_id_array){
#   print(parameter_id)
#   work_dir <- "D:/Abheepsa_and_rashi/Stochastic Simulations v2/input_files"
#   setwd(work_dir)
#   
#   topoFile <- "ERMS.topo"
#   
#   simSwitching(topoFile, ParID = parameter_id,noiseLevel = 1,scaledNoise=T,nInit = 20, tMax = 100,workers = 10)
#   
#   dfTraj <- read_csv(paste0("./ERMS/",parameter_id, "_simDat.csv"))
#   
#   dir.create(file.path("./../", paste0("scaled_noise_zoom_1/Parameter_id=",parameter_id)),recursive = TRUE)
#   folder<- paste0("./../scaled_noise_zoom_1/Parameter_id=",parameter_id,"/")
#   
#   ids <- unique(dfTraj$Id)
#   sapply(ids, function(x){
#     trajectoriesSwitching(dfTraj, x, tMin = 20, tMax = 60,folder) # use tMin and tMax for zoomed plots (make sure tMax is not NULL)
#   })
#   
# }
#################################################################################################
## Run simulations for all parameters in a paramter files and create csv file (with all initial conditions) for each parameter set.
## Then make landscape for each parameter set with different trajectories in one plot.
# work_dir <- "D:/Abheepsa_and_rashi/Stochastic Simulations v2/input_files"
# setwd(work_dir)
# simulateFull("ERMS.topo", paramIndex=c(118,124),scaledNoise = T, deltaT = 0.1,noiseLevel=0.5 ,tMax = 200, nInit = 100, parallel = T, workers = 10)
# landscape_full_folder("./ERMS/")

#################################################################################################






















