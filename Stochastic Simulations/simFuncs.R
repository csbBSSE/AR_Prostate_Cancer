library(tidyverse)
options(stringsAsFactors = F, warn = -1)
library(future.apply)

### For a full description of the parameters, see the "simulateFull" function.

## calculates the hill shift function for the given parameters
hillShift <- function(x, x0y, n, l)
{
    # print(l)
    hP <- x^n/(x^n + x0y^n)
    hS <- (1-hP) + l*hP
    if (l>1)
    {
        hS <- hS/l
    }
    return(hS)
}

# Reads the RACIPE parameters corresponding to a topo File
paramRead <- function(topoFile, prs = T)
{
    net <- str_remove(topoFile, ".topo")
    if (!prs)
        parameters <- read.delim(paste0(net, "_parameters.dat"))
    else
    {
        prs <- read.delim(paste0(net, ".prs"), stringsAsFactors = F)
        paramOrder <- prs$Parameter
        parameters <- read.delim(paste0(net, "_parameters.dat"), header = F)
        colnames(parameters) <- c("ID", "nStates", paramOrder)
    }
    
    return(parameters)
}


# Generates an discrete time dynamic state updater for the network in the topo file
#   Takes in the input of a topofile, a set of named parameters, a dummy state vector, 
#noiseLevel and scaledNoise (boolean). Constructs a function that evaluates the next state
#using the provided parameter sets given a state vector.
# Note: the generated function retains the parameter values and the equations in its environment and
#   only needs a state vector and deltaT as input.
emFunc <- function(topoFile, Params, state, deltaT, noiseLevel, scaledNoise = F)
{
    print("reached emFunc")
    df <- read.delim(topoFile, sep = "")
    nodes <- unique(c(df$Source, df$Target)) %>% sort
    update <- function(state, deltaT)
    {
        stateNext <- sapply(nodes, function(node){
            s0 <- state[node]
            noise <- noiseLevel*rnorm(1)*sqrt(deltaT)
            if (scaledNoise)
            {
                noise <- noise*s0
            }
            s1 <- s0 + noise - Params[paste0("Deg_of_", node)]*s0*deltaT
            links <- which(df$Target == node)
            hills <- apply(df[links,], 1, function(link){
                source <- link[1]
                target <- link[2]
                type <- link[3]
                x <- state[source]
                n <- Params[paste0("Num_of_", source, "To", target)]
                x0y <- Params[paste0("Trd_of_", source, "To", target)]
                lKey <- ifelse(type == "1", "Act_of_", "Inh_of_") %>% paste0(., source, "To", target)
                l <- Params[lKey]
                hillShift(x, x0y, n,l)
            })
            # print(node)
            production <- Params[paste0("Prod_of_", node)]*prod(hills)
            s1 <- s1 + production*deltaT
            return(s1)
        })
        names(stateNext) <- nodes
        return(stateNext)
    }
    return(update)
}


# simulates the network in a topo file for a specific set of parameters and initial conditions given the update function
# returns a trajectory dataframe with time and node expression levels
simulateInitCond <- function(topoFile, Params, initState, deltaT, tMax, updFunc)
{
    print("reached simulateInitCond")
    # browser()
    df <- read.delim(topoFile, sep = "")
    nodes <- unique(c(df$Source, df$Target)) %>% sort
    tVec <- seq(0, tMax, deltaT)
    tSteps <- 2:length(tVec)
    state <- initState
    trajectory <- sapply(tSteps, function(t){
        state <<- updFunc(state, deltaT)
        return(state)
    }) %>% t
    trajectory <- rbind.data.frame(initState, trajectory) %>% set_names(nodes) %>%
        mutate(Time = tVec)
    return(trajectory)
    
}

# Simulates the network for a parameter set for multiple initial conditions
# For a given topofile and parameter set, simulates the network stochastically for nInit random initial conditions

## Steps:
#1. read the topo file into a data frame
#2. obtain the maximum value for each node expression level from which 
#   the initial conditions are to be sampled (sMax). The level used is 1.5*production rate/degradation rate
#3. generate the update function corresponding to the parameter set, using a dummy state vector
#4. generate a random initial condition sampled uniformly from the range [0, sMax].
#5. calculate the trajectory for that initial condition using the function simulateInitCond
#6. Repeat 4 and 5 nInit times. 
#7. rbind all the trajectories
#8. log-normalize the expression values using mean and SD from RACIPE.
#9. calculate EMscore and Resistance score for each row. store the outcome

## Note: For normalization, you need the mean and SD for each node expression in 
### RACIPE. In this code, I supply the information via a txt file. 
## Note: At the end of the code, there is a mutate function being called, which is 
### specific to the system you are running. Either modify it or comment out if you 
### don't need it. 
soloParSim <- function(topoFile, Params, ParID, 
                            noiseLevel = 0.1, scaledNoise = F,
                            deltaT = 0.1, tMax = 200, nInit = 100, parallel = T, workers = 1)
{
    print("reached soloParSim")
    net <- topoFile %>% str_remove(".topo")
    # browser()
    if(!dir.exists(net))
    {
        dir.create(net)
    }
    ParID <- paste0(net, "/", ParID)
    df <- read.delim(topoFile, sep = "")
    nodes <- unique(c(df$Source, df$Target)) %>% sort
    sMax <- sapply(nodes, function(node){
        1.5*Params[paste0("Prod_of_",node)]/Params[paste0("Deg_of_", node)]
    })
    names(sMax) <- nodes
    state <- sMax*runif(length(nodes))
    updFun <- emFunc(topoFile, Params, state, deltaT, noiseLevel, scaledNoise)
    if(parallel)
        plan(multiprocess, workers = workers)
    applyfunc <- ifelse(parallel, future_lapply, lapply)
    dfTraj <- applyfunc(1:nInit, function(id){
        initState <- sMax*runif(length(nodes))
        trajectory <- simulateInitCond(topoFile, Params, initState, deltaT, tMax, updFun)
        trajectory$Id <- id
        print(id)
        return(trajectory)
    })
    if(parallel)
        future:::ClusterRegistry("stop")
    dfTraj <- dfTraj %>% reduce(rbind.data.frame) %>% mutate(Id = factor(Id)) %>% 
        set_names(c(nodes, "Time", "Id"))
    meanSDFile <- read.delim("Means_and_Std_Dev.txt", sep = "", header = F)
    Means <- meanSDFile[2,] %>% unlist %>% as.numeric
    nodeNames <- topoFile %>% str_replace(".topo", ".prs") %>% read.delim %>%
        filter(str_detect(Parameter, "Prod_of")) %>% select(Parameter) %>% 
        unlist %>% str_remove("Prod_of_") 
    names(Means) <- nodeNames
    Means <- Means[nodes]
    SDs <- meanSDFile[3, ] %>% unlist %>% as.numeric
    names(SDs) <- nodeNames
    SDs <- SDs[nodes]
    
    dummy <- sapply(nodes, function(x){#browser()
        d <- dfTraj[[x]]
        d[d<0] <- 0
        d <- d + 0.1*min(d[d!=0])*runif(length(d))
        dfTraj[[x]] <<- (log2(d) - Means[x])/SDs[x]
    })
    dfTraj <- dfTraj %>% 
        set_names(c(nodes %>% str_remove("\\W+"), "Time", "Id")) %>%
        mutate(EMTScore = ZEB1-miR200, Sensitivity = -(AR+ARv7))
    
    write.csv(dfTraj, paste0(ParID, "_simDat.csv"), row.names = F)
}


# simulates the network for multiple parameters
#Inputs:
#   topoFile: the name of the topofile that contains the network. 
#           Note that a parameter set file with similar name must be available in the same folder
#   paramIndex: A vector of row numbers of the parameter file at which the network needs to be simulated
#   noiseLevel: A numeric value to be multiplied to the noise term in the equation. Default : 1.
#   scaledNoise: Boolean, whether to scale noise to the expression levels
#   detaT: time step for the simulations
#   tMax: maximum simualation time
#   nInit: number of initial conditions
#   parallel: Boolean, whether to impliment parallelization
#   workers: If parallel, then how many cores
simulateFull <- function(topoFile, paramIndex = NULL, 
                         noiseLevel = 1, scaledNoise = F,
                         deltaT = 0.1, tMax = 200, nInit = 100, parallel = T, workers = 1)
{
    print("working in simulate full")
    df <- read.delim(topoFile, sep = "")
    nodes <- unique(c(df$Source, df$Target)) %>% sort
    parameters <- paramRead(topoFile)
    # browser()
    if(is.null(paramIndex))
    {
        paramIndex <- parameters$ID
    }
    dummy <- sapply(paramIndex, function(i){
        Params <- parameters %>% filter(ID == i) %>% unlist
        soloParSim(topoFile, Params, i, 
                        noiseLevel, scaledNoise ,
                        deltaT, tMax, nInit, parallel, workers)
        print(i)
    })
    
}

# To get switching dynamics for the network for a particular parameter set specified by ParID
simSwitching <- function(topoFile, ParID, 
                       noiseLevel = 1, scaledNoise = F,
                       deltaT = 0.1, tMax = 200, nInit = 100, parallel = T, workers = 1)
{
    parameters <- paramRead(topoFile)
    Params <- parameters %>% filter(ID == ParID) %>% unlist
    net <- topoFile %>% str_remove(".topo")
    if(!dir.exists(net))
    {
        dir.create(net)
    }
    ParID <- paste0(net, "/", ParID)
    df <- read.delim(topoFile, sep = "")
    nodes <- unique(c(df$Source, df$Target)) %>% sort
    sMax <- sapply(nodes, function(node){
        1.5*Params[paste0("Prod_of_",node)]/Params[paste0("Deg_of_", node)]
    })
    names(sMax) <- nodes
    state <- sMax*runif(length(nodes))
    updFun <- emFunc(topoFile, Params, state, deltaT, noiseLevel, scaledNoise)
    meanSDFile <- read.delim("Means_and_Std_Dev.txt", header = F) #### change to Means and Std Dev
    Means <- meanSDFile[2,] %>% unlist %>% as.numeric
    names(Means) <- meanSDFile[1,] %>% unlist %>% str_remove("_.*")
    Means <- Means[nodes]
    SDs <- meanSDFile[3,] %>% unlist %>% as.numeric
    names(SDs) <- meanSDFile[1, ] %>% str_remove("_.*")
    SDs <- SDs[nodes]
    if(parallel)
        plan(multiprocess, workers = workers)
    dfTraj <- future_lapply(1:nInit, function(id){
        # browser()
        initState <- sMax*runif(length(nodes))
        iClog <- (log2(initState) - Means)/SDs
        # print(iClog)
        names(iClog) <- nodes
        # E <- iClog["ZEB1"] - iClog["MIR200"]
        # R <- iClog["ERa36"] - iClog["ERa66"]
        # while((E>1 || E < -1) || (R<0))
        # {
        #     initState <- runif(length(nodes))*sMax
        #     iClog <- (log2(initState) - Means)/SDs
        #     # print(iClog)
        #     names(iClog) <- nodes
        #     E <- iClog["ZEB1"] - iClog["MIR200"]
        #     R <- iClog["ERa36"] - iClog["ERa66"]
        # }
        trajectory <- simulateInitCond(topoFile, Params, initState, deltaT, tMax, updFun)
        trajectory$Id <- id
        print(id)
        return(trajectory)
    })
    if(parallel)
        future:::ClusterRegistry("stop")
    dfTraj <- dfTraj %>% reduce(rbind.data.frame) %>% mutate(Id = factor(Id)) %>% 
        set_names(c(nodes, "Time", "Id"))
    
    dummy <- sapply(nodes, function(x){#browser()
        d <- dfTraj[[x]]
        d[d<0] <- 0
        d <- d + 0.1*min(d[d!=0])*runif(length(d))
        dfTraj[[x]] <<- (log2(d) - Means[x])/SDs[x]
    })
    dfTraj <- dfTraj %>% 
        set_names(c(nodes %>% str_remove("\\W+"), "Time", "Id")) %>%
        mutate(EMTScore = ZEB1-miR200, Sensitivity = -(AR+ARv7))
    print(paste0(ParID, "_simDat.csv"))
    
    write.csv(dfTraj, paste0(ParID, "_simDat.csv"), row.names = F)
}


