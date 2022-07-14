library(tidyverse)
options(stringsAsFactors = F)
library(plotly)
source("D:/Abheepsa_and_rashi/Stochastic Simulations v3/plot_theme.R")

landScape <- function(dfTraj, nam, pThresh = -2.5)
{
    kd <- with(dfTraj, MASS::kde2d(EMTScore,Sensitivity,n = 200))
    names(kd) <- c("EMTScore", "Sensitivity", "Probability")
    kd$Probability <- log10(kd$Probability)
    kd$Probability[kd$Probability < pThresh] <- NA
    colorAxis <- kd$Probability
    colorAxis[kd$EMTScore < 0.3, ] <- -1
    colorAxis[kd$EMTScore > 0.3, ] <- 1
    colorAxis[kd$EMTScore == 0.3] <- 0
    colorAxis <- apply(colorAxis, 2, function(x){
        as.character(x)
    })
    
    colorAxis2 <- kd$Probability
    colorAxis2[kd$Sensitivity >= -0.5, ] <- -1
    colorAxis2[kd$Sensitivity < -0.5, ] <- 1
    colorAxis2 <- apply(colorAxis2, 2, function(x){
        as.character(x)
    })
    p <- plot_ly(x = kd$EMTScore, y = kd$Sensitivity, z = -kd$Probability,
                     colors = c("grey", "orange", "red")) %>% 
            add_surface(surfacecolor = colorAxis, 
                        contours = list(
                            z = list(
                                show=TRUE,
                                usecolormap=FALSE,
                                highlightcolor="#ff0000",
                                project=list(z=TRUE)
                            )
                        )) %>%
            layout(title = nam, xaxis = list(title = 'EMT score'),yaxis = list(title = 'Sensitivity score'),scene = list(camera=list(eye = list(x=1.1, y=-1, z=0.64))) )
        
    #orca(p, paste0(nam, "_EMTCol.png"))#, width = 1280, height = 720
    #orca(p, width = 1280, height = 720, format = "png", scale = 1,out_file = paste0(nam, "_EMTCol.png"))
    #plotly_IMAGE(p, format = "png", out_file =paste0(nam, "_EMTCol.png"))
    htmlwidgets::saveWidget(as_widget(p), paste0(nam, "_EMTCol.html"))
}

# Generated trajectory plots for multiple initial conditions
trajectoriesnIC <- function(dfTraj, nam,folder)
{
    dfInit <- dfTraj %>% filter(Time == max(dfTraj$Time))
    dfInit$State <- "Epithelial"
    dfInit$State[dfInit$EMTScore > 0.3] <- "Mesenchymal"
    dfInit$State[dfInit$EMTScore == 0.3 ] <- "Hybrid"
    dfTraj <- merge(dfTraj, 
                    dfInit %>% select(Id, State), by = "Id",all = T)
    dfNew <- dfTraj %>% select(Id, Time, EMTScore, Sensitivity, State) %>% 
        gather(key = "Metric", value = "Score", - Time, -Id, -State)
    dfNew$State <- factor(dfNew$State, levels = c("Epithelial", "Hybrid", "Mesenchymal"))
    p <- ggplot()
    for (i in dfNew$Id %>% unique)
    {
        p <- p + geom_line(data = dfNew %>% filter(Id == i), mapping = aes(x = Time, y = Score, color = State))
    }
    p <- p + facet_grid(rows = vars(Metric)) + 
        theme_Publication()  + 
        scale_color_manual(values= c("red", "green", "blue")) +
        theme(legend.position = "top", panel.grid = element_blank())
            
    ggsave(plot = p, filename = paste0(folder,nam, "_trajecs.png"))
}

trajPlots <- function(topoFile, workers)
{
    net <- topoFile %>% str_remove(".topo")
    setwd(net)
    filz <- list.files(".",".csv")
    plan(multiprocess, workers = workers)
    dummy <- future_lapply(filz, function(x){
        nam <- x %>% str_remove("_simDat.csv")
        dfTraj <- read.csv(x)
        trajectoriesnIC(dfTraj, nam)
    })
    future:::ClusterRegistry("stop")
    setwd("..")
    
}

#Generates trajectory plots for switching.
trajectoriesSwitching <- function(dfTraj, id, tMin = 0, tMax = NULL,folder)
{
    d <- dfTraj %>% filter(Id == id) %>% select(Time, EMTScore, Sensitivity) %>% 
        gather(key = "Metric", value = "Score", -Time)
    p <- ggplot(d, aes(x = Time, y = Score, color = Metric)) + geom_line() + theme_Publication() + 
        theme(legend.position = "top")
    ggsave(paste0(folder,id,"_Normal.png"), width = 6.5, height = 4.5)
    print(p)
    if (!is.null(tMax))
    { print("saving zoomed image")
        p <- ggplot(d %>% filter(Time >=tMin, Time<=tMax), aes(x = Time, y = Score, color = Metric)) + 
            geom_line() + theme_Publication() + 
            theme(legend.position = "top")
        ggsave(paste0(folder,id, "_Zoom.png"), width = 6.5, height = 4.5)
        print(p)
    }
}


## Generates trajectories for all files in folder:
trajectories_full_folder<- function(folder)
{
    filenames <- list.files(folder, pattern="*.csv", full.names=FALSE)
    for (file in filenames){
        file
        id=str_remove(file, "_simDat.csv")
        print(paste(folder,file,sep=""))
        dfTraj <- read.csv(paste(folder,file,sep=""))
        trajectoriesnIC(dfTraj, id,folder)
        
    }
    
}

## Generates landscape for all files in folder:
landscape_full_folder<- function(folder,filenames_list=NULL )
{
    
    print("Making lanscape for all files in a folder")
    filenames <- list.files(folder, pattern="*.csv", full.names=FALSE)
    # browser()
    if(is.null(filenames))
    {
        filenames_list <- filenames
    }
    
    dummy <- sapply(filenames_list, function(i){
        
        print(paste(folder,i,"_simDat.csv", sep=""))
        dfTraj <- read.csv(paste(folder,i,"_simDat.csv", sep=""))
        landScape(dfTraj, i)
    })

}
