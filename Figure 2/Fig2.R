##########################################################################

##Correlation plots and scatter plots

##Author: Maalavika Pillai

##Created: 26.08.21

##Inputs

##1. Correlation plot: 
#correlation matrix function
#df = dataframe  for sequencing data, genes/scores in rows, samples in columns 
#list = list of genes/scores to be plotted
#method = "pearson", "spearman"
#save = TRUE/FALSE, to save image 

##2. Scatter plot:
#df = dataframe  for sequencing data, genes/scores in rows, samples in columns 


##########################################################################


#Correlation plot (Triangular)
correlation <- function(df, list, method, save = FALSE, title){
  require(ggplot2)
  require(ggcorrplot) 
  df<- df[list,]
  df<- df[complete.cases(df),]  # Removing NAs, in case they exist
  corr <- cor(t(df) , method = method)  #Correlation matrix
  corr.p <- cor_pmat(t(df), method = method ) #p-values
  if(save == T){
    jpeg(paste0(title,"_corr_",method,".jpeg"))
    print(ggcorrplot(corr,
                     hc.order = FALSE,
                     type = "lower",
                     outline.color = "white", p.mat = as.matrix(corr.p), sig.level = 0.05, insig = "pch",
                     show.diag = TRUE, show.legend = FALSE, title = paste0(title," (n=",ncol(df),")"),
                     pch.cex = 10, tl.cex = 20) + theme(text = element_text(size=20)))
    dev.off()
  }else{
    return(ggcorrplot(corr,
                      hc.order = FALSE,
                      type = "lower",
                      outline.color = "white", p.mat = as.matrix(corr.p), sig.level = 0.05, insig = "pch",
                      show.diag = TRUE, show.legend = FALSE, title = paste0(title," (n=",ncol(df),")"),
                      pch.cex = 10, tl.cex = 20) + theme(text = element_text(size=20)))
  }
}


#Scatter plot: ZEB1 v/s WANG PCA Androgen Independent
#This can be changed according to x and y-axis
 
ggscatter(as.data.frame(t(df)), x ="ZEB1", y="WANG PCA Androgen Independent
", add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman")
