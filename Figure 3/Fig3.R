##########################################################################

##Barplots for EMT score comparison

##Author: Maalavika Pillai

##Created: 26.08.21

##Inputs

#df = dataframe  for GSE168688
#Rows correspond to samples
#Columns correspond to EMT scores (KS and GS) and sample class (DMSO, EN48 and Res)
#Column names/order: KS, GS, samples

##########################################################################

require(dplyr)
require(rstatix)
require(ggpubr)
require(readxl)


plot <- reshape2::melt(plot, id.vars="samples")   #stacking data
names(plot)[2:3] <- c("EMT", "score")   #Rename column

#t-test on data
stat.test <- plot %>%   
  group_by(EMT) %>%
  t_test(score ~ samples)%>%
  add_significance("p")
stat.test <- stat.test %>%
  add_xy_position(x = "EMT", dodge = 0.8)

#Plot data
ggbarplot(
  plot, y = "scores", x = "samples", add = "mean_sd", 
  fill= "samples", palette = c("#00AFBB", "#E7B800", "#ABB000"),title = paste0(i," ",names(plot)[1]),
  position = position_dodge(0.8)) +
  xlab("Treatment")+
  ylab("EMT score")+
  stat_pvalue_manual(stat.test,
                     label = "p.signif" , tip.length = 0)+
  theme(text = element_text(size=14))
