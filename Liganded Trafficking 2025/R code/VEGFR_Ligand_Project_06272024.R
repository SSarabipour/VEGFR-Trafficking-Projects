##########################################################################################################################################
######### R code for Main MS and Supplementary file FIGUREs visualized for:
########## Sarabipour et al. Impact of ligand binding on VEGFR1, VEGFR2, and NRP1 localization in human endothelial cells. bioRxiv (2024)
##########################################################################################################################################

##############
#Figure 2C
##############

#whole_cell_VEGF165a_barplots

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig2C_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Receptors <- factor(df1$Receptors, levels=c("0min NT VEGF Exp 1",
                                                "15min NT VEGF Exp 1",
                                                "30min NT VEGF Exp 1",
                                                "60min NT VEGF Exp 1",
                                                "120min NT VEGF Exp 1",
                                                "240min NT VEGF Exp 1",
                                                "0min NT VEGF Exp 2",
                                                "15min NT VEGF Exp 2",
                                                "30min NT VEGF Exp 2",
                                                "60min NT VEGF Exp 2",
                                                "120min NT VEGF Exp 2",
                                                "240min NT VEGF Exp 2",
                                                "0min NT VEGF Exp 3",
                                                "15min NT VEGF Exp 3",
                                                "30min NT VEGF Exp 3",
                                                "60min NT VEGF Exp 3",
                                                "120min NT VEGF Exp 3",
                                                "240min NT VEGF Exp 3"))

bp <- ggplot(df1, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd, color = Receptors),width=0.5, alpha=2, size=0.5)+   
  
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 220), breaks = seq(0,220, by = 20),expand = c(0, 0))+labs(y = "Whole cell receptors\n (normalized to Tubulin)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#CCCCCC",
                                "#0066CC",  
                                "#0066CC", 
                                "#0066CC",
                                "#0066CC",
                                "#0066CC", 
                                
                                "#CCCCCC",
                                "#FF6666", 
                                "#FF6666", 
                                "#FF6666",
                                "#FF6666",
                                "#FF6666",
                                
                                "#CCCCCC",
                                "#339966", 
                                "#339966", 
                                "#339966",
                                "#339966",
                                "#339966"))+
  
  scale_fill_manual(values = c("#CCCCCC",
                               "#0066CC",  
                               "#0066CC", 
                               "#0066CC",
                               "#0066CC",
                               "#0066CC", 
                               
                               "#CCCCCC",
                               "#FF6666", 
                               "#FF6666", 
                               "#FF6666",
                               "#FF6666",
                               "#FF6666",
                               
                               "#CCCCCC",
                               "#339966", 
                               "#339966", 
                               "#339966",
                               "#339966",
                               "#339966"))

bp+scale_x_discrete(expand = c(0.04, 0.04))+ theme(   panel.background = element_blank(),
                                                      panel.grid.major.x = element_blank(),
                                                      panel.grid.major.y = element_blank(),
                                                      panel.grid.minor.x = element_blank(),
                                                      panel.grid.minor.y = element_blank(),
                                                      legend.title = element_blank(),
                                                      legend.text  = element_blank(),
                                                      legend.position = "none",
                                                      plot.title = element_blank(),
                                                      axis.title.x = element_text(size = 6, color = "black"),
                                                      axis.title.y = element_text(size = 7, color = "black"),
                                                      axis.ticks.x = element_blank(),
                                                      axis.ticks.length = unit(0.7, "mm"),
                                                      # width of tick marks in mm
                                                      axis.ticks.y = element_line(size = .2, colour = "black"),
                                                      axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                      axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 2C.png", height = 1.9, width = 3.8, dpi=400)

##########################################################################################################################################
##############
#Figure 2D
##############

#whole_cell_normalization_PLGF1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig2D_data.csv", sep=",", header=T)
df2<- data.frame(df)

head(df2)

df2$Receptors <- factor(df2$Receptors, levels=c("0min NT PLGF Exp 1",
                                                "15min NT PLGF Exp 1",
                                                "30min NT PLGF Exp 1",
                                                "60min NT PLGF Exp 1",
                                                "120min NT PLGF Exp 1",
                                                "240min NT PLGF Exp 1",
                                                "0min NT PLGF Exp 2",
                                                "15min NT PLGF Exp 2",
                                                "30min NT PLGF Exp 2",
                                                "60min NT PLGF Exp 2",
                                                "120min NT PLGF Exp 2",
                                                "240min NT PLGF Exp 2",
                                                "0min NT PLGF Exp 3",
                                                "15min NT PLGF Exp 3",
                                                "30min NT PLGF Exp 3",
                                                "60min NT PLGF Exp 3",
                                                "120min NT PLGF Exp 3",
                                                "240min NT PLGF Exp 3"))

bp2 <- ggplot(df2, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd, color = Receptors),width=0.5, alpha=2, size=0.5)+   
  
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 202), breaks = seq(0,202, by = 20),expand = c(0, 0))+labs(y = "Whole cell receptors\n (normalized to Tubulin)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#CCCCCC",
                                "#6699CC",  
                                "#6699CC", 
                                "#6699CC",
                                "#6699CC",
                                "#6699CC", 
                                
                                "#CCCCCC",
                                "#CC6666", 
                                "#CC6666", 
                                "#CC6666",
                                "#CC6666",
                                "#CC6666",
                                
                                "#CCCCCC",
                                "#336666", 
                                "#336666", 
                                "#336666",
                                "#336666",
                                "#336666"))+
  
  scale_fill_manual(values = c("#CCCCCC",
                               "#6699CC",  
                               "#6699CC", 
                               "#6699CC",
                               "#6699CC",
                               "#6699CC", 
                               
                               "#CCCCCC",
                               "#CC6666", 
                               "#CC6666", 
                               "#CC6666",
                               "#CC6666",
                               "#CC6666",
                               
                               "#CCCCCC",
                               "#336666", 
                               "#336666", 
                               "#336666",
                               "#336666",
                               "#336666"))

bp2+scale_x_discrete(expand = c(0.04, 0.04))+ theme( panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 2D.png", height = 2.5, width = 3, dpi=400)

##########################################################################################################################################
#############
#Figure 2G
############

#Surf_Int_normalization_VEGF165a_PLGF_1h

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig2G_data.csv", sep=",", header=T)
df3<- data.frame(df)

head(df3)

df3$Receptors <- factor(df3$Receptors, levels=c("1hr NT Exp 1s",
                                                "1hr NT VEGF Exp 1s",
                                                "1hr NT PLGF Exp 1s",
                                                "1hr NT Exp 1int",
                                                "1hr NT VEGF Exp 1int",
                                                "1hr NT PLGF Exp 1int",
                                                "1hr NT Exp 1tot",
                                                "1hr NT VEGF Exp 1tot",
                                                "1hr NT PLGF Exp 1tot",
                                                "1hr NT Exp 2s",
                                                "1hr NT VEGF Exp 2s",
                                                "1hr NT PLGF Exp 2s",
                                                "1hr NT Exp 2int",
                                                "1hr NT VEGF Exp 2int",
                                                "1hr NT PLGF Exp 2int",
                                                "1hr NT Exp 2tot",
                                                "1hr NT VEGF Exp 2tot",
                                                "1hr NT PLGF Exp 2tot",
                                                "1hr NT Exp 3s",
                                                "1hr NT VEGF Exp 3s",
                                                "1hr NT PLGF Exp 3s",
                                                "1hr NT Exp 3int",
                                                "1hr NT VEGF Exp 3int",
                                                "1hr NT PLGF Exp 3int",
                                                "1hr NT Exp 3tot",
                                                "1hr NT VEGF Exp 3tot",
                                                "1hr NT PLGF Exp 3tot"))

bp2 <- ggplot(df3, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd, color = Receptors),width=0.5, alpha=2, size=0.5)+   
  
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 160), breaks = seq(0,160, by = 20),expand = c(0, 0))+labs(y = "Receptors\n (normalized to Tubulin)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#CCCCCC",
                                "#0066CC",
                                "#6699CC",
                                "#CCCCCC",
                                "#0066CC",
                                "#6699CC",
                                "#CCCCCC",
                                "#0066CC",
                                "#6699CC",
                                
                                "#CCCCCC",
                                "#FF6666", 
                                "#CC6666",
                                "#CCCCCC",
                                "#FF6666", 
                                "#CC6666",
                                "#CCCCCC",
                                "#FF6666", 
                                "#CC6666",
                                
                                "#CCCCCC",
                                "#339966",
                                "#336666",
                                "#CCCCCC",
                                "#339966",
                                "#336666",
                                "#CCCCCC",
                                "#339966",
                                "#336666"))+
  
  scale_fill_manual(values = c("#CCCCCC",
                               "#0066CC",
                               "#6699CC",
                               "#CCCCCC",
                               "#0066CC",
                               "#6699CC",
                               "#CCCCCC",
                               "#0066CC",
                               "#6699CC",
                               
                               "#CCCCCC",
                               "#FF6666", 
                               "#CC6666",
                               "#CCCCCC",
                               "#FF6666", 
                               "#CC6666",
                               "#CCCCCC",
                               "#FF6666", 
                               "#CC6666",
                               
                               "#CCCCCC",
                               "#339966",
                               "#336666",
                               "#CCCCCC",
                               "#339966",
                               "#336666",
                               "#CCCCCC",
                               "#339966",
                               "#336666"))

bp2+scale_x_discrete(expand = c(0.03, 0.03))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 2G.png", height = 2.5, width = 3.5, dpi=400)

##########################################################################################################################################
###########
#Figure 2H
###########

#Surf_Int_normalization_VEGF165a_PLGF_4h

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig2H_data.csv", sep=",", header=T)
df4<- data.frame(df)

head(df4)

df4$Receptors <- factor(df4$Receptors, levels=c("4hr NT Exp 1s",
                                                "4hr NT VEGF Exp 1s",
                                                "4hr NT PLGF Exp 1s",
                                                "4hr NT Exp 1int",
                                                "4hr NT VEGF Exp 1int",
                                                "4hr NT PLGF Exp 1int",
                                                "4hr NT Exp 1tot",
                                                "4hr NT VEGF Exp 1tot",
                                                "4hr NT PLGF Exp 1tot",
                                                "4hr NT Exp 2s",
                                                "4hr NT VEGF Exp 2s",
                                                "4hr NT PLGF Exp 2s",
                                                "4hr NT Exp 2int",
                                                "4hr NT VEGF Exp 2int",
                                                "4hr NT PLGF Exp 2int",
                                                "4hr NT Exp 2tot",
                                                "4hr NT VEGF Exp 2tot",
                                                "4hr NT PLGF Exp 2tot",
                                                "4hr NT Exp 3s",
                                                "4hr NT VEGF Exp 3s",
                                                "4hr NT PLGF Exp 3s",
                                                "4hr NT Exp 3int",
                                                "4hr NT VEGF Exp 3int",
                                                "4hr NT PLGF Exp 3int",
                                                "4hr NT Exp 3tot",
                                                "4hr NT VEGF Exp 3tot",
                                                "4hr NT PLGF Exp 3tot"))


bp2 <- ggplot(df4, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd, color = Receptors),width=0.5, alpha=2, size=0.5)+   
  
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 220), breaks = seq(0,220, by = 20),expand = c(0, 0))+labs(y = "Receptors\n (normalized to Tubulin)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#CCCCCC",
                                "#0066CC",
                                "#6699CC",
                                "#CCCCCC",
                                "#0066CC",
                                "#6699CC",
                                "#CCCCCC",
                                "#0066CC",
                                "#6699CC",
                                
                                "#CCCCCC",
                                "#FF6666", 
                                "#CC6666",
                                "#CCCCCC",
                                "#FF6666", 
                                "#CC6666",
                                "#CCCCCC",
                                "#FF6666", 
                                "#CC6666",
                                
                                "#CCCCCC",
                                "#339966",
                                "#336666",
                                "#CCCCCC",
                                "#339966",
                                "#336666",
                                "#CCCCCC",
                                "#339966",
                                "#336666"))+
  
  
  scale_fill_manual(values = c("#CCCCCC",
                               "#0066CC",
                               "#6699CC",
                               "#CCCCCC",
                               "#0066CC",
                               "#6699CC",
                               "#CCCCCC",
                               "#0066CC",
                               "#6699CC",
                               
                               "#CCCCCC",
                               "#FF6666", 
                               "#CC6666",
                               "#CCCCCC",
                               "#FF6666", 
                               "#CC6666",
                               "#CCCCCC",
                               "#FF6666", 
                               "#CC6666",
                               
                               "#CCCCCC",
                               "#339966",
                               "#336666",
                               "#CCCCCC",
                               "#339966",
                               "#336666",
                               "#CCCCCC",
                               "#339966",
                               "#336666"))


bp2+scale_x_discrete(expand = c(0.03, 0.03))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 2H.png", height = 1.5, width = 5, dpi=400)

##########################################################################################################################################
#######################
#Figure 3A
#######################

#Surface R2 receptors fit to kint variation
#165a_Surf_perc_multi_kint_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3A_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333","#CC3333", "#CC3333",  "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 130), breaks = seq(0, 130, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure3A.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure 3B
#######################

#Internal R2 receptors fit to kint variation
#165a_Internal_perc_multi_kint_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3B_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333","#CC3333", "#CC3333",  "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Internal Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure 3B.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure 3C
#######################

#Whole cell R2 receptors fit to kint variation
#165a_Tot_perc_multi_kint_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3C_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333","#CC3333", "#CC3333",  "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Whole cell Receptors (% of control)",limits = c(0, 101), breaks = seq(0, 101, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure 3C.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure 3D
#######################

#Surface R2 receptors fit to k_deg variation
#165a_Surf_perc_multi_kDeg_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3D_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333","#CC3333", "#CC3333",  "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 101), breaks = seq(0, 101, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure 3D.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure 3E
#######################

#Internal R2 receptors fit to k_deg variation
#165a_Internal_perc_multi_kDeg_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3E_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333","#CC3333", "#CC3333",  "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Internal Receptors (% of control)",limits = c(0, 180), breaks = seq(0, 180, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure 3E.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure 3F
#######################

#Whole cell R2 receptors fit to k_deg variation
#165a_Tot_perc_multi_kDeg_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3F_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333","#CC3333", "#CC3333",  "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Whole cell Receptors (% of control)",limits = c(0, 120), breaks = seq(0, 120, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure 3F.png", height = 3, width = 3.5, dpi=400)

##############################
########Figure 3G
#############################

#Percent Surface VEGFR1
#165a_Surf_R1_perc_re-Opt1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3G_data.csv", sep=",", header=T)
df200<- data.frame(df)

p = ggplot(data=df200, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-20),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#0066CC"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 130), breaks = seq(0, 130, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure 3G.png", height = 3, width = 3.5, dpi=400)

##############################
########Figure 3H
#############################

#Percent Surface VEGFR2
#165a_Surf_R2_perc_re-Opt1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3H_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-20),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#990000"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 102), breaks = seq(0, 102, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure 3H.png", height = 3, width = 3.5, dpi=400)

##############################
########Figure 3I
#############################

#Percent Surface Neuropilin-1
#165a_Surf_N1_perc_re-Opt1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3I_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-20),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#009966"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 130), breaks = seq(0, 130, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure 3I.png", height = 3, width = 3.5, dpi=400)

################################
#Figure 3J 
################################

#R1_Surf_Int_ratio_barplots

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3J_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("1hr NT H2O Exp", 
                                                "1hr NT H2O Sim",
                                                "1hr NT VEGF Exp", 
                                                "1hr NT VEGF Sim",
                                                "1hr NT PLGF Exp", 
                                                "1hr NT PLGF Sim",
                                                "4hr NT H2O Exp", 
                                                "4hr NT H2O Sim",
                                                "4hr NT VEGF Exp", 
                                                "4hr NT VEGF Sim",
                                                "4hr NT PLGF Exp", 
                                                "4hr NT PLGF Sim"))


bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd, color = Receptors),width=0.5, alpha=2, size=0.5)+   
  
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 160), breaks = seq(0,160, by = 20),expand = c(0, 0))+labs(y = "Surface to Internal ratio\n (normalized to Tubulin)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#CCCCCC","#999999",
                                "#0066CC", "#003366",
                                "#6699CC", "#336699",
                                
                                "#CCCCCC","#999999",
                                "#0066CC", "#003366",
                                "#6699CC", "#336699"))+
  
  scale_fill_manual(values = c("#CCCCCC","#999999",
                               "#0066CC", "#003366",
                               "#6699CC", "#336699",
                               
                               "#CCCCCC","#999999",
                               "#0066CC", "#003366",
                               "#6699CC", "#336699"))


bp2+scale_x_discrete(expand = c(0.05, 0.05))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 3J.png", height = 2.5, width = 3, dpi=400)

################################
#Figure 3K
################################

#R2_Surf_Int_ratio_barplots

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3K_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("1hr NT H2O Exp", 
                                                "1hr NT H2O Sim",
                                                "1hr NT VEGF Exp", 
                                                "1hr NT VEGF Sim",
                                                "1hr NT PLGF Exp", 
                                                "1hr NT PLGF Sim",
                                                "4hr NT H2O Exp", 
                                                "4hr NT H2O Sim",
                                                "4hr NT VEGF Exp", 
                                                "4hr NT VEGF Sim",
                                                "4hr NT PLGF Exp", 
                                                "4hr NT PLGF Sim"))


bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd, color = Receptors),width=0.5, alpha=2, size=0.5)+   
  
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 140), breaks = seq(0,140, by = 20),expand = c(0, 0))+labs(y = "Surface to Internal ratio\n (normalized to Tubulin)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#CCCCCC","#999999",
                                "#FF6666", "#990000",
                                "#CC6666", "#330000",
                                
                                "#CCCCCC","#999999",
                                "#FF6666", "#990000",
                                "#CC6666", "#330000"))+
  
  scale_fill_manual(values = c("#CCCCCC","#999999",
                               "#FF6666", "#990000",
                               "#CC6666", "#330000",
                               
                               "#CCCCCC","#999999",
                               "#FF6666", "#990000",
                               "#CC6666", "#330000"))

bp2+scale_x_discrete(expand = c(0.05, 0.05))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 3K.png", height = 2.5, width = 3, dpi=400)

################################
#Figure 3L
################################

#N1_Surf_Int_ratio_barplots

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3L_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("1hr NT H2O Exp", 
                                                "1hr NT H2O Sim",
                                                "1hr NT VEGF Exp", 
                                                "1hr NT VEGF Sim",
                                                "1hr NT PLGF Exp", 
                                                "1hr NT PLGF Sim",
                                                "4hr NT H2O Exp", 
                                                "4hr NT H2O Sim",
                                                "4hr NT VEGF Exp", 
                                                "4hr NT VEGF Sim",
                                                "4hr NT PLGF Exp", 
                                                "4hr NT PLGF Sim"))


bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd, color = Receptors),width=0.5, alpha=2, size=0.5)+   
  
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0,160), breaks = seq(0,160, by = 20),expand = c(0, 0))+labs(y = "Surface to Internal ratio\n (normalized to Tubulin)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#CCCCCC","#999999",
                                "#339966", "#006600",
                                "#336666", "#003333",
                                
                                "#CCCCCC","#999999",
                                "#339966", "#006600",
                                "#336666", "#003333"))+
  
  scale_fill_manual(values = c("#CCCCCC","#999999",
                               "#339966", "#006600",
                               "#336666", "#003333",
                               
                               "#CCCCCC","#999999",
                               "#339966", "#006600",
                               "#336666", "#003333"))

bp2+scale_x_discrete(expand = c(0.05, 0.05))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 3L.png", height = 2.5, width = 3, dpi=400)

###########################################
############## Figure 5 A
#########################################

#all_lig_Surf_Internal_1h

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig5A_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J",
                                        "R2 Sim K",
                                        "R2 Sim L"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+  
  
  #geom_errorbar(aes(ymin=Rvalue-RSD, ymax=Rvalue+RSD, color = Rcomp),width=0.5, alpha=2, size=0.5)+    
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 16500), breaks = seq(0,16500, by = 2000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#FF9999",
                                "#990000",
                                "#FF9999",
                                "#990000",
                                "#99CCFF",
                                "#003366",
                                "#FF9999",
                                "#990000",
                                "#FF9999",
                                "#990000",
                                "#99CCFF",
                                "#003366" ))+
  
  scale_fill_manual(values = c("#FF9999",
                               "#990000",
                               "#FF9999",
                               "#990000",
                               "#99CCFF",
                               "#003366",
                               "#FF9999",
                               "#990000",
                               "#FF9999",
                               "#990000",
                               "#99CCFF",
                               "#003366"))

bp+scale_x_discrete(expand = c(0.05, 0.05))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 5A.png", height = 2, width = 3, dpi=400)

###########################################
############## Figure 5 B
#########################################

#all_lig_Surf_Internal_4h

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig5B_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J",
                                        "R2 Sim K",
                                        "R2 Sim L"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+  
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 16500), breaks = seq(0,16500, by = 2000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#FF9999",
                                "#990000",
                                "#FF9999",
                                "#990000",
                                "#99CCFF",
                                "#003366",
                                "#FF9999",
                                "#990000",
                                "#FF9999",
                                "#990000",
                                "#99CCFF",
                                "#003366" ))+
  
  scale_fill_manual(values = c("#FF9999",
                               "#990000",
                               "#FF9999",
                               "#990000",
                               "#99CCFF",
                               "#003366",
                               "#FF9999",
                               "#990000",
                               "#FF9999",
                               "#990000",
                               "#99CCFF",
                               "#003366"))

bp+scale_x_discrete(expand = c(0.05, 0.05))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 5B.png", height = 2, width = 3, dpi=400)

###########################################
############## Figure 5 C
#########################################

#all_lig_1h

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig5C_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.85)+  
  theme_classic()+theme(aspect.ratio = .4)+xlab("")+ 
  scale_y_continuous(limits = c(0, 1.8), breaks = seq(0,1.8, by = 0.2),expand = c(0.005, 0.005))+
  labs(y = "Ligand Concentration (nM)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF",
                                "#003366",
                                "#FF9999",
                                "#990000",
                                "#99CCFF",
                                "#003366",
                                "#FF9999",
                                "#990000"))+
  
  scale_fill_manual(values = c("#99CCFF",
                               "#003366",
                               "#FF9999",
                               "#990000",
                               "#99CCFF",
                               "#003366",
                               "#FF9999",
                               "#990000"))

bp+scale_x_discrete(expand = c(0.07, 0.07))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 5C.png", height = 2, width = 3, dpi=400)

###########################################
############## Figure 5 D
#########################################

#all_lig_4h

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig5D_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.85)+  
  theme_classic()+theme(aspect.ratio = .4)+xlab("")+ 
  scale_y_continuous(limits = c(0, 1.8), breaks = seq(0,1.8, by = 0.2),expand = c(0.005, 0.005))+
  labs(y = "Ligand Concentration (nM)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF",
                                "#003366",
                                "#FF9999",
                                "#990000",
                                "#99CCFF",
                                "#003366",
                                "#FF9999",
                                "#990000"))+
  
  scale_fill_manual(values = c("#99CCFF",
                               "#003366",
                               "#FF9999",
                               "#990000",
                               "#99CCFF",
                               "#003366",
                               "#FF9999",
                               "#990000"))

bp+scale_x_discrete(expand = c(0.07, 0.07))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure 5D.png", height = 2, width = 3, dpi=400)

######################################################################################
########################### Supplementary File ################
######################################################################################

#######################
#Figure S2 A
#######################

#Surface R2 receptors fit to kint variation
#165a_Surf_perc_multi_kint_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2A_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 130), breaks = seq(0, 130, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2A.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 B
#######################

#Internal R2 receptors fit to kint variation
#165a_Internal_perc_multi_kint_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2B_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2B.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 C
#######################

#Whole cell R2 receptors fit to kint variation
#165a_Tot_perc_multi_kint_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2C_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 120), breaks = seq(0, 120, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2C.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 D
#######################

#Surface R2 receptors fit to k_deg variation
#165a_Surf_perc_multi_kDeg_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2D_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 100), breaks = seq(0, 100, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2D.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 E
#######################

#Internal R2 receptors fit to k_deg variation
#165a_Internal_perc_multi_kDeg_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2E_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dotted","dashed","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 200), breaks = seq(0, 200, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2E.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 F
#######################

#Whole cell R2 receptors fit to k_deg variation
#165a_Tot_perc_multi_kDeg_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2F_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dotted","dashed","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 120), breaks = seq(0, 120, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2F.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 G
#######################

#Surface R2 receptors fit to k4a variation
#165a_Surf_perc_multi_k4tos_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2G_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 103), breaks = seq(0, 103, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure SG.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 H
#######################

#Internal R2 receptors fit to k4tos variation
#165a_Internal_perc_multi_k4tos_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2H_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 140), breaks = seq(0, 140, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2H.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 I
#######################

#Whole cell R2 receptors fit to kint variation
#165a_Tot_perc_multi_k4tos_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2I_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 103), breaks = seq(0, 103, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2H.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 J
#######################

#Surface R2 receptors fit to k4to11 variation
#165a_Surf_perc_multi_k4to11_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2J_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 103), breaks = seq(0, 103, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2J.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 K
#######################

#Internal R2 receptors fit to k4to11 variation
#165a_Internal_perc_multi_k4to11_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2K_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 140), breaks = seq(0, 140, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2K.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 L
#######################

#Whole cell R2 receptors fit to k4to11 variation
#165a_Tot_perc_multi_k4to11_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2L_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 103), breaks = seq(0, 103, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2L.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 M
#######################

#Surface R2 receptors fit to k11tos variation
#165a_Surf_perc_multi_k11tos_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2M_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 103), breaks = seq(0, 103, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2M.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 N
#######################

#Internal R2 receptors fit to k11tos variation
#165a_Internal_perc_multi_k11tos_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2N_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 140), breaks = seq(0, 140, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2N.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S2 O
#######################

#Whole cell R2 receptors fit to k11tos variation
#165a_Tot_perc_multi_k11tos_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS2O_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","solid","dashed","longdash","dotted"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#CC3333","#CC3333", "#CC3333", "#CC3333", "#CC3333"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 103), breaks = seq(0, 103, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S2O.png", height = 3, width = 3.5, dpi=400)

######################################
####### Figure S3 A
######################################

#Surface R1 receptors fit to kint variation
#165a_Surf_perc_multi_kint_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3A_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 210), breaks = seq(0, 210, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3A.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 B
#######################

#Internal R1 receptors fit to kint variation
#165a_Internal_perc_multi_kint_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3B_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Internal Receptors (% of control)",limits = c(0, 125), breaks = seq(0, 125, by = 30),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3B.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 C
#######################

#Whole cell R1 receptors fit to kint variation
#165a_Tot_perc_multi_kint_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3C_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Whole Cell Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3C.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 D
#######################

#Surface R1 receptors fit to k_deg variation
#165a_Surf_perc_multi_kDeg_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3D_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 190), breaks = seq(0, 190, by = 30),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3D.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 E
#######################

#Internal R1 receptors fit to k_deg variation
#165a_Internal_perc_multi_kDeg_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3E_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Internal Receptors (% of control)",limits = c(0, 210), breaks = seq(0, 210, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3E.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 F
#######################

#Whole cell R1 receptors fit to k_deg variation
#165a_Tot_perc_multi_kDeg_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3F_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Whole Cell Receptors (% of control)",limits = c(0, 210), breaks = seq(0, 210, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3F.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 G
#######################

#Surface R1 receptors fit to k4a variation
#165a_Surf_perc_multi_k4tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3G_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 200), breaks = seq(0, 200, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3G.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 H
#######################

#Internal R1 receptors fit to k4tos variation
#165a_Internal_perc_multi_k4tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3H_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Internal Receptors (% of control)",limits = c(0, 125), breaks = seq(0, 125, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3H.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 I
#######################

#Whole cell R1 receptors fit to kint variation
#165a_Tot_perc_multi_k4tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3I_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Whole cell Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3I.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 J
#######################

#Surface R1 receptors fit to k4to11 variation
#165a_Surf_perc_multi_k4to11_R1

library(ggplot2)
library(dplyr)
library(Cairo)
library(reshape2)

df<-read.csv("FigS3J_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 200), breaks = seq(0, 200, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3J.png",  height = 3, width = 3.5, type = "cairo", dpi=400)

#######################
#Figure S3 K
#######################

#Internal R1 receptors fit to k4to11 variation
#165a_Internal_perc_multi_k4to11_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3K_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Internal Receptors (% of control)",limits = c(0, 125), breaks = seq(0, 125, by = 25),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3K.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 L
#######################

#Whole cell R1 receptors fit to k4to11 variation
#165a_Tot_perc_multi_k4to11_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3L_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Whole cell Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3L.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 M
#######################

#Surface R1 receptors fit to k11tos variation
#165a_Surf_perc_multi_k11tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3M_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 130), breaks = seq(0, 130, by = 25),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3M.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 N
#######################

#Internal R1 receptors fit to k11tos variation
#165a_Internal_perc_multi_k11tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3N_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Internal Receptors (% of control)",limits = c(0, 125), breaks = seq(0, 125, by = 25),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3N.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S3 O
#######################

#Whole cell R1 receptors fit to k11tos variation
#165a_Tot_perc_multi_k11tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3O_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#3399FF", "#3399FF","#3399FF","#3399FF","#3399FF"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Whole cell Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S3O.png", height = 3, width = 3.5, dpi=400)

########################################################################
#Figure S4 A
#######################

#Surface R1 receptors fit to kint variation
#PLGF1_Surf_perc_multi_kint_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4A_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 210), breaks = seq(0, 210, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4A.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 B
#######################

#Internal R1 receptors fit to kint variation
#PLGF1_Internal_perc_multi_kint_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4B_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Internal Receptors (% of control)",limits = c(0, 125), breaks = seq(0, 125, by = 30),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4B.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 C
#######################

#Whole cell R1 receptors fit to kint variation
#PLGF1_Tot_perc_multi_kint_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4C_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Whole Cell Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4C.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 D
#######################

#Surface R1 receptors fit to k_deg variation
#PLGF1_Surf_perc_multi_kDeg_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4D_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 200), breaks = seq(0, 200, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4D.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 E
#######################

#Internal R1 receptors fit to k_deg variation
#PLGF1_Internal_perc_multi_kDeg_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4E_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Internal Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4E.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 F
#######################

#Whole cell R1 receptors fit to k_deg variation
#PLGF1_Tot_perc_multi_kDeg_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4F_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("whole Cell Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4F.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 G
#######################

#Surface R1 receptors fit to k4a variation
#PLGF1_Surf_perc_multi_k4tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4G_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 200), breaks = seq(0, 200, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4G.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 H
#######################

#Internal R1 receptors fit to k4tos variation
#PLGF1_Internal_perc_multi_k4tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4H_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 125), breaks = seq(0, 125, by = 25),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4H.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 I
#######################

#Whole cell R1 receptors fit to kint variation
#PLGF1_Tot_perc_multi_k4tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4I_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4I.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 J
#######################

#Surface R1 receptors fit to k4to11 variation
#PLGF1_Surf_perc_multi_k4to11_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4J_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 200), breaks = seq(0, 200, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4J.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 K
#######################

#Internal R1 receptors fit to k4to11 variation
#PLGF1_Internal_perc_multi_k4to11_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4K_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 125), breaks = seq(0, 125, by = 25),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4K.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 L
#######################

#Whole cell R1 receptors fit to k4to11 variation
#PLGF1_Tot_perc_multi_k4to11_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4L_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4L.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 M
#######################

#Surface R1 receptors fit to k11tos variation
#PLGF1_Surf_perc_multi_k11tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4M_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 200), breaks = seq(0, 200, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4M.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 N
#######################

#Internal R1 receptors fit to k11tos variation
#PLGF1_Internal_perc_multi_k11tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4N_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 125), breaks = seq(0, 125, by = 25),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("Figure S4N.png", height = 3, width = 3.5, dpi=400)

#######################
#Figure S4 O
#######################

#Whole cell R1 receptors fit to k11tos variation
#PLGF1_Tot_perc_multi_k11tos_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4O_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("twodash","dashed","dotted","longdash","solid"))+
  
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=10, alpha=2, size=0.9,position=position_dodge(2))+
  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(axis.title.x = element_text(colour = "black", size = 9),
        axis.title.y = element_text(colour = "black", size = 9),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=7, color = "black"),
        
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#000066","#000066","#000066","#000066","#000066"))+
  labs(title = "", x = "Time (min)", y = "")+ theme(legend.title=element_blank())+
  scale_y_continuous("Surface Receptors (% of control)",limits = c(0, 160), breaks = seq(0, 160, by = 40),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 245), breaks = seq(0, 245, by = 30),expand = c(0, 0))
p

ggsave ("FigSure S4O.png", height = 3, width = 3.5, dpi=400)

#####################################################################################################################
###################################FIGURE S7 Flux plots for VEGFR1, VEGFR2, Neuropilin1##############################
#####################################################################################################################

####################################################################################
###FIGURE S7A Flux plots for VEGFR1 at Steady State (no ligand treatment)
####################################################################################

library(ggplot2)

df <- data.frame(Complex=rep(c("recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                "transfer from Rab4a to Rab7a", "production to surface", "internalization from surface to Rab4a"),each=6),
  
  Cargo_Type= rep(c("enters surface", "leaves surface", "enters Rab4a/5a",  "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"), each = 1),
  Conc = c(37.74, 0,0,37.74,0,0,       32.33,0,0,0,0,32.33,      0,0,0,32.34,32.34,0,
           0,0,0,3.76,0,0,         3.8,0,0,0,0,0,        0,73.85,73.85,0,0,0))


df$Complex <- factor(df$Complex, levels = c("transfer from Rab4a to Rab7a","production to surface","recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                            "internalization from surface to Rab4a"))

df$Cargo_Type <- factor(df$Cargo_Type, levels=c("enters surface", "leaves surface", "enters Rab4a/5a",  
                                                "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"))

head(df)

p <- ggplot(df, aes(x = Cargo_Type, y = Conc, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c("#0066CC", "#6699CC","#003366","#336699", "#99CCFF","#000033"))

p+theme_classic()+scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 20),expand = c(0.015, 0.015))+
  scale_x_discrete("",expand = c(0.08, 0.08))+labs(title="")+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=8, color = "black"),
        legend.text  = element_text(size=7, color = "black"),
        plot.title = element_text(size=12, color = "black", hjust=0.6),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        legend.position = "right") + ylab(expression(paste("Receptors per cell per second")))+
  coord_flip()

ggsave("Fig7A_R1_Flux_Steady_State.png", width=5, height=4, dpi=400)

##################################################################################
###FIGURE S7D Flux plots for VEGFR1 after 1 hour of 50 ng/mL VEGF-A165a Treatment
##################################################################################

library(ggplot2)

df <- data.frame(Complex=rep(c("recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                "transfer from Rab4a to Rab7a", "production to surface", "internalization from surface to Rab4a"),each=6),
  
  Cargo_Type = rep(c("enters surface", "leaves surface", "enters Rab4a/5a",  "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"), each = 1),
        Conc = c(9.598,0,0,9.598,0,0,       10.568,0,0,0,0,10.568,       0,0,0,10.568,10.568,0,
                 0,0,0,4.101,0,0,         4.101,0,0,0,0,0,           0,24.27,24.27,0,0,0))

df$Complex <- factor(df$Complex, levels = c("transfer from Rab4a to Rab7a","production to surface","recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                            "internalization from surface to Rab4a"))

df$Cargo_Type <- factor(df$Cargo_Type, levels=c("enters surface", "leaves surface", "enters Rab4a/5a",  
                                                "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"))
head(df)

p <- ggplot(df, aes(x = Cargo_Type, y = Conc, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c("#0066CC", "#6699CC","#003366","#336699", "#99CCFF","#000033"))

p+theme_classic()+scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5),expand = c(0.015, 0.015))+
  scale_x_discrete("",expand = c(0.08, 0.08))+labs(title="")+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=8, color = "black"),
        legend.text  = element_text(size=7, color = "black"),
        plot.title = element_text(size=12, color = "black", hjust=0.6),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        legend.position = 'none') + ylab(expression(paste("Receptors per cell per second")))+
  coord_flip()

ggsave("FigS7D_R1_Flux_1h.png", width=5, height=4, dpi=400)

####################################################################################
###FIGURE S7G Flux plots for VEGFR1 after 4 hours of 50 ng/mL VEGF-A165a Treatment
####################################################################################

library(ggplot2)

df <- data.frame(Complex=rep(c("recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                "transfer from Rab4a to Rab7a", "production to surface", "internalization from surface to Rab4a"),each=6),
  
  Cargo_Type = rep(c("enters surface", "leaves surface", "enters Rab4a/5a",  "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"), each = 1),
        Conc = c(9.598,0,0,9.598,0,0,       10.568,0,0,0,0,10.568,       0,0,0,10.568,10.568,0,
                 0,0,0,4.101,0,0,         4.101,0,0,0,0,0,           0,24.27,24.27,0,0,0))

df$Complex <- factor(df$Complex, levels = c("transfer from Rab4a to Rab7a","production to surface","recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                            "internalization from surface to Rab4a"))

df$Cargo_Type <- factor(df$Cargo_Type, levels=c("enters surface", "leaves surface", "enters Rab4a/5a",  
                                                "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"))
head(df)

p <- ggplot(df, aes(x = Cargo_Type, y = Conc, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c("#0066CC", "#6699CC","#003366","#336699", "#99CCFF","#000033"))

p+theme_classic()+scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5),expand = c(0.015, 0.015))+
  scale_x_discrete("",expand = c(0.08, 0.08))+labs(title="")+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=8, color = "black"),
        legend.text  = element_text(size=7, color = "black"),
        plot.title = element_text(size=12, color = "black", hjust=0.6),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        legend.position = 'none') + ylab(expression(paste("Receptors per cell per second")))+
  coord_flip()

ggsave("FigS7G_R1_Flux_4h.png", width=5, height=4, dpi=400)

####################################################################################
###FIGURE S7B Flux plots for VEGFR2 at Steady State (no ligand treatment)
####################################################################################

library(ggplot2)

df11 <- data.frame(Complex=rep(c("recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                "transfer from Rab4a to Rab7a", "production to surface", "internalization from surface to Rab4a"),each=6),
  
  Cargo_Type= rep(c("enters surface", "leaves surface", "enters Rab4a/5a",  "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"), each = 1),
  Conc = c(12.53, 0,0,12.53,0,0,       3.65,0,0,0,0,3.65,    0,0,0,3.65,3.65,0,
           0,0,0,1.5,0,0,         1.5,0,0,0,0,0,     0,17.7,17.7,0,0,0))

df11$Complex <- factor(df11$Complex, levels = c("transfer from Rab4a to Rab7a","production to surface","recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                                "internalization from surface to Rab4a"))

df11$Cargo_Type <- factor(df11$Cargo_Type, levels=c("enters surface", "leaves surface", "enters Rab4a/5a",  
                                                    "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"))

head(df11)

p <- ggplot(df11, aes(x = Cargo_Type, y = Conc, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c( "#993333","#660000","#663333", "#996666","#CC9999","#330000"))

p+theme_classic()+scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5),expand = c(0.015, 0.015))+
  scale_x_discrete("",expand = c(0.08, 0.08))+labs(title="")+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=8, color = "black"),
        legend.text  = element_text(size=7, color = "black"),
        plot.title = element_text(size=11, color = "black", hjust=0.6, vjust=5),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        legend.position  = "right")+ ylab(expression(paste("Receptors per cell per second")))+
  coord_flip()

ggsave("Fig7B_R2_Flux_Steady_State.png", width=5, height=4, dpi=400)

####################################################################################
###FIGURE S7E Flux plots for VEGFR2 after 1 hours of 50 ng/mL VEGF-A165a Treatment
####################################################################################
library(ggplot2)

df11 <- data.frame(Complex=rep(c("recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                "transfer from Rab4a to Rab7a", "production to surface", "internalization from surface to Rab4a"),each=6),
  
  Cargo_Type = rep(c("enters surface", "leaves surface", "enters Rab4a/5a",  "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"), each = 1),
        Conc = c(0.00739,0,0,0.00739,0,0,       0.008933,0,0,0,0,0.008933,       0,0,0,0.00925,0.00925,0,
                 0,0,0,1.451,0,0,         1.114,0,0,0,0,0,           0,1.3192,1.3192,0,0,0))

df11$Complex <- factor(df11$Complex, levels = c("transfer from Rab4a to Rab7a","production to surface","recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                                "internalization from surface to Rab4a"))

df11$Cargo_Type <- factor(df11$Cargo_Type, levels=c("enters surface", "leaves surface", "enters Rab4a/5a",  
                                                    "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"))

head(df11)

p <- ggplot(df11, aes(x = Cargo_Type, y = Conc, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c( "#993333","#660000","#663333", "#996666","#CC9999","#330000"))

p+theme_classic()+scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, by = 0.5),expand = c(0.01, 0.01))+
  scale_x_discrete("",expand = c(0.08, 0.08))+labs(title="")+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=8, color = "black"),
        legend.text  = element_text(size=7, color = "black"),
        plot.title = element_text(size=11, color = "black", hjust=0.6, vjust=5),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        legend.position  = 'none')+ ylab(expression(paste("Receptors per cell per second")))+
  coord_flip()

ggsave("FigS7E_R2_Flux_1h.png", width=5, height=4, dpi=400)

####################################################################################
###FIGURE S7H Flux plots for VEGFR2 after 4 hours of 50 ng/mL VEGF-A165a Treatment
####################################################################################
library(ggplot2)

df11 <- data.frame(Complex=rep(c("recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                 "transfer from Rab4a to Rab7a", "production to surface", "internalization from surface to Rab4a"),each=6),
                   
                   Cargo_Type = rep(c("enters surface", "leaves surface", "enters Rab4a/5a",  "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"), each = 1),
                         Conc = c(0.00586,0,0,0.00586,0,0,       0.00715,0,0,0,0,0.00715,       0,0,0,0.007338,0.007338,0,
                                  0,0,0,1.1503,0,0,         1.114,0,0,0,0,0,           0,1.1274,1.1274,0,0,0))

df11$Complex <- factor(df11$Complex, levels = c("transfer from Rab4a to Rab7a","production to surface","recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                                "internalization from surface to Rab4a"))

df11$Cargo_Type <- factor(df11$Cargo_Type, levels=c("enters surface", "leaves surface", "enters Rab4a/5a",  
                                                    "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"))

head(df11)

p <- ggplot(df11, aes(x = Cargo_Type, y = Conc, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c( "#993333","#660000","#663333", "#996666","#CC9999","#330000"))

p+theme_classic()+scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, by = 0.5),expand = c(0.01, 0.01))+
  scale_x_discrete("",expand = c(0.08, 0.08))+labs(title="")+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=8, color = "black"),
        legend.text  = element_text(size=7, color = "black"),
        plot.title = element_text(size=11, color = "black", hjust=0.6, vjust=5),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        legend.position  = 'none')+ ylab(expression(paste("Receptors per cell per second")))+
  coord_flip()

ggsave("FigS7H_R2_Flux_4h.png", width=5, height=4, dpi=400)

####################################################################################
###FIGURE S7C Flux plots for NRP1 at Steady State (no ligand treatment)
####################################################################################

library(ggplot2)

df12 <- data.frame(
  Complex=rep(c("recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                "transfer from Rab4a to Rab7a", "production to surface", "internalization from surface to Rab4a"),each=6),
  
  Cargo_Type= rep(c("enters surface", "leaves surface", "enters Rab4a/5a",  "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"), each = 1),
  Conc = c(123,0,0,123,0,0,       117.8,0,0,0,0,117.8,       0,0,0,121.06,121.06,0,
           0,0,0,3.74,0,0,         3.7,0,0,0,0,0,           0,244.2,244.2,0,0,0))

df12$Complex <- factor(df12$Complex, levels = c("transfer from Rab4a to Rab7a","production to surface","recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                                "internalization from surface to Rab4a"))

df12$Cargo_Type <- factor(df12$Cargo_Type, levels=c("enters surface", "leaves surface", "enters Rab4a/5a",  
                                                    "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"))
head(df12)

p <- ggplot(df12, aes(x = Cargo_Type, y = Conc, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c( "#339966", "#669966","#006633","#66CC99","#99CC99","#003300"))

p+theme_classic()+scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 50),expand = c(0.015, 0.015))+
  scale_x_discrete("",expand = c(0.08, 0.08))+labs(title="")+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=8, color = "black"),
        legend.text  = element_text(size=7, color = "black"),
        plot.title = element_text(size=11, color = "black", hjust=0.6, vjust=5),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        legend.position  = "right")+
  ylab(expression(paste("Receptors per cell per second")))+
  coord_flip()

ggsave("Fig7C_N1_Flux_Steady_State.png", width=5, height=4, dpi=400)

####################################################################################
###FIGURE S7F Flux plots for NRP1 after 1 hours of 50 ng/mL VEGF-A165a Treatment
####################################################################################

library(ggplot2)

df12 <- data.frame(Complex=rep(c("recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                "transfer from Rab4a to Rab7a", "production to surface", "internalization from surface to Rab4a"),each=6),
  
  Cargo_Type = rep(c("enters surface", "leaves surface", "enters Rab4a/5a",  "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"), each = 1),
        Conc = c(6.3473,0,0,6.3473,0,0,       20.464,0,0,0,0,20.464,       0,0,0,20.179,20.179,0,
                 0,0,0,0.19839,0,0,         0.459,0,0,0,0,0,           0,26.674,26.674,0,0,0))

df12$Complex <- factor(df12$Complex, levels = c("transfer from Rab4a to Rab7a","production to surface","recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                                "internalization from surface to Rab4a"))

df12$Cargo_Type <- factor(df12$Cargo_Type, levels=c("enters surface", "leaves surface", "enters Rab4a/5a",  
                                                    "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"))
head(df12)

p <- ggplot(df12, aes(x = Cargo_Type, y = Conc, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c( "#339966", "#669966","#006633","#66CC99","#99CC99","#003300"))

p+theme_classic()+scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by =10),expand = c(0.015, 0.015))+
  scale_x_discrete("",expand = c(0.08, 0.08))+labs(title="")+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=8, color = "black"),
        legend.text  = element_text(size=7, color = "black"),
        plot.title = element_text(size=11, color = "black", hjust=0.6, vjust=5),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        legend.position  = 'none')+
  ylab(expression(paste("Receptors per cell per second")))+
  coord_flip()

ggsave("FigS7F_N1_Flux_1h.png", width=5, height=4, dpi=400)

####################################################################################
###FIGURE S7I Flux plots for NRP1 after 4 hours of 50 ng/mL VEGF-A165a Treatment
####################################################################################

library(ggplot2)

df12 <- data.frame(Complex=rep(c("recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                "transfer from Rab4a to Rab7a", "production to surface", "internalization from surface to Rab4a"),each=6),
  
  Cargo_Type = rep(c("enters surface", "leaves surface", "enters Rab4a/5a",  "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"), each = 1),
        Conc = c(6.4176,0,0,6.4176,0,0,       20.401,0,0,0,0,20.401,       0,0,0,20.464,20.464,0,
                 0,0,0,0.18084,0,0,         0.459,0,0,0,0,0,           0,27.064,27.064,0,0,0))

df12$Complex <- factor(df12$Complex, levels = c("transfer from Rab4a to Rab7a","production to surface","recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                                "internalization from surface to Rab4a"))

df12$Cargo_Type <- factor(df12$Cargo_Type, levels=c("enters surface", "leaves surface", "enters Rab4a/5a",  
                                                    "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"))
head(df12)

p <- ggplot(df12, aes(x = Cargo_Type, y = Conc, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c( "#339966", "#669966","#006633","#66CC99","#99CC99","#003300"))

p+theme_classic()+scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by =10),expand = c(0.015, 0.015))+
  scale_x_discrete("",expand = c(0.08, 0.08))+labs(title="")+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=8, color = "black"),
        legend.text  = element_text(size=7, color = "black"),
        plot.title = element_text(size=11, color = "black", hjust=0.6, vjust=5),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        legend.position  = 'none')+
  ylab(expression(paste("Receptors per cell per second")))+
  coord_flip()

ggsave("FigS7I_N1_Flux_4h.png", width=5, height=4, dpi=400)

############################################################## Rab KD Simulations ########################################
#################
## Figure S8 A 
#################

#Surf_1h165a

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8A_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 140), breaks = seq(0,140, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8A.png", height = 3, width = 3, dpi=400)

#########################################
## Figure S8 B
#######################################

#knock down experiments and simulation barplots
#Internal_1h165a

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8B_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 140), breaks = seq(0,140, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8B.png", height = 3, width = 3, dpi=400)

####################
## Figure S8 C
#####################

#Total_1h165a

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8C_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 140), breaks = seq(0,140, by = 20),expand = c(0, 0))+labs(y = "Whole cell receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8C.png", height = 3, width = 3, dpi=400)

#########################################
## Figure S8 D
#######################################

#Surf_4h165a

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8D_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 140), breaks = seq(0,140, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8D.png", height = 3, width = 3, dpi=400)

#########################################
## Figure S8 E
#######################################

#Internal_4h165a

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8E_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 140), breaks = seq(0,140, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8E.png", height = 3, width = 3, dpi=400)

#########################################
## Figure S8 F
#######################################

#Total_4h165a

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8F_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 140), breaks = seq(0,140, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8F.png", height = 3, width = 3, dpi=400)

#########################################
## Figure S8 G
#######################################

#Surface_1hPLGF1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8G_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 101), breaks = seq(0,101, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8G.png", height = 3, width = 3, dpi=400)

#########################################
## Figure S8 H
#######################################

#Internal_1hPLGF1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8H_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 101), breaks = seq(0,101, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8H.png", height = 3, width = 3, dpi=400)

#########################################
## Figure S8 I
#######################################

#Total_1hPLGF1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8I_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 101), breaks = seq(0,101, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8I.png", height = 3, width = 3, dpi=400)

#########################################
## Figure S8 J
#######################################

#Surface_4hPLGF1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8J_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 101), breaks = seq(0,101, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8J.png", height = 3, width = 3, dpi=400)

#########################################
## Figure S8 K
#######################################

#Internal_4hPLGF1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8K_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 101), breaks = seq(0,101, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8K.png", height = 3, width = 3, dpi=400)

#########################################
## Figure S8 L
#######################################

#Total_4hPLGF1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS8L_data.csv", sep=",", header=T)
df7<- data.frame(df)

head(df7)

df7$Receptors <- factor(df7$Receptors, levels=c("R1", "R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df7, aes(x = Receptors, y = values, color = Receptors),binwidth=2)+
  geom_col(aes(y = values, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 101), breaks = seq(0,101, by = 20),expand = c(0, 0))+labs(y = "Surface receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#003399",
                                "#FF9999", "#660000",
                                "#99CC99", "#003300"))+
  
  scale_fill_manual(values = c("#99CCFF","#003399",
                               "#FF9999", "#660000",
                               "#99CC99", "#003300"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(  panel.background = element_blank(),
                                                     panel.grid.major.x = element_blank(),
                                                     panel.grid.major.y = element_blank(),
                                                     panel.grid.minor.x = element_blank(),
                                                     panel.grid.minor.y = element_blank(),
                                                     legend.title = element_blank(),
                                                     legend.text  = element_blank(),
                                                     legend.position = "none",
                                                     plot.title = element_blank(),
                                                     axis.title.x = element_text(size = 6, color = "black"),
                                                     axis.title.y = element_text(size = 7, color = "black"),
                                                     axis.ticks.x = element_blank(),
                                                     axis.ticks.length = unit(0.7, "mm"),
                                                     # width of tick marks in mm
                                                     axis.ticks.y = element_line(size = .2, colour = "black"),
                                                     axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                     axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S8L.png", height = 3, width = 3, dpi=400)

##############################
########## Figure S11 A
##############################

#Ligand dose gradient barplots 1h 4h
#121a_dosing_Surf_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS11A_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0,7000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999" ))+
  
  
  
  scale_fill_manual(values = c("#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999" ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S11A.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S11 B
##############################

#Ligand dose gradient barplots 1h 4h
#121a_dosing_Internal_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS11B_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))


bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0,7000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999" ))+
  
  scale_fill_manual(values = c("#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999" ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S11B.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S11 C
##############################

#Ligand dose gradient barplots 1h 4h
#121a_dosing_Total_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS11C_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0,7000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999" ))+
  
  scale_fill_manual(values = c("#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999" ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S11C.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S11 D
##############################

#Ligand dose gradient barplots 1h 4h
#165a_dosing_Surf_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS11D_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0,7000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000" ))+

  
  scale_fill_manual(values = c("#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000"))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S11D.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S11 E
##############################

#Ligand dose gradient barplots 1h 4h
#165a_dosing_Internal_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS11E_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0,7000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000" ))+
  
  scale_fill_manual(values = c("#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000"))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S11E.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S11 F
##############################

#Ligand dose gradient barplots 1h 4h
#165a_dosing_Total_R2

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS11F_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))


bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0,7000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000"))+

  scale_fill_manual(values = c("#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000"))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S11F.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 A
##############################

#Ligand dose gradient barplots 1h 4h
#121a_dosing_Surf_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12A_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 18000), breaks = seq(0,18000, by = 2000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999"))+

  scale_fill_manual(values = c("#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999"))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("FigSure S12A.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 B
##############################

#Ligand dose gradient barplots 1h 4h
#121a_dosing_Internal_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12B_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 18000), breaks = seq(0,18000, by = 2000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999" ))+
  
  scale_fill_manual(values = c("#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999" ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12B.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 C
##############################

#Ligand dose gradient barplots 1h 4h
#121a_dosing_Total_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12C_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 18000), breaks = seq(0,18000, by = 2000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999",
                                "#FF9999" ))+
  
  scale_fill_manual(values = c("#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999",
                               "#FF9999" ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12C.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 D
##############################

#Ligand dose gradient barplots 1h 4h
#165a_dosing_Surf_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12D_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 18000), breaks = seq(0,18000, by = 2000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000" ))+
  
  scale_fill_manual(values = c("#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000"  ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12D.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 E
##############################

#Ligand dose gradient barplots 1h 4h
#165a_dosing_Internal_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12E_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 18000), breaks = seq(0,18000, by = 2000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000"  ))+
  
  scale_fill_manual(values = c("#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000"  ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12E.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 F
##############################

#Ligand dose gradient barplots 1h 4h
#12L_165a_dosing_Total_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12F_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 18000), breaks = seq(0,18000, by = 2000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000",
                                "#990000"  ))+
  
  scale_fill_manual(values = c("#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000",
                               "#990000"  ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12F.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 G
##############################

#Ligand dose gradient barplots 1h 4h
#PLGF1_dosing_Surf_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12G_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 15000), breaks = seq(0,15000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF" ))+
  
  scale_fill_manual(values = c("#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF" ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12G.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 H
##############################

#Ligand dose gradient barplots 1h 4h
#PLGF1_dosing_Internal_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12H_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 15000), breaks = seq(0,15000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF" ))+
  
  scale_fill_manual(values = c("#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF" ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12H.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 I
##############################

#Ligand dose gradient barplots 1h 4h
#PLGF1_dosing_Total_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12I_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 15000), breaks = seq(0,15000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF",
                                "#99CCFF"  ))+
  
  scale_fill_manual(values = c("#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF",
                               "#99CCFF" ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12I.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 J
##############################

#Ligand dose gradient barplots 1h 4h
#PLGF2_dosing_Surf_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12J_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0,7000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366" ))+
  
  scale_fill_manual(values = c("#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366"))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12J.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 K
##############################

#Ligand dose gradient barplots 1h 4h
#PLGF2_dosing_Internal_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12K_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0,7000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366"))+
  
  scale_fill_manual(values = c("#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366"))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12K.png", height = 2, width = 3, dpi=400)

##############################
########## Figure S12 L
##############################

#Ligand dose gradient barplots 1h 4h
#PLGF2_dosing_Total_R1

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS12L_data.csv", sep=",", header=T)
df1<- data.frame(df)

head(df1)

df1$Rtype <- factor(df1$Rtype, levels=c("R2 Sim A",
                                        "R2 Sim B",
                                        "R2 Sim C",
                                        "R2 Sim D",
                                        "R2 Sim E",
                                        "R2 Sim F",
                                        "R2 Sim G",
                                        "R2 Sim H",
                                        "R2 Sim I",
                                        "R2 Sim J"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0,7000, by = 1000),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366",
                                "#003366"))+
  
  scale_fill_manual(values = c("#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366",
                               "#003366" ))

bp+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
                                                  panel.grid.major.x = element_blank(),
                                                  panel.grid.major.y = element_blank(),
                                                  panel.grid.minor.x = element_blank(),
                                                  panel.grid.minor.y = element_blank(),
                                                  legend.title = element_blank(),
                                                  legend.text  = element_blank(),
                                                  legend.position = "none",
                                                  plot.title = element_blank(),
                                                  axis.title.x = element_text(size = 8, color = "black"),
                                                  axis.title.y = element_text(size = 8, color = "black"),
                                                  axis.ticks.x = element_blank(),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Figure S12L.png", height = 2, width = 3, dpi=400)

########## THE End ###################################################
