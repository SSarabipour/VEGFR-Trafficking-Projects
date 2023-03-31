################################################################################################################################
######### R code for FIGUREs visualized for Sarabipour et al (2022) bioRxiv preprint https://doi.org/10.1101/2022.09.30.510412
################################################################################################################################
# Highlight code for desired figures and use ctrl+enter to run


##############################
######## FIGURE 3
#############################

# Figure 3 uses simulation results based on parameters from 
# 2015-2017 data from Clegg LE & Mac Gabhann F (2015) and (2017)
# https://doi.org/10.1371/journal.pcbi.1004158
# https://doi.org/10.1371/journal.pcbi.1005445

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig3B_data.csv", sep=",", header=T)
df1<- data.frame(df)

df1$Rtype <- factor(df1$Rtype, levels=c("R1 Expt","R1 Sim","R2 Expt","R2 Sim", "N1 Expt", "N1 Sim"))

bp <- ggplot(df1, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  geom_errorbar(aes(ymin=Rvalue-RSD, ymax=Rvalue+RSD, color = Rcomp),width=0.5, alpha=2, size=0.5)+    
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100, by = 25),expand = c(0.005, 0.005))+
  labs(y = "Receptors by location (% of total)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#66CC99","#006633","#003300","#003300",
                                "#9999FF","#6666FF","#0000FF","#0000FF", 
                                "#FFCCCC","#CC6666","#CC0000","#CC0000" ))+
  
  scale_fill_manual(values = c("#66CC99","#006633","#003300","#003300",
                               "#9999FF","#6666FF","#0000FF","#0000FF",
                               "#FFCCCC","#CC6666","#CC0000","#CC0000" ))

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

ggsave ("Fig3B.png", height = 3, width = 4.5, dpi=400)


df<-read.csv("Fig3C_data.csv", sep=",", header=T)
df2<- data.frame(df)

p = ggplot(data=df2, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType, linetype=RecType),size=1)+
  scale_linetype_manual(values=c("solid","solid","solid"))+ #dashed
  
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
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=8, color = "black"),
        
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(0.8, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-6,-6,-6,-10),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#339966","#0000FF","#CC0000"))+
  labs(title = "", x = "Time (min)", y = "Total Receptors (% of control)")+ theme(legend.title=element_blank())+
  scale_y_continuous("Total Receptors (% of control)",limits = c(0, 130), breaks = seq(0, 130, by = 20),expand = c(0, 0))+
  scale_x_continuous("Time (min)",limits = c(0, 250), breaks = seq(0, 240, by = 30),expand = c(0, 0))

ggsave ("Fig3C.png", height = 3, width = 3.5, dpi=400)


########################################################################################
###### FIGURE 4 
########################################################################################

#Box plot optimized unligated parameters

library(ggplot2)
library(dplyr)
library(reshape2)

Fits<-read.csv("Fig4_data.csv", sep=",", header=T)
df3<- data.frame(Fits)

# lock in factor level order
df3$parameter <- as.character(df3$parameter)
df3$parameter <- factor(df3$parameter , levels = unique(df3$parameter))

df3$parameter <- factor(df3$parameter,levels=c( "kprod(N1)",  "kprod(R2)",  "kprod(R1)",
                                                 "k4to11(N1)", "k4to11(R2)", "k4to11(R1)", 
                                                 "krec11(N1)", "krec11(R2)", "krec11(R1)", 
                                                 "krec4(N1)",  "krec4(R2)",  "krec4(R1)", 
                                                 "kdeg(N1)",   "kdeg(R2)",   "kdeg(R1)", 
                                                 "kint(N1)",   "kint(R2)",   "kint(R1)"))

table_of_colors = c( "kprod(N1)" = "dark green","kprod(R2)"= "dark red",  "kprod(R1)"= "blue",
                    "k4to11(N1)"= "dark green", "k4to11(R2)"= "dark red", "k4to11(R1)"="blue", 
                    "krec11(N1)"= "dark green", "krec11(R2)"= "dark red", "krec11(R1)"="blue", 
                    "krec4(N1)"= "dark green",  "krec4(R2)"= "dark red",  "krec4(R1)"="blue", 
                    "kdeg(N1)"= "dark green",   "kdeg(R2)"= "dark red",   "kdeg(R1)"="blue", 
                    "kint(N1)"= "dark green",   "kint(R2)"= "dark red",   "kint(R1)"="blue")

#box plot
p <- ggplot(df3, aes(y = parameter, x = value)) + 
  geom_boxplot(width=0.7, color = table_of_colors)+
  theme(plot.margin = margin(30, 30, 30, 30),
        plot.title = element_text(hjust = 0.5, vjust = 5, size = 14), 
        axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = -0.5, size = 10),
        axis.title.x = element_text(margin = margin(9, 9, 9, 9),  size = 11),
        axis.text.y = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.background = element_blank(), 
        legend.position = "none")+
  theme(strip.background=element_rect(fill="white"))+
  labs(title = "", x = 'log(k)', y = "") +labs(title="Unligated receptor trafficking parameters")+
  scale_x_continuous(limits = c(-6.5, -0.5), breaks = seq(-6.5, -0.5, by = 1),expand = c(0, 0))

ggsave ("Fig4.png", height = 7, width = 6, dpi=400)



#################################################################
####### FIGURE 5B-C
################################################################

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig5B_data.csv", sep=",", header=T)
df4<- data.frame(df)

df4$Rtype <- factor(df4$Rtype, levels=c("R1 Expt","R1 Sim","R2 Expt","R2 Sim", "N1 Expt", "N1 Sim"))

bp2 <- ggplot(df4, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  geom_errorbar(aes(ymin=Rvalue-RSD, ymax=Rvalue+RSD, color = Rcomp),width=0.5, alpha=2, size=0.5)+    
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100, by = 10),expand = c(0.005, 0.005))+
  labs(y = "Receptors by location (% of total)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#66CC99","#006633","#003300", 
                                "#003300",
                                "#9999FF","#6666FF","#0000FF",
                                "#0000FF", 
                                "#FFCCCC", "#CC6666","#CC0000",
                                "#CC0000"))+
  
  scale_fill_manual(values = c(  "#66CC99","#006633","#003300", 
                                 "#003300",
                                 "#9999FF","#6666FF","#0000FF",
                                 "#0000FF",
                                 "#FFCCCC", "#CC6666","#CC0000",
                                 "#CC0000"))

bp2+scale_x_discrete(expand = c(0.08, 0.08))+theme(panel.background = element_blank(),
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
                                                   axis.ticks.length = unit(0.7, "mm"),# width of tick marks in mm
                                                   axis.ticks.y = element_line(size = .2, colour = "black"),
                                                   axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                   axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Fig5B.png", height = 3, width = 4.5, dpi=400)



df<-read.csv("Fig5C_data.csv", sep=",", header=T)
df5<- data.frame(df)

# lock in factor level order
p = ggplot(data=df5, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType),size=1,position=position_jitter(w=0, h=0))+
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType),width=10, alpha=2, size=0.9,
                                                                     position=position_dodge(2))+
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
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=8, color = "black"),
      
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(0.8, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-6,-6,-6,-10),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#339966","#0000FF","#CC0000"))+
  labs(title = "", x = "Time (min)", y = "Total Receptors (% of control)")+ theme(legend.title=element_blank())+
scale_y_continuous("Total Receptors (% of control)",limits = c(0, 130), breaks = seq(0, 130, by = 20),expand = c(0, 0))+
scale_x_continuous("Time (min)",limits = c(0, 250), breaks = seq(0, 240, by = 30),expand = c(0, 0))

ggsave ("Fig5C.png", height = 3, width = 3.5, dpi=400)

#######################################################################
########### FIGURE 5E-F
#######################################################################

#Rab4a and/or Rab11a knock-downs - surface receptor levels, unligated HUVECs

library(ggplot2)
library(dplyr)
library(reshape2)

#Surface Receptor levels at time = 0 normalized to no treatment zero time point receptor

df<-read.csv("Fig5EF_data.csv", sep=",", header=T)
df6<- data.frame(df)

df6$Receptors <- factor(df6$Receptors, levels=c("R1", "R1 siRNA Rab4a","R1 siRNA Rab11a","R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a", "R2 siRNA Rab11a", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a", "N1 siRNA Rab11a", "N1 siRNA Rab4a11a"))

bp <- ggplot(df6, aes(x = Receptors, y = surfvalues, color = Receptors),binwidth=2)+
  geom_col(aes(y = surfvalues, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100, by = 20), expand = c(0, 0))+labs(y = "Surface receptors\n(% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#6699CC","#336699","#003399",
                                "#FF9999","#CC6666","#990000", "#660000",
                                "#99CC99","#669966", "#006600", "#003300"))+
  scale_fill_manual(values = c("#99CCFF","#6699CC","#336699","#003399",
                               "#FF9999","#CC6666","#990000", "#660000",
                               "#99CC99","#669966", "#006600", "#003300"))

bp+scale_x_discrete(expand = c(0.05, 0.05))+theme(panel.background = element_blank(),
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
                                                  axis.ticks.length = unit(0.7, "mm"),# width of tick marks in mm 
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("Fig5E.png", height = 3, width = 3.5, dpi=400)


bp2 <- ggplot(df6, aes(x = Receptors, y = wcvalues, color = Receptors),binwidth=2)+
  geom_col(aes(y = wcvalues, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100, by = 20),expand = c(0, 0))+labs(y = "Total receptors\n(% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#6699CC","#336699","#003399",
                                "#FF9999","#CC6666","#990000", "#660000",
                                "#99CC99","#669966", "#006600", "#003300"))+
  scale_fill_manual(values = c("#99CCFF","#6699CC","#336699","#003399",
                               "#FF9999","#CC6666","#990000", "#660000",
                               "#99CC99","#669966", "#006600", "#003300"))

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

ggsave ("Fig5F.png", height = 3, width = 3.5, dpi=400)



######################################################################################
######## FIGURE 6 
######################################################################################

#Heatmap of Local Sensitivity Analysis results

library(ggplot2)
library(reshape2)
library(plyr)
library(Cairo)

data9 <- read.csv("Fig6_data.csv", header = TRUE, stringsAsFactors = FALSE)

data10 <- melt(data9,id.vars = c("Parameters"))

#rename column names
colnames(data10) <- c("Parameters", "variable","value")

head(data10)
str(data10)
data10$Parameters <- as.character(data10$Parameters)
data10$Parameters <- as.factor(data10$Parameters)
levels(data10$Parameters)
data10$Parameters <- factor(data10$Parameters, levels = levels(data10$Parameters)[c(2,4,3,5,6,1,14,16,15,17,18,13,8,10,9,11,12,7)])

p <- ggplot(data10,aes(x=variable,y=Parameters,fill=value))+
  geom_tile()+ 
  scale_fill_gradientn(colors=c("purple","white","dark green"),guide="colorbar",
                       #same midpoint for plots (mean of the range)
                       breaks=c(-1,-0.5,0,0.5,1), #breaks in the scale bar
                       limits=c(-1,1))+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=3),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.title = element_text(size = 5, color = "black"),
        axis.title.x = element_text(size = 5, color = "black"),
        axis.title.y = element_text(size = 5, color = "black"),
        axis.text.x=element_text(size=4),
        axis.text.y=element_text(size=4))+
  
  labs(title= "Sensitivity analysis", y="Parameters", x = "Output")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1),
        plot.title = element_text(size = 5.5,  vjust = 1.05, hjust=0.5, colour = "black"))

ggsave("Fig6.png", width=3, height=3, dpi=400)


##########################################################################################
################# FIGURE 7A-C 
##########################################################################################

#Flux plots (stacked barplots) for VEGFR1, VEGFR2, Neuropilin1 at Steady State cell culture condition


library(ggplot2)

data10 <- read.csv("Fig7_data.csv", header = TRUE, stringsAsFactors = FALSE)
df10 <- data.frame(data10)
df10$Complex <- factor(df10$Complex, levels = c("transfer from Rab4a to Rab7a","production to surface","recycling from Rab4a to surface","recycling from Rab11a to surf", "transfer from Rab4a to Rab11a",
                                            "internalization from surface to Rab4a"))
df10$Cargo_Type <- factor(df10$Cargo_Type, levels=c("enters surface", "leaves surface", "enters Rab4a/5a",  
                                                "leaves Rab4a/5a", "enters Rab11a", "leaves Rab11a"))

p <- ggplot(df10, aes(x = Cargo_Type, y = R1rate, fill=Complex), binwidth=0)+
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
        legend.position = "none") + ylab(expression(paste("Receptors per cell per second"))
        )+
        coord_flip()

ggsave("Fig7_R1_Flux.png", width=5, height=4, dpi=400)


p2 <- ggplot(df10, aes(x = Cargo_Type, y = R2rate, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c( "#993333","#660000","#663333", "#996666","#CC9999","#330000"))

p2+theme_classic()+scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5),expand = c(0.015, 0.015))+
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
        legend.position  = "none")+ ylab(expression(paste("Receptors per cell per second")))+
        coord_flip()

ggsave("Fig7_R2_Flux.png", width=5, height=4, dpi=400)


p3 <- ggplot(df10, aes(x = Cargo_Type, y = N1rate, fill=Complex), binwidth=0)+
  geom_col(aes(fill = Complex), width = 0.7)+ theme(aspect.ratio = .8) + 
  scale_fill_manual(values = c( "#339966", "#669966","#006633","#66CC99","#99CC99","#003300"))

p3+theme_classic()+scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 50),expand = c(0.015, 0.015))+
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
        legend.position  = "none")+
        ylab(expression(paste("Receptors per cell per second")))+
        coord_flip()

ggsave("Fig7_N1_Flux.png", width=5, height=4, dpi=400)


##################################################################
######## Figure 8B and 8D
#################################################################

#CHQ time course for unligated HUVECs

library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("Fig8B_data.csv", sep=",", header=T)
df13<- data.frame(df)

# lock in factor level order
p = ggplot(data=df13, aes(x=Time)) + geom_point(aes(y=TotalRecExp, color=RecType),size=2)+ 
  geom_line(aes(y=TotalRecComp, color=RecType),size=1,position=position_jitter(w=0, h=0))+
  geom_errorbar(aes(ymin=TotalRecExp-sd, ymax=TotalRecExp+sd,color=RecType), 
                width=1, alpha=1, size=0.9,
                position=position_dodge(2))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 0.1, vjust = 0.15, size = 9))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  
  theme(plot.margin = margin(20, 20, 20, 20),
        axis.title.x = element_text(colour = "black", size = 8),
        axis.title.y = element_text(colour = "black", size = 8),
        legend.position = "top",
        legend.text  = element_text(margin = margin(r = 15, unit = "pt"),size=8, color = "black"),
        
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(0.8, 'cm'), #change legend key width
        legend.spacing.x = unit(0.01, 'cm'),
        legend.box.margin=margin(-6,-6,-6,-10),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(size = .2, colour = "black"),
        axis.ticks.x = element_line(size = .2, colour = "black"),
        axis.ticks.length = unit(0.7, "mm"))+scale_color_manual(values=c("#0000FF"))+
        labs(title = "", x = "Time (hours)", y = "% Total Receptor")+ theme(legend.title=element_blank())+
        scale_y_continuous("Total receptors (% of control)",limits = c(0, 250), breaks = seq(0, 250, by = 50),expand = c(0, 0))+
        scale_x_continuous("Time (hours)",limits = c(0, 24.8), breaks = seq(0, 24.8, by = 2),expand = c(0, 0))

ggsave ("Fig8B.png", height = 3, width = 3.5, dpi=400)



df<-read.csv("Fig8D_data.csv", sep=",", header=T)
df14<- data.frame(df)

df14$Rtype <- factor(df14$Rtype, levels=c("R1 Expt A",
                                          "R1 Expt B",
                                          "R1 Sim A",
                                          "R1 Sim B",
                                          "R2 Expt C",
                                          "R2 Expt D",
                                          "R2 Sim C",
                                          "R2 Sim D"))

bp <- ggplot(df14, aes(x = Rtype, y = Rvalue, color = Rcomp),binwidth=2)+
  geom_col(aes(y = Rvalue, fill = Rcomp), width = 0.7)+ 
  geom_errorbar(aes(ymin=Rvalue-RSD, ymax=Rvalue+RSD, color = Rcomp),width=0.5, alpha=2, size=0.5)+    
  
  theme_classic()+xlab("")+ 
  scale_y_continuous(limits = c(0, 250), breaks = seq(0,250, by = 50),expand = c(0.005, 0.005))+
  labs(y = "Receptors (% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#0066CC",
                                "#99CCFF",
                                "#0066CC",
                                "#99CCFF",
                                "#0066CC",
                                "#99CCFF",
                                "#0066CC",
                                "#99CCFF"))+
                                
  
  scale_fill_manual(values = c("#0066CC",
                               "#99CCFF",
                               "#0066CC",
                               "#99CCFF",
                               "#0066CC",
                               "#99CCFF",
                               "#0066CC",
                               "#99CCFF"))

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

ggsave ("Fig8D_4h18hCHQ_R1.png", height = 2, width = 3, dpi=400)




################################################################################
########   Supplementary Figures    #####################################
################################################################################



#######################################################################
########  FIGURE S2B
#######################################################################

#Rab4a and/or Rab11a knock-downs then CHX time course unligated HUVECs

library(ggplot2)
library(dplyr)
library(reshape2)

#Surface Receptor levels at time = 0 normalized to no treatment zero time point receptor, also kR1N1_on=0

df<-read.csv("FigS2BC_data.csv", sep=",", header=T)
df6<- data.frame(df)

df6$Receptors <- factor(df6$Receptors, levels=c("R1", "R1 siRNA Rab4a","R1 siRNA Rab11a","R1 siRNA Rab4a11a",
                                                "R2", "R2 siRNA Rab4a", "R2 siRNA Rab11a", "R2 siRNA Rab4a11a",
                                                "N1", "N1 siRNA Rab4a", "N1 siRNA Rab11a", "N1 siRNA Rab4a11a"))

bp2 <- ggplot(df6, aes(x = Receptors, y = surfvalues, color = Receptors),binwidth=2)+
  geom_col(aes(y = surfvalues, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100, by = 25),expand = c(0, 0))+labs(y = "Surface receptors\n(% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#6699CC","#336699","#003399",
                                "#FF9999","#CC6666","#990000", "#660000",
                                "#99CC99","#669966", "#006600", "#003300"))+
  scale_fill_manual(values = c("#99CCFF","#6699CC","#336699","#003399",
                               "#FF9999","#CC6666","#990000", "#660000",
                               "#99CC99","#669966", "#006600", "#003300"))

bp2+scale_x_discrete(expand = c(0.05, 0.05))+theme(panel.background = element_blank(),
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
                                                   axis.ticks.length = unit(0.7, "mm"),# width of tick marks in mm 
                                                   axis.ticks.y = element_line(size = .2, colour = "black"),
                                                   axis.text.x = element_text(color = "black",angle = 45, hjust = 0.9, vjust = 0.85, size = 6),
                                                   axis.text.y = element_text(size = 6, color = "black"))

ggsave ("SupplFigS2B.png", height = 3, width = 3.5, dpi=400)

#Total Receptor levels at time = 0 normalized to no treatment zero time point receptor, also kR1N1_on=0

bp3 <- ggplot(df6, aes(x = Receptors, y = wcvalues, color = Receptors),binwidth=2)+
  geom_col(aes(y = wcvalues, fill = Receptors), width = 0.7)+ 
  theme_classic()+theme(aspect.ratio = .6)+xlab("")+ 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100, by = 25),expand = c(0, 0))+labs(y = "Total receptors\n(% of control)")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#99CCFF","#6699CC","#336699","#003399",
                                "#FF9999","#CC6666","#990000", "#660000",
                                "#99CC99","#669966", "#006600", "#003300"))+
  scale_fill_manual(values = c("#99CCFF","#6699CC","#336699","#003399",
                               "#FF9999","#CC6666","#990000", "#660000",
                               "#99CC99","#669966", "#006600", "#003300"))

bp3+scale_x_discrete(expand = c(0.05, 0.05))+theme(  panel.background = element_blank(),
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

ggsave ("SupplFigS2C.png", height = 3, width = 3.5, dpi=400)




######################################################################################
######################  Supplementary FIGURE 3
######################################################################################


library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS3_data.csv", sep=",", header=T)
df15<- data.frame(df)

bp <- ggplot(df15, aes(x = kR1N1, y = values),binwidth=2)+
  geom_point(aes(x = kR1N1, y = values), color='#006600', fill='#006600', shape=21)+ 
  
  theme_classic()+xlab("log(R1-N1 coupling rate constant)")+ 
  scale_y_continuous(limits = c(0.000001, 0.005), breaks = seq(0.000001, 0.005, by = 0.001),expand = c(0.0001, 0.0001))+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(limits = c(-7, -0.5), breaks = seq(-7,-0.5, by = 0.5),expand = c(0.005, 0.005))+
  labs(y = "kN1 prod")+
  theme( legend.position = "none")+
  scale_color_manual(values = c("#006600"))+scale_fill_manual(values = c("#006600"))

bp+theme(panel.background = element_blank(),
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
                                                  
                                                  axis.ticks.x = element_line(size = .2, colour = "black"),
                                                  axis.ticks.y = element_line(size = .2, colour = "black"),
                                                  axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
                                                  
                                                  axis.text.x = element_text(color = "black", size = 6),
                                                  axis.text.y = element_text(size = 6, color = "black"))

ggsave ("SupplFigS3.png", height = 2, width = 3.2, dpi=400)



######################################################################################
########## Supplementary FIGURE 4 A-B
######################################################################################


library(ggplot2)
library(dplyr)
library(reshape2)

df<-read.csv("FigS4_data.csv", sep=",", header=T)
df15<- data.frame(df)

bp <- ggplot(df15, aes(x = kR1N1, y = wcvalues),binwidth=2)+
  geom_point(aes(x = kR1N1, y = wcvalues), color='#006600', fill='#006600', shape=21)+ 
  geom_point(aes(x = kR1N1, y = rabvalues), color='#006600', fill='#FFFFFF', shape=21)+ 
  
  theme_classic()+xlab("log(R1-N1 coupling rate constant)")+ 
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 20),expand = c(0.0001, 0.0001))+
  scale_x_continuous(limits = c(-7, -0.5), breaks = seq(-7,-0.5, by = 1.0),expand = c(0.005, 0.005))+
  labs(y = "Whole cell NRP1 levels\n(,000 receptors per cell)")+
  theme( legend.position = "none")
#+  scale_color_manual(values = c("#006600"))+scale_fill_manual(values = c("#006600"))

bp+theme(panel.background = element_blank(),
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
         
         axis.ticks.x = element_line(size = .2, colour = "black"),
         axis.ticks.y = element_line(size = .2, colour = "black"),
         axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
         
         axis.text.x = element_text(color = "black", size = 6),
         axis.text.y = element_text(size = 6, color = "black"))

ggsave ("SupplFigS4A.png", height = 2, width = 3.2, dpi=400)

bp2 <- ggplot(df15, aes(x = kR1N1, y = ratiovalues),binwidth=2)+
  geom_point(aes(x = kR1N1, y = ratiovalues), color='#006600', fill='#006600', shape=21)+ 

  theme_classic()+xlab("log(R1-N1 coupling rate constant)")+ 
  scale_y_continuous(limits = c(0.6, 1.4), breaks = seq(0.6, 1.4, by = 0.2),expand = c(0.0001, 0.0001))+
  scale_x_continuous(limits = c(-7, -0.5), breaks = seq(-7,-0.5, by = 1.0),expand = c(0.005, 0.005))+
  labs(y = "Ratio of whole cell NRP1 levels\n(post-Rab-knockdown vs control)")+
  theme( legend.position = "none")

bp2+theme(panel.background = element_blank(),
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
         
         axis.ticks.x = element_line(size = .2, colour = "black"),
         axis.ticks.y = element_line(size = .2, colour = "black"),
         axis.ticks.length = unit(0.7, "mm"), # width of tick marks in mm
         
         axis.text.x = element_text(color = "black", size = 6),
         axis.text.y = element_text(size = 6, color = "black"))

ggsave ("SupplFigS4B.png", height = 2, width = 3.2, dpi=400)



######################################################################################
######## Supplementary FIGURE 5
######################################################################################

#Heatmap of Local Sensitivity Analysis results with kN1R1_on set to zero

library(ggplot2)
library(reshape2)
library(plyr)
library(Cairo)

data9 <- read.csv("FigS5_data.csv", header = TRUE, stringsAsFactors = FALSE)

data10 <- melt(data9,id.vars = c("Parameters"))

#rename column names
colnames(data10) <- c("Parameters", "variable","value")

head(data10)
str(data10)
data10$Parameters <- as.character(data10$Parameters)
data10$Parameters <- as.factor(data10$Parameters)
levels(data10$Parameters)
data10$Parameters <- factor(data10$Parameters, levels = levels(data10$Parameters)[c(2,4,3,5,6,1,14,16,15,17,18,13,8,10,9,11,12,7)])

p <- ggplot(data10,aes(x=variable,y=Parameters,fill=value))+
  geom_tile()+ 
  scale_fill_gradientn(colors=c("purple","white","dark green"),guide="colorbar",
                       breaks=c(-1,-0.5,0,0.5,1), #breaks in the scale bar
                       limits=c(-1,1))+
  
  theme(plot.margin = margin(30, 30, 30, 30),
        axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=3),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.title = element_text(size = 5, color = "black"),
        axis.title.x = element_text(size = 5, color = "black"),
        axis.title.y = element_text(size = 5, color = "black"),
        axis.text.x=element_text(size=4),
        axis.text.y=element_text(size=4))+
  
  labs(title= "Sensitivity analysis", y="Parameters", x = "Output")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1),
        plot.title = element_text(size = 5.5,  vjust = 1.05, hjust=0.5, colour = "black"))

p

ggsave("SupplFigS5.png", width=3, height=3, dpi=400)

########## THE END ###################################################
