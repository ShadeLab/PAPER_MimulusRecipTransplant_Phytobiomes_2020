## Plant mass plots and statistics

#load metadata
#remove "#" from first line before importing
#remove "BarcodeSequence" column before importing.
#remove "LinkerPrimerSequence" column before importing.
full_mimulus_data<-read.table("Mimulus_metadata_FULL.txt", sep='\t', header=T, row.names=1)

#remove bulk soil samples from metadata list
full_mimulus_data<-full_mimulus_data[!(full_mimulus_data$Rhizo_Bulk=="Bulk"),]

#subset to contain only shoot mass data and variables of interest
mimulus_shoots <- subset(full_mimulus_data, select=c("Genotype", "Origin", "Site_planted", "Shoot_Mass_g"))

#melt and plot shoot mass data
library(reshape2)
library(ggplot2)
mimulus_shoots_m<-melt(mimulus_shoots)

#change order of factors to plot them appropriately
mimulus_shoots_m$Genotype <- factor(mimulus_shoots_m$Genotype,levels = c('LMC','OCC','MRR','SWB'),ordered = TRUE)

shootfig<-ggplot(mimulus_shoots_m, aes(Genotype, value, fill=Origin))+
  geom_boxplot(outlier.shape=NA)+
  facet_grid(.~Site_planted)+
  theme_bw()+
  scale_fill_discrete(name="Ecotype")+
  scale_fill_manual(values=c("royalblue1","green"))+
  ylab("Shoot Mass (g)")+
  coord_cartesian(ylim=c(0,3))+
  ggtitle("Coastal site                                 Inland site")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(axis.text.x=element_blank(),axis.text=element_text(size=15, colour="black"), axis.title=element_text(size=20, colour="black"), strip.text = element_text(size = 12),legend.position="none")+
  xlab("")
shootfig


#Now plotting root mass data.

#subset full dataset to contain only root mass data and variables of interest
mimulus_roots<- subset(full_mimulus_data, select=c("Genotype", "Origin", "Site_planted","Root_Mass_mg"))

#melt and plot root mass data
library(reshape2)
library(ggplot2)
mimulus_roots_m<-melt(mimulus_roots)

#change order of factors to plot them appropriately
mimulus_roots_m$Genotype <- factor(mimulus_roots_m$Genotype,levels = c('LMC','OCC','MRR','SWB'),ordered = TRUE)

rootfig<-ggplot(mimulus_roots_m, aes(Genotype, value, fill=Origin))+
  geom_boxplot(outlier.shape=NA)+
  #geom_point()+
  facet_grid(. ~ Site_planted)+
  theme_bw()+
  scale_fill_manual(labels = c("Coastal ecotype", "Inland ecotype"), values=c("royalblue1","green"))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  ylab("Root Mass (mg)")+
  coord_cartesian(ylim=c(0,125))+
  theme(axis.text=element_text(size=15, colour="black"), axis.title=element_text(size=20, colour="black"), strip.text = element_text(size = 12),legend.position="bottom")+
  xlab("")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  theme(legend.title=element_blank())
rootfig

#combine the two figures and save as a single file
library("cowplot")
FigureS1<-plot_grid(shootfig, rootfig, align = "v", nrow = 2, rel_heights=c(1,1.35))
FigureS1
ggsave("FigureS1.tiff",FigureS1,width=6,height=6, units="in",dpi=300)






##Statistical tests

#ANOVA for shoot mass

#transform shoot mass
mimulus_shoots$trans_Shoot_Mass_g<- log10(mimulus_shoots$Shoot_Mass_g)
#Do ANOVA to compare shoot mass
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
shoot_fit <- aov(trans_Shoot_Mass_g ~ Origin*Site_planted, data=mimulus_shoots, na.action=na.exclude)
drop1(shoot_fit,~.,test="F") # type III SS and F Tests

#TESTING ASSUMPTIONS
#check whether residuals are normally distributed
shoot_resids <- residuals(shoot_fit)
shoot_preds <- predict(shoot_fit)
plot(shoot_resids ~ shoot_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Test for normality
shapiro.test(shoot_resids)
qqnorm(shoot_resids)
qqline(shoot_resids)
library("car")
leveneTest(trans_Shoot_Mass_g ~ Origin*Site_planted, data=mimulus_shoots, na.action=na.exclude)
TukeyHSD(shoot_fit)






#ANOVA for root mass



#transform root mass
mimulus_roots$trans_Root_Mass_mg<- (mimulus_roots$Root_Mass_mg)
#Do ANOVA to compare Root mass
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Root_fit <- aov(trans_Root_Mass_mg ~ Origin*Site_planted, data=mimulus_roots, na.action=na.exclude)
drop1(Root_fit,~.,test="F") # type III SS and F Tests

#TESTING ASSUMPTIONS
#check whether residuals are normally distributed
Root_resids <- residuals(Root_fit)
Root_preds <- predict(Root_fit)
plot(Root_resids ~ Root_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)

text(Root_resids ~ Root_preds, labels=rownames(mimulus_roots),data=mimulus_roots, cex=0.9, font=2)

#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Test for normality
shapiro.test(Root_resids)
qqnorm(Root_resids)
qqline(Root_resids)
library("car")
leveneTest(trans_Root_Mass_mg ~ Origin*Site_planted, data=mimulus_roots, na.action=na.exclude)
TukeyHSD(Root_fit)



