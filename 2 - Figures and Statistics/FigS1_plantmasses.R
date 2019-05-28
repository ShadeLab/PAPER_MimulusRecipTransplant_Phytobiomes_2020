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

#make new variable called GenotypeSite
mimulus_shoots$GenotypeSite<-paste(mimulus_shoots$Genotype,mimulus_shoots$Site_planted)
mimulus_roots$GenotypeSite<-paste(mimulus_roots$Genotype,mimulus_roots$Site_planted)




#shoot mass:  all combos of pairwise t-tests, with FDR correction
shoot_t.tests<-pairwise.t.test(mimulus_shoots$Shoot_Mass_g, mimulus_shoots$GenotypeSite, p.adjust.method = 'fdr')
View(shoot_t.tests$p.value)
shoot_wilcox.tests<-pairwise.wilcox.test(mimulus_shoots$Shoot_Mass_g, mimulus_shoots$GenotypeSite, p.adjust.method = 'fdr')
View(shoot_wilcox.tests$p.value)

#root mass:   all combos of pairwise t-tests, with FDR correction
root_t.tests<-pairwise.t.test(mimulus_roots$Root_Mass_mg, mimulus_roots$GenotypeSite, p.adjust.method = 'fdr')
View(root_t.tests$p.value)
root_wilcox.tests<-pairwise.wilcox.test(mimulus_roots$Root_Mass_mg, mimulus_roots$GenotypeSite, p.adjust.method = 'fdr')
View(root_wilcox.tests$p.value)





#ANOVA for shoot mass

#transform shoot mass
mimulus_shoots$trans_Shoot_Mass_g<- sqrt(mimulus_shoots$Shoot_Mass_g)

#Do ANOVA to compare shoot mass
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
shoot_fit <- aov(trans_Shoot_Mass_g ~ GenotypeSite, data=mimulus_shoots, na.action=na.exclude)
drop1(shoot_fit,~.,test="F") # type III SS and F Tests

#TESTING ASSUMPTIONS
shoot_resids <- residuals(shoot_fit)
shoot_preds <- predict(shoot_fit)
plot(shoot_resids ~ shoot_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(trans_Shoot_Mass_g ~ GenotypeSite, data=mimulus_shoots, na.action=na.exclude)
library("car")
leveneTest(trans_Shoot_Mass_g ~ GenotypeSite, data=mimulus_shoots, na.action=na.exclude)
#Test for normality
shapiro.test(shoot_resids)
qqnorm(shoot_resids)
qqline(shoot_resids)
TukeyHSD(shoot_fit)





#ANOVA for root mass

#transform root mass
mimulus_roots$trans_root_Mass_g<- sqrt(mimulus_roots$root_Mass_g)

#Do ANOVA to compare root mass
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
root_fit <- aov(Root_Mass_mg ~ GenotypeSite, data=mimulus_roots, na.action=na.exclude)
drop1(root_fit,~.,test="F") # type III SS and F Tests

#TESTING ASSUMPTIONS
root_resids <- residuals(root_fit)
root_preds <- predict(root_fit)
plot(root_resids ~ root_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(Root_Mass_mg ~ GenotypeSite, data=mimulus_roots, na.action=na.exclude)
library("car")
leveneTest(Root_Mass_mg ~ GenotypeSite, data=mimulus_roots, na.action=na.exclude)
#Test for normality
shapiro.test(root_resids)
qqnorm(root_resids)
qqline(root_resids)
TukeyHSD(root_fit)



