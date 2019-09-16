## Plots and statistics: alpha diversity

library(ggplot2)
library(reshape2)

#import otu table (remove "#Constructed from biom file " from first line before importing).
mim_otu<-read.table("zotu_table_forR.txt", row.names=1, header=T)
#calculate richness (number of OTUs) in each sample
richness <-as.data.frame(colSums(mim_otu != 0))
colnames(richness)<-c("Richness")
richness<- data.frame(names = row.names(richness), richness)

#load metadata
#remove "#" from first line before importing
#remove "BarcodeSequence" column before importing.
#remove "LinkerPrimerSequence" column before importing.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, sep='\t', row.names=1)

#Verify that rows are in same order, then combine tables.
row.names(richness)
row.names(mim_envir)
mim_diversity<-cbind(mim_envir,richness)

#format alpha diversity values as numeric
mim_diversity$PD_whole_tree_alpha<-as.numeric(mim_envir$PD_whole_tree_alpha)
mim_diversity$shannon_alpha<-as.numeric(mim_envir$shannon_alpha)
mim_diversity$Richness<-as.numeric(mim_diversity$Richness)

mim_diversity$Ecotype<- ifelse(mim_diversity$SampleType =="Bulk Soil", paste("Bulk soil"), paste(mim_diversity$Origin," ","ecotype"))

#Retain only data of interest
mim_diversity_richness <- subset(mim_diversity, select=c("Genotype", "Ecotype","Origin", "Site_planted","Richness"))
#reshape the data
mim_diversity_richness_m<-melt(mim_diversity_richness)
mim_diversity_richness_m$Genotype <- factor(mim_diversity_richness_m$Genotype,levels = c('BulkSoil','LMC','OCC','MRR','SWB'),ordered = TRUE)

#plot the data
Rich<-ggplot(mim_diversity_richness_m, aes(Genotype, value, fill=Ecotype))+
  facet_grid(.~Site_planted)+
  theme_bw()+
  ylab("Species richness")+
  ylim(3000,4900)+
  scale_fill_manual(values=c("gray10","royalblue1","green"))+
  geom_boxplot(outlier.shape=NA)+
  theme(axis.text=element_text(size=10, colour="black"), axis.title=element_text(size=12, colour="black"), strip.text = element_text(size = 12),axis.text.x=element_blank(),legend.position="none")+
  xlab("")+
  ggtitle("Coastal site                                  Inland site")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())
Rich

#Retain only data of interest
mim_diversity_PD <- subset(mim_diversity, select=c("Genotype", "Ecotype","Origin", "Site_planted","PD_whole_tree_alpha"))
#reshape the data
mim_diversity_PD_m<-melt(mim_diversity_PD)
mim_diversity_PD_m$Genotype <- factor(mim_diversity_PD_m$Genotype,levels = c('BulkSoil','LMC','OCC','MRR','SWB'),ordered = TRUE)

#plot the data
PD<-ggplot(mim_diversity_PD_m, aes(Genotype, value, fill=Ecotype))+
  facet_grid(.~Site_planted)+
  theme_bw()+
  ylab("Phylogenetic Diversity")+
  ylim(140,230)+
  scale_fill_manual(values=c("gray10","royalblue1","green"))+
  geom_boxplot(outlier.shape=NA)+
  theme(axis.text=element_text(size=10, colour="black"), axis.title=element_text(size=12, colour="black"), strip.text = element_text(size = 12),axis.text.x=element_blank(),legend.position="none")+
  xlab("")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())
PD

#Retain only data of interest
mim_diversity_Shannon <- subset(mim_diversity, select=c("Genotype", "Ecotype","Origin", "Site_planted","shannon_alpha"))
#reshape the data
mim_diversity_Shannon_m<-melt(mim_diversity_Shannon)
mim_diversity_Shannon_m$Genotype <- factor(mim_diversity_Shannon_m$Genotype,levels = c('BulkSoil','LMC','OCC','MRR','SWB'),ordered = TRUE)

#plot the data
Shannon<-ggplot(mim_diversity_Shannon_m, aes(Genotype, value, fill=Ecotype))+
  facet_grid(.~Site_planted)+
  theme_bw()+
  ylab("Shannon Diversity")+
  scale_fill_manual(labels=c("Bulk soil","Coastal ecotype","Inland ecotype"),values=c("gray10","royalblue1","green"))+
  ylim(9.5, 11.5)+
  geom_boxplot(outlier.shape=NA)+
  theme(axis.text=element_text(size=10, colour="black"), axis.title=element_text(size=12, colour="black"), strip.text = element_text(size = 12),legend.position="bottom")+
  xlab("")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  #guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  theme(legend.title=element_blank())
Shannon

library(cowplot)
Fig2<-plot_grid(Rich,PD,Shannon,nrow=3, rel_heights =  c(1.75,1.65,2.5))
Fig2
ggsave("Figure2.tiff",Fig2,width=6,height=6, units="in",dpi=300)



















#First, remove bulk soils from each alpha diversity table
mim_diversity_richness<-mim_diversity_richness[!(mim_diversity_richness$Genotype=="BulkSoil"),]
mim_diversity_PD<-mim_diversity_PD[!(mim_diversity_PD$Genotype=="BulkSoil"),]
mim_diversity_Shannon<-mim_diversity_Shannon[!(mim_diversity_Shannon$Genotype=="BulkSoil"),]





##Statistical tests

#ANOVA for richness


#transform richness
mim_diversity_richness$trans_Richness<- (mim_diversity_richness$Richness)
#Do ANOVA to compare richness
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
rich_fit <- aov(trans_Richness ~ Origin*Site_planted, data=mim_diversity_richness, na.action=na.exclude)
drop1(rich_fit,~.,test="F") # type III SS and F Tests

#TESTING ASSUMPTIONS
#check whether residuals are normally distributed
rich_resids <- residuals(rich_fit)
rich_preds <- predict(rich_fit)
plot(rich_resids ~ rich_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Test for normality
shapiro.test(rich_resids)
qqnorm(rich_resids)
qqline(rich_resids)
library("car")
leveneTest(trans_Richness ~ Origin*Site_planted, data=mim_diversity_richness, na.action=na.exclude)
TukeyHSD(rich_fit)




#transform PD
mim_diversity_PD$trans_PD<- (mim_diversity_PD$PD_whole_tree_alpha)
#Do ANOVA 
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
PD_fit <- aov(trans_PD ~ Origin*Site_planted, data=mim_diversity_PD, na.action=na.exclude)
drop1(PD_fit,~.,test="F") # type III SS and F Tests

#TESTING ASSUMPTIONS
#check whether residuals are normally distributed
PD_resids <- residuals(PD_fit)
PD_preds <- predict(PD_fit)
plot(PD_resids ~ PD_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Test for normality
shapiro.test(PD_resids)
qqnorm(PD_resids)
qqline(PD_resids)
library("car")
leveneTest(trans_PD ~ Origin*Site_planted, data=mim_diversity_PD, na.action=na.exclude)
TukeyHSD(PD_fit)







# NOTE THAT tukey results ARE NOT SIGNIFICANT WHEN OUTLIER IS RETAINED.

#Mim228 appears to be an outlier (compared to the other samples of that ecotype/site).
#Here, we do the calculations to show that it's greater than 3 standard deviations from the mean.
test1<-mim_diversity_Shannon[(mim_diversity_Shannon$Origin=="Inland"),]
test2<-test1[(test1$Site_planted=="Coastal"),]
test3<-test2[(rownames(test2)!="Mim228"),]
mean(test3$shannon_alpha) #10.53
sd(test3$shannon_alpha)  #.302
mean(test3$shannon_alpha) + 3*sd(test3$shannon_alpha) #11.44
mean(test3$shannon_alpha) - 3*sd(test3$shannon_alpha) #9.62


#remove the outlier
mim_diversity_Shannon <- subset(mim_diversity_Shannon, rownames(mim_diversity_Shannon)!="Mim228")


#transform Shannon
mim_diversity_Shannon$trans_Shannon<- (mim_diversity_Shannon$shannon_alpha)
#Do ANOVA 
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Shannon_fit <- aov(trans_Shannon ~ Origin*Site_planted, data=mim_diversity_Shannon, na.action=na.exclude)
drop1(Shannon_fit,~.,test="F") # type III SS and F Tests

#TESTING ASSUMPTIONS
#check whether residuals are normally distributed
Shannon_resids <- residuals(Shannon_fit)
Shannon_preds <- predict(Shannon_fit)
plot(Shannon_resids ~ Shannon_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Test for normality
shapiro.test(Shannon_resids)
qqnorm(Shannon_resids)
qqline(Shannon_resids)
library("car")
leveneTest(trans_Shannon ~ Origin*Site_planted, data=mim_diversity_Shannon, na.action=na.exclude)
TukeyHSD(Shannon_fit)





