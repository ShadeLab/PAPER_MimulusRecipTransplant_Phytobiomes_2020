## Dot plot of class relative abundance for each genotype.
library(ggplot2)
library(plyr)
library(reshape2)

#read in class relative abundance data
mim_class <- read.delim("single_rare_L3.txt", row.names=1)

#find average relative abundance for each class (across samples)
average_relabund <- as.data.frame(rowMeans(mim_class))
colnames(average_relabund)<-c("average_relabund")
mim_class <-cbind(mim_class, average_relabund)

#sort classes by relative abundance
mim_class_sorted <- mim_class[order(-average_relabund),] 

#Select top20 most abundant classes.
mim_class_top20<-mim_class_sorted[1:20,]

#sum all the other (less abundant) classes for each sample.
mim_class_others<-mim_class_sorted[21:153,]
others_colsums <- as.data.frame(colSums(mim_class_others))
colnames(others_colsums)<-c("Less Abundant Classes (Pooled)")
#transpose the data frame to make in same orientation as main table
others_colsums<-as.data.frame(t(others_colsums))

#put the two tables together (20 most abundant classes, plus column of all other pooled classes)
final_data <- rbind(mim_class_top20,others_colsums)
#remove column of average relative abundance data.
final_data <-subset(final_data, select=-c(average_relabund))




#panel A = at coastal site, coastal genotype1 vs coastal genotype2

#first, subset dataframe for just coastal genotype SWB planted at coastal site.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted=="Coastal" & Origin == "Coastal" & Genotype=="SWB"))
mim_envir  
coastalSWB_coastalsite <- subset(final_data, select=c("Mim151","Mim152","Mim172","Mim176"))

#next subset dataframe for just coastal genotype MRR planted at coastal site.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted=="Coastal" & Origin == "Coastal" & Genotype=="MRR"))
mim_envir  
coastalMRR_coastalsite <- subset(final_data, select=c("Mim141","Mim143","Mim161","Mim162","Mim212"))

#make a new dataframe for means and standard deviation
coastalgeno_coastalsitesummary <- data.frame(row.names=rownames(coastalSWB_coastalsite))
coastalgeno_coastalsitesummary$coastalSWB_coastalsite_mean <- rowMeans(coastalSWB_coastalsite)
coastalgeno_coastalsitesummary$coastalSWB_coastalsite_SD<- apply(coastalSWB_coastalsite,1,sd) 
coastalgeno_coastalsitesummary$coastalMRR_coastalsite_mean <- rowMeans(coastalMRR_coastalsite)
coastalgeno_coastalsitesummary$coastalMRR_coastalsite_SD<- apply(coastalMRR_coastalsite,1,sd) 
#make class names into a new column in that dataframe
coastalgeno_coastalsitesummary <- data.frame(names = row.names(coastalgeno_coastalsitesummary), coastalgeno_coastalsitesummary)

#remove things other than names from class names
coastalgeno_coastalsitesummary$names<-gsub('D_0__Bacteria;', '', coastalgeno_coastalsitesummary$names)
coastalgeno_coastalsitesummary$names<-gsub('D_1__', '', coastalgeno_coastalsitesummary$names)
coastalgeno_coastalsitesummary$names<-gsub(';D_2__', '-', coastalgeno_coastalsitesummary$names)

#reverse order of factors to graph in alphabetical order
coastalgeno_coastalsitesummary$names = forcats::fct_rev(factor(coastalgeno_coastalsitesummary$names))


#T-tests comparing Coastal genotype SWB vs coastal genotype MRR at Coastal site
ttest.out <- data.frame(matrix(nrow = nrow(final_data), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(final_data)
for (i in seq(nrow(coastalSWB_coastalsite))){
  test<-c()
  test[[i]]<-t.test(coastalSWB_coastalsite[i,], coastalMRR_coastalsite[i,], paired=FALSE, var.equal = TRUE)
  ttest.out[i,1] <- row.names(final_data)[i]
  ttest.out[i,2] <-(test[[i]]$statistic)
  ttest.out[i,3] <-(test[[i]]$parameter)
  ttest.out[i,4] <-(test[[i]]$p.value)
}
#Calculate adjusted p-values and add to table.
ttest.out[,5] <-p.adjust(ttest.out$p.value, method = "fdr")
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value", "adj.p.value")

#dot plot with means and standard deviations
cols <- c("LINE1"="green","LINE2"="green4")
Fig5_coastalgeno_coastalsite<-ggplot(coastalgeno_coastalsitesummary, aes(names))+
  geom_point(aes(y=coastalSWB_coastalsite_mean,colour="LINE1"), size=3)+
  geom_errorbar(aes(ymax=coastalSWB_coastalsite_mean+coastalSWB_coastalsite_SD, ymin=coastalSWB_coastalsite_mean-coastalSWB_coastalsite_SD, width=0.4))+
  geom_point(aes(y=coastalMRR_coastalsite_mean,colour="LINE2"),size=3)+
  geom_errorbar(aes(ymax=coastalMRR_coastalsite_mean+coastalMRR_coastalsite_SD, ymin=coastalMRR_coastalsite_mean-coastalMRR_coastalsite_SD, width=0.4))+
  scale_color_manual(name="Coastal Genotype",values=cols,labels=c("SWB","MRR"))+
  ylab("")+
  ggtitle("Coastal Site")+
  theme_bw()+
  coord_flip(ylim=c(0.00,0.1999))+
  xlab("")+
  theme(axis.text.x=element_blank())+
  theme(text = element_text(size=12),legend.position="none")
Fig5_coastalgeno_coastalsite



#Panel B = at inland site, coastal genotype1 vs coastal genotype2

#first, subset dataframe for just coastal genotype SWB planted at inland site.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted=="Inland" & Origin == "Coastal" & Genotype=="SWB"))
mim_envir 
coastalSWB_inlandsite <- subset(final_data, select=c("Mim12","Mim31","Mim47","Mim57","Mim97"))

#next subset dataframe for just coastal genotype MRR planted at inland site.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted=="Inland" & Origin == "Coastal" & Genotype=="MRR"))
mim_envir  
coastalMRR_inlandsite <- subset(final_data, select=c("Mim11","Mim21","Mim2","Mim46","Mim58"))

#make a new dataframe for means and standard deviation
coastalgeno_inlandsitesummary <- data.frame(row.names=rownames(coastalSWB_inlandsite))
coastalgeno_inlandsitesummary$coastalSWB_inlandsite_mean <- rowMeans(coastalSWB_inlandsite)
coastalgeno_inlandsitesummary$coastalSWB_inlandsite_SD<- apply(coastalSWB_inlandsite,1,sd) 
coastalgeno_inlandsitesummary$coastalMRR_inlandsite_mean <- rowMeans(coastalMRR_inlandsite)
coastalgeno_inlandsitesummary$coastalMRR_inlandsite_SD<- apply(coastalMRR_inlandsite,1,sd) 
#make class names into a new column in that dataframe
coastalgeno_inlandsitesummary <- data.frame(names = row.names(coastalgeno_inlandsitesummary), coastalgeno_inlandsitesummary)

#remove things other than names from class names
coastalgeno_inlandsitesummary$names<-gsub('D_0__Bacteria;', '', coastalgeno_inlandsitesummary$names)
coastalgeno_inlandsitesummary$names<-gsub('D_1__', '', coastalgeno_inlandsitesummary$names)
coastalgeno_inlandsitesummary$names<-gsub(';D_2__', '-', coastalgeno_inlandsitesummary$names)

#reverse order of factors to graph in alphabetical order
coastalgeno_inlandsitesummary$names = forcats::fct_rev(factor(coastalgeno_inlandsitesummary$names))

#T-tests comparing Coastal genotype SWB vs coastal genotype MRR at inland site
ttest.out <- data.frame(matrix(nrow = nrow(final_data), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(final_data)
for (i in seq(nrow(coastalSWB_inlandsite))){
  test<-c()
  test[[i]]<-t.test(coastalSWB_inlandsite[i,], coastalMRR_inlandsite[i,], paired=FALSE, var.equal = TRUE)
  ttest.out[i,1] <- row.names(final_data)[i]
  ttest.out[i,2] <-(test[[i]]$statistic)
  ttest.out[i,3] <-(test[[i]]$parameter)
  ttest.out[i,4] <-(test[[i]]$p.value)
}
#Calculate adjusted p-values and add to table.
ttest.out[,5] <-p.adjust(ttest.out$p.value, method = "fdr")
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value", "adj.p.value")

cols <- c("LINE1"="green","LINE2"="green4")
Fig5_coastalgeno_inlandsite<-ggplot(coastalgeno_inlandsitesummary, aes(names))+
  geom_point(aes(y=coastalSWB_inlandsite_mean,colour="LINE1"), size=3)+
  geom_errorbar(aes(ymax=coastalSWB_inlandsite_mean+coastalSWB_inlandsite_SD, ymin=coastalSWB_inlandsite_mean-coastalSWB_inlandsite_SD, width=0.4))+
  geom_point(aes(y=coastalMRR_inlandsite_mean,colour="LINE2"),size=3)+
  geom_errorbar(aes(ymax=coastalMRR_inlandsite_mean+coastalMRR_inlandsite_SD, ymin=coastalMRR_inlandsite_mean-coastalMRR_inlandsite_SD, width=0.4))+
  scale_color_manual(name="Coastal Genotype",values=cols,labels=c("SWB","MRR"))+
  ylab("")+
  ggtitle("Inland Site")+
  theme_bw()+
  coord_flip(ylim=c(0.00,0.1999))+
  xlab("")+
  theme(axis.text.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(text = element_text(size=12),legend.position="right")
Fig5_coastalgeno_inlandsite



#PanelC = at coastal site, inland genotype1 vs inland genotype2

#first, subset dataframe for just inland genotype LMC planted at coastal site.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted=="Coastal" & Origin == "Inland" & Genotype=="LMC"))
mim_envir 
inlandLMC_coastalsite <- subset(final_data, select=c("Mim191","Mim237","Mim238","Mim255"))

#next subset dataframe for just inland genotype OCC planted at coastal site.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted=="Coastal" & Origin == "Inland" & Genotype=="OCC"))
mim_envir  
inlandOCC_coastalsite <- subset(final_data, select=c("Mim181","Mim182","Mim190","Mim224","Mim228","Mim250"))

#make a new dataframe for means and standard deviation
inlandgeno_coastalsitesummary <- data.frame(row.names=rownames(inlandLMC_coastalsite))
inlandgeno_coastalsitesummary$inlandLMC_coastalsite_mean <- rowMeans(inlandLMC_coastalsite)
inlandgeno_coastalsitesummary$inlandLMC_coastalsite_SD<- apply(inlandLMC_coastalsite,1,sd) 
inlandgeno_coastalsitesummary$inlandOCC_coastalsite_mean <- rowMeans(inlandOCC_coastalsite)
inlandgeno_coastalsitesummary$inlandOCC_coastalsite_SD<- apply(inlandOCC_coastalsite,1,sd) 
#make class names into a new column in that dataframe
inlandgeno_coastalsitesummary <- data.frame(names = row.names(inlandgeno_coastalsitesummary), inlandgeno_coastalsitesummary)

#remove things other than names from class names
inlandgeno_coastalsitesummary$names<-gsub('D_0__Bacteria;', '', inlandgeno_coastalsitesummary$names)
inlandgeno_coastalsitesummary$names<-gsub('D_1__', '', inlandgeno_coastalsitesummary$names)
inlandgeno_coastalsitesummary$names<-gsub(';D_2__', '-', inlandgeno_coastalsitesummary$names)

#reverse order of factors to graph in alphabtical order
inlandgeno_coastalsitesummary$names = forcats::fct_rev(factor(inlandgeno_coastalsitesummary$names))

#T-tests comparing Inland genotype LMC vs inland genotype OCC at Coastal site
ttest.out <- data.frame(matrix(nrow = nrow(final_data), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(final_data)
for (i in seq(nrow(inlandLMC_coastalsite))){
  test<-c()
  test[[i]]<-t.test(inlandLMC_coastalsite[i,], inlandOCC_coastalsite[i,], paired=FALSE, var.equal = TRUE)
  ttest.out[i,1] <- row.names(final_data)[i]
  ttest.out[i,2] <-(test[[i]]$statistic)
  ttest.out[i,3] <-(test[[i]]$parameter)
  ttest.out[i,4] <-(test[[i]]$p.value)
}
#Calculate adjusted p-values and add to table.
ttest.out[,5] <-p.adjust(ttest.out$p.value, method = "fdr")
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value", "adj.p.value")

cols <- c("LINE1"="royalblue1","LINE2"="royalblue4")
Fig5_inlandgeno_coastalsite<-ggplot(inlandgeno_coastalsitesummary, aes(names))+
  geom_point(aes(y=inlandLMC_coastalsite_mean,colour="LINE1"), size=3)+
  geom_errorbar(aes(ymax=inlandLMC_coastalsite_mean+inlandLMC_coastalsite_SD, ymin=inlandLMC_coastalsite_mean-inlandLMC_coastalsite_SD, width=0.4))+
  geom_point(aes(y=inlandOCC_coastalsite_mean,colour="LINE2"),size=3)+
  geom_errorbar(aes(ymax=inlandOCC_coastalsite_mean+inlandOCC_coastalsite_SD, ymin=inlandOCC_coastalsite_mean-inlandOCC_coastalsite_SD, width=0.4))+
  scale_color_manual(name="Inland Genotype",values=cols,labels=c("LMC","OCC"))+
  ylab("Relative Abundance")+
  ggtitle("")+
  theme_bw()+
  coord_flip(ylim=c(0.00,0.1999))+
  xlab("")+
  theme(text = element_text(size=12),legend.position="none")
Fig5_inlandgeno_coastalsite



#PanelD = at inland site, inland genotype1 vs inland genotype2

#first, subset dataframe for just inland genotype LMC planted at inland site.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted=="Inland" & Origin == "Inland" & Genotype=="LMC"))
mim_envir 
inlandLMC_inlandsite <- subset(final_data, select=c("Mim109","Mim22","Mim83","Mim87","Mim93"))

#next subset dataframe for just inland genotype OCC planted at inland site.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted=="Inland" & Origin == "Inland" & Genotype=="OCC"))
mim_envir  
inlandOCC_inlandsite <- subset(final_data, select=c("Mim1","Mim32","Mim66","Mim67","Mim98"))

#make a new dataframe for means and standard deviation
inlandgeno_inlandsitesummary <- data.frame(row.names=rownames(inlandLMC_inlandsite))
inlandgeno_inlandsitesummary$inlandLMC_inlandsite_mean <- rowMeans(inlandLMC_inlandsite)
inlandgeno_inlandsitesummary$inlandLMC_inlandsite_SD<- apply(inlandLMC_inlandsite,1,sd) 
inlandgeno_inlandsitesummary$inlandOCC_inlandsite_mean <- rowMeans(inlandOCC_inlandsite)
inlandgeno_inlandsitesummary$inlandOCC_inlandsite_SD<- apply(inlandOCC_inlandsite,1,sd) 
#make class names into a new column in that dataframe
inlandgeno_inlandsitesummary <- data.frame(names = row.names(inlandgeno_inlandsitesummary), inlandgeno_inlandsitesummary)

#remove things other than names from class names
inlandgeno_inlandsitesummary$names<-gsub('D_0__Bacteria;', '', inlandgeno_inlandsitesummary$names)
inlandgeno_inlandsitesummary$names<-gsub('D_1__', '', inlandgeno_inlandsitesummary$names)
inlandgeno_inlandsitesummary$names<-gsub(';D_2__', '-', inlandgeno_inlandsitesummary$names)

#reverse order of factors to graph in alphabetical order
inlandgeno_inlandsitesummary$names = forcats::fct_rev(factor(inlandgeno_inlandsitesummary$names))

#T-tests comparing Inland genotype LMC vs inland genotype OCC at inland site
ttest.out <- data.frame(matrix(nrow = nrow(final_data), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(final_data)
for (i in seq(nrow(inlandLMC_inlandsite))){
  test<-c()
  test[[i]]<-t.test(inlandLMC_inlandsite[i,], inlandOCC_inlandsite[i,], paired=FALSE, var.equal = TRUE)
  ttest.out[i,1] <- row.names(final_data)[i]
  ttest.out[i,2] <-(test[[i]]$statistic)
  ttest.out[i,3] <-(test[[i]]$parameter)
  ttest.out[i,4] <-(test[[i]]$p.value)
}
#Calculate adjusted p-values and add to table. 
ttest.out[,5] <-p.adjust(ttest.out$p.value, method = "fdr")
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value", "adj.p.value")

cols <- c("LINE1"="royalblue1","LINE2"="royalblue4")
Fig5_inlandgeno_inlandsite<-ggplot(inlandgeno_inlandsitesummary, aes(names))+
  geom_point(aes(y=inlandLMC_inlandsite_mean,colour="LINE1"), size=3)+
  geom_errorbar(aes(ymax=inlandLMC_inlandsite_mean+inlandLMC_inlandsite_SD, ymin=inlandLMC_inlandsite_mean-inlandLMC_inlandsite_SD, width=0.4))+
  geom_point(aes(y=inlandOCC_inlandsite_mean,colour="LINE2"),size=3)+
  geom_errorbar(aes(ymax=inlandOCC_inlandsite_mean+inlandOCC_inlandsite_SD, ymin=inlandOCC_inlandsite_mean-inlandOCC_inlandsite_SD, width=0.4))+
  scale_color_manual(name="Inland Genotype",values=cols,labels=c("LMC","OCC"))+
  ylab("Relative Abundance")+
  ggtitle("")+
  theme_bw()+
  coord_flip(ylim=c(0.00,0.1999))+
  xlab("")+
  theme(axis.text.y=element_blank())+
  theme(text = element_text(size=12),legend.position="right")
Fig5_inlandgeno_inlandsite




#Plotting the graphs together
library("gridExtra")
Figure5 <- grid.arrange(Fig5_coastalgeno_coastalsite,Fig5_coastalgeno_inlandsite,Fig5_inlandgeno_coastalsite,Fig5_inlandgeno_inlandsite, ncol = 2,nrow=2, widths=c(3.3,2.8),heights = c(3, 3))
Figure5
ggsave("Figure5.pdf",Figure5,width=300,height=174, units="mm")



