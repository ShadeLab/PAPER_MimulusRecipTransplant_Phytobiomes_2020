## Dot plot of class relative abundance for coastal vs inland ecotypes.
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


####
#Prepare data: subset into coastal versus inland sites
#Compare coastal versus inland ecotypes at each site


#First, compare ecotypes at the Coastal site.

#Subset the dataframe for coastal ecotypes (planted coastal)
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Origin == "Coastal" & Site_planted =="Coastal" & Rhizo_Bulk =="Rhizo"))
mim_envir 
coastaleco_coastal <- subset(final_data, select=c("Mim141","Mim143","Mim151","Mim152","Mim161","Mim162","Mim172","Mim176","Mim212"))

#next subset dataframe for inland ecotypes (planted coastal)
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Origin == "Inland" & Site_planted =="Coastal" & Rhizo_Bulk =="Rhizo"))
mim_envir 
inlandeco_coastal <- subset(final_data, select=c("Mim181","Mim182","Mim190","Mim191","Mim224","Mim228","Mim237","Mim238","Mim250","Mim255"))

#make a new dataframe for means and standard deviation
coastalsummary <- data.frame(row.names=rownames(final_data))
coastalsummary$CoastalEco_mean <- rowMeans(coastaleco_coastal)
coastalsummary$CoastalEco_SD<- apply(coastaleco_coastal,1,sd) 
coastalsummary$InlandEco_mean <- rowMeans(inlandeco_coastal)
coastalsummary$InlandEco_SD<- apply(inlandeco_coastal,1,sd) 
#make class names into a new column in that dataframe
coastalsummary <- data.frame(names = row.names(coastalsummary), coastalsummary)

#remove things other than names from class names
coastalsummary$names<-gsub('D_0__Bacteria;', '', coastalsummary$names)
coastalsummary$names<-gsub('D_1__', '', coastalsummary$names)
coastalsummary$names<-gsub(';D_2__', '-', coastalsummary$names)

#reverse order of factors to graph in alphabtical order
coastalsummary$names = forcats::fct_rev(factor(coastalsummary$names))

#T-tests to compare coastal ecotypes versus inland ecotypes at Coastal Site
#Prepare output file into which the t-test results will go.
ttest.out <- data.frame(matrix(nrow = nrow(final_data), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(final_data)
#Perform the t-tests.
for (i in seq(nrow(coastaleco_coastal))){
  test<-c()
  test[[i]]<-t.test(coastaleco_coastal[i,], inlandeco_coastal[i,], paired=FALSE, var.equal = TRUE)
  ttest.out[i,1] <- row.names(final_data)[i]
  ttest.out[i,2] <-(test[[i]]$statistic)
  ttest.out[i,3] <-(test[[i]]$parameter)
  ttest.out[i,4] <-(test[[i]]$p.value)
}
#Calculate adjusted p-values and add to table.
ttest.out[,5] <-p.adjust(ttest.out$p.value, method = "fdr")
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value", "adj.p.value")


#dot plot with means and standard deviation
cols <- c("LINE1"="royalblue1","LINE2"="green")
Fig4_coastal<-ggplot(coastalsummary, aes(names))+
  geom_point(aes(y=CoastalEco_mean,colour="LINE1"),size=3)+
  geom_errorbar(aes(ymax=CoastalEco_mean+CoastalEco_SD, ymin=CoastalEco_mean-CoastalEco_SD, width=0.4))+
  geom_point(aes(y=InlandEco_mean,colour="LINE2"),size=3)+
  geom_errorbar(aes(ymax=InlandEco_mean+InlandEco_SD, ymin=InlandEco_mean-InlandEco_SD, width=0.4))+
  scale_color_manual(name="Ecotype",values=cols,labels=c("Coastal ecotype","Inland ecotype"))+
  ylab("Relative Abundance")+
  ggtitle("Coastal site")+
  theme_bw()+
  coord_flip()+
  xlab("")+
  theme(text = element_text(size=12),legend.position="none")
Fig4_coastal




#Second, compare ecotypes at the Inland site.

#Subset the dataframe for coastal ecotypes (planted inland)
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Origin == "Coastal" & Site_planted =="Inland" & Rhizo_Bulk =="Rhizo"))
mim_envir  
coastaleco_inland <- subset(final_data, select=c("Mim11","Mim12","Mim21","Mim2","Mim31","Mim46","Mim47","Mim57","Mim58","Mim97"))

#next subset dataframe for inland ecotypes (planted inland)
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Origin == "Inland" & Site_planted =="Inland" & Rhizo_Bulk =="Rhizo"))
mim_envir  
inlandeco_inland <- subset(final_data, select=c("Mim109","Mim1","Mim22","Mim32","Mim66","Mim67","Mim83","Mim87","Mim93","Mim98"))

#make a new dataframe for means and standard deviation
inlandsummary <- data.frame(row.names=rownames(final_data))
inlandsummary$CoastalEco_mean <- rowMeans(coastaleco_inland)
inlandsummary$CoastalEco_SD<- apply(coastaleco_inland,1,sd) 
inlandsummary$InlandEco_mean <- rowMeans(inlandeco_inland)
inlandsummary$InlandEco_SD<- apply(inlandeco_inland,1,sd) 
#make class names into a new column in that dataframe
inlandsummary <- data.frame(names = row.names(inlandsummary), inlandsummary)

#remove things other than names from class names
inlandsummary$names<-gsub('D_0__Bacteria;', '', inlandsummary$names)
inlandsummary$names<-gsub('D_1__', '', inlandsummary$names)
inlandsummary$names<-gsub(';D_2__', '-', inlandsummary$names)

#reverse order of factors to graph in alphabtical order
inlandsummary$names = forcats::fct_rev(factor(inlandsummary$names))


###Conduct t-tests: Coastal ecotypes versus inland ecotypes at inland Site
ttest.out <- data.frame(matrix(nrow = nrow(final_data), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(final_data)
for (i in seq(nrow(coastaleco_inland))){
  test<-c()
  test[[i]]<-t.test(coastaleco_inland[i,], inlandeco_inland[i,], paired=FALSE, var.equal = TRUE)
  ttest.out[i,1] <- row.names(final_data)[i]
  ttest.out[i,2] <-(test[[i]]$statistic)
  ttest.out[i,3] <-(test[[i]]$parameter)
  ttest.out[i,4] <-(test[[i]]$p.value)
}
#Calculate adjusted p-values and add to table.
ttest.out[,5] <-p.adjust(ttest.out$p.value, method = "fdr")
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value", "adj.p.value")
#remove things other than names from class names
ttest.out$Taxa<-gsub('D_0__Bacteria;', '', ttest.out$Taxa)
ttest.out$Taxa<-gsub('D_1__', '', ttest.out$Taxa)
ttest.out$Taxa<-gsub(';D_2__', '-', ttest.out$Taxa)
#get list of taxa which significantly differed
sigdifftaxa<-(subset(ttest.out, adj.p.value<=0.05))
sigdifftaxa<-subset(sigdifftaxa,select=c("Taxa"))


#First, set up color palette
cols <- c("LINE1"="royalblue1","LINE2"="green")
#Next, set up vector to plot askterisks for classes that differed between ecotypes.
label.df<-data.frame(names=sigdifftaxa$Taxa, value=c(0.105,0.085,0.022,0.013,0.004))

Fig4_inland<-ggplot(inlandsummary, aes(names))+
  geom_point(aes(y=CoastalEco_mean, colour="LINE1"),size=3)+
  geom_errorbar(aes(ymax=CoastalEco_mean+CoastalEco_SD, ymin=CoastalEco_mean-CoastalEco_SD, width=0.4))+
  geom_point(aes(y=InlandEco_mean,colour="LINE2"),size=3)+
  geom_errorbar(aes(ymax=InlandEco_mean+InlandEco_SD, ymin=InlandEco_mean-InlandEco_SD, width=0.4))+
  scale_color_manual(name="",values=cols,labels=c("Coastal ecotype","Inland ecotype"))+
  ylab("Relative Abundance")+
  ggtitle("Inland site")+
  theme_bw()+
  coord_flip()+
  xlab("")+
  theme(axis.text.y=element_blank())+
  theme(text = element_text(size=12),legend.position="right")+
  geom_text(data=label.df, aes(x=label.df$names,y=label.df$value,label = "*"),size=6,fontface="bold",position=position_nudge(x=0.15))
Fig4_inland




#Plotting the two together
library("grid")
tiff("Figure4.tiff", height = 6, width = 14,units="in",res=300)
grid.newpage()
grid.draw(cbind(ggplotGrob(Fig4_coastal), ggplotGrob(Fig4_inland), size="last"))
dev.off()
