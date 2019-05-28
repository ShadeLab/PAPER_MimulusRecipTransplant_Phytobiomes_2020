## Dot plot of class relative abundance between rhizosphere and bulk soil.

library(ggplot2)
library(plyr)
library(reshape2)

#read in class relative abundance data
mim_class <- read.delim("single_rare_L3.txt", row.names=1)

#find average relative abundance for each class (across all samples)
average_relabund <- as.data.frame(rowMeans(mim_class))
colnames(average_relabund)<-c("average_relabund")
mim_class <-cbind(mim_class, average_relabund)

#sort classes by average relative abundance
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
#Prepare data to plot Rhizosphere versus Bulk soil samples at Coastal site.

#first subset the full dataframe for just coastal site coastal ecotype samples.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted =="Coastal" & Origin=="Coastal" & Rhizo_Bulk =="Rhizo"))
mim_envir  
coastalrhizo_coastaleco <- subset(final_data, select=c("Mim141", "Mim143", "Mim151", "Mim152", "Mim161", "Mim162", "Mim172", "Mim176", "Mim212"))

#next subset the full dataframe for just coastal site inland ecotype samples.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted =="Coastal" & Origin=="Inland" & Rhizo_Bulk =="Rhizo"))
mim_envir  
coastalrhizo_inlandeco <- subset(final_data, select=c("Mim181", "Mim182", "Mim190", "Mim191", "Mim224", "Mim228", "Mim237", "Mim238", "Mim250", "Mim255"))

#lastly, subset the full dataframe for just coastal site bulk soils.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Origin =="Coastal" & Rhizo_Bulk =="Bulk"))
mim_envir 
coastalbulk <- subset(final_data, select=c("Mim261","Mim266","Mim271","Mim277","Mim286"))
  
#make a new dataframe for means and standard deviations
coastalsummary <- data.frame(row.names=rownames(final_data))
#calculate mean and SD for rhizosphere samples and for bulk samples
coastalsummary$CoastalEco_mean <- rowMeans(coastalrhizo_coastaleco)
coastalsummary$CoastalEco_SD<- apply(coastalrhizo_coastaleco,1,sd) 
coastalsummary$InlandEco_mean <- rowMeans(coastalrhizo_inlandeco)
coastalsummary$InlandEco_SD<- apply(coastalrhizo_inlandeco,1,sd) 
coastalsummary$Bulk_mean <- rowMeans(coastalbulk)
coastalsummary$Bulk_SD<- apply(coastalbulk,1,sd) 
#make class names into a new column in that dataframe
coastalsummary <- data.frame(names = row.names(coastalsummary), coastalsummary)

#remove things other than names from class names
coastalsummary$names<-gsub('D_0__Bacteria;', '', coastalsummary$names)
coastalsummary$names<-gsub('D_1__', '', coastalsummary$names)
coastalsummary$names<-gsub(';D_2__', '-', coastalsummary$names)

#reverse order of factors to graph in alphabtical order
coastalsummary$names = forcats::fct_rev(factor(coastalsummary$names))



####
#Prepare data to plot Rhizosphere versus Bulk soil samples at Inland site.

#First, subset the full dataframe for just inland site coastal ecotype samples.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted =="Inland" & Origin=="Coastal" & Rhizo_Bulk =="Rhizo"))
mim_envir 
inlandrhizo_coastaleco <- subset(final_data, select=c("Mim11", "Mim12", "Mim21", "Mim2",  "Mim31", "Mim46", "Mim47", "Mim57", "Mim58", "Mim97"))

#next, subset the full dataframe for just inland site inland ecotype samples.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted =="Inland" & Origin=="Inland" & Rhizo_Bulk =="Rhizo"))
mim_envir 
inlandrhizo_inlandeco <- subset(final_data, select=c( "Mim109", "Mim1",   "Mim22",  "Mim32",  "Mim66",  "Mim67",  "Mim83",  "Mim87",  "Mim93",  "Mim98"))

#lastly subset the full dataframe for just inland site bulk soils.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Origin =="Inland" & Rhizo_Bulk =="Bulk"))
mim_envir 
inlandbulk <- subset(final_data, select=c("Mim116","Mim124","Mim129","Mim132","Mim139"))

#make a new dataframe for means and standard deviations
inlandsummary <- data.frame(row.names=rownames(final_data))
inlandsummary$CoastalEco_mean <- rowMeans(inlandrhizo_coastaleco)
inlandsummary$CoastalEco_SD<- apply(inlandrhizo_coastaleco,1,sd)
inlandsummary$InlandEco_mean <- rowMeans(inlandrhizo_inlandeco)
inlandsummary$InlandEco_SD<- apply(inlandrhizo_inlandeco,1,sd)
inlandsummary$Bulk_mean <- rowMeans(inlandbulk)
inlandsummary$Bulk_SD<- apply(inlandbulk,1,sd) 
#make class names into a new column in that dataframe
inlandsummary <- data.frame(names = row.names(inlandsummary), inlandsummary)

#remove things other than names from class names
inlandsummary$names<-gsub('D_0__Bacteria;', '', inlandsummary$names)
inlandsummary$names<-gsub('D_1__', '', inlandsummary$names)
inlandsummary$names<-gsub(';D_2__', '-', inlandsummary$names)

#reverse order of factors to graph in alphabtical order
inlandsummary$names = forcats::fct_rev(factor(inlandsummary$names))



####
#Perform T-tests to compare rhizosphere versus bulk samples

#First, at the Coastal Site:  coastal ecotypes vs bulk
#Prepare output file into which the t-test results will go.
ttest.out <- data.frame(matrix(nrow = nrow(final_data), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(final_data)
#Perform the t-tests.
for(i in seq(nrow(coastalrhizo_coastaleco))){
  test<-c()
  test[[i]]<-t.test(coastalrhizo_coastaleco[i,], coastalbulk[i,], paired=FALSE, var.equal = TRUE)
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
sigdifftaxa_coastalsite_coastaleco<-(subset(ttest.out, adj.p.value<=0.05))
sigdifftaxa_coastalsite_coastaleco<-subset(sigdifftaxa_coastalsite_coastaleco,select=c("Taxa"))



#Second, at the Coastal Site:  inland ecotypes vs bulk
#Prepare output file into which the t-test results will go.
ttest.out <- data.frame(matrix(nrow = nrow(final_data), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(final_data)
#Perform the t-tests.
for(i in seq(nrow(coastalrhizo_inlandeco))){
  test<-c()
  test[[i]]<-t.test(coastalrhizo_inlandeco[i,], coastalbulk[i,], paired=FALSE, var.equal = TRUE)
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
sigdifftaxa_coastalsite_inlandeco<-(subset(ttest.out, adj.p.value<=0.05))
sigdifftaxa_coastalsite_inlandeco<-subset(sigdifftaxa_coastalsite_inlandeco,select=c("Taxa"))




#Third, at the Inland Site:  coastal ecotypes vs bulk
#Prepare output file into which the t-test results will go.
ttest.out <- data.frame(matrix(nrow = nrow(final_data), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(final_data)
#Perform the t-tests.
for(i in seq(nrow(inlandrhizo_coastaleco))){
  test<-c()
  test[[i]]<-t.test(inlandrhizo_coastaleco[i,], inlandbulk[i,], paired=FALSE, var.equal = TRUE)
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
sigdifftaxa_inlandsite_coastaleco<-(subset(ttest.out, adj.p.value<=0.05))
sigdifftaxa_inlandsite_coastaleco<-subset(sigdifftaxa_inlandsite_coastaleco,select=c("Taxa"))



#Fourth, at the Inland Site:  inland ecotypes vs bulk
#Prepare output file into which the t-test results will go.
ttest.out <- data.frame(matrix(nrow = nrow(final_data), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(final_data)
#Perform the t-tests.
for(i in seq(nrow(inlandrhizo_inlandeco))){
  test<-c()
  test[[i]]<-t.test(inlandrhizo_inlandeco[i,], inlandbulk[i,], paired=FALSE, var.equal = TRUE)
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
sigdifftaxa_inlandsite_inlandeco<-(subset(ttest.out, adj.p.value<=0.05))
sigdifftaxa_inlandsite_inlandeco<-subset(sigdifftaxa_inlandsite_inlandeco,select=c("Taxa"))








####
#Plotting 


#COASTAL DATA. 

#First, set up color palette
cols <- c("LINE1"="royalblue1","LINE2"="green", "LINE3"="gray10")
#Next, set up vector to plot askterisks for classes that differed between rhizo and bulk.
label_coastalsite_coastaleco.df<-merge(x=sigdifftaxa_coastalsite_coastaleco,y=coastalsummary[,c("names","CoastalEco_mean")],by.x='Taxa',by.y='names',all.x=FALSE,all.y=FALSE)
label_coastalsite_inlandeco.df<-merge(x=sigdifftaxa_coastalsite_inlandeco, y=coastalsummary[,c("names","InlandEco_mean")],by.x='Taxa',by.y='names',all.x=FALSE,all.y=FALSE)

#dot plot with means and standard deviation
Fig3_coastal<-ggplot(coastalsummary, aes(names))+
  geom_point(aes(y=CoastalEco_mean,colour="LINE1"), size=3)+
  geom_errorbar(aes(ymax=CoastalEco_mean+CoastalEco_SD, ymin=CoastalEco_mean-CoastalEco_SD, width=0.4))+
  geom_point(aes(y=InlandEco_mean,colour="LINE2"), size=3)+
  geom_errorbar(aes(ymax=InlandEco_mean+InlandEco_SD, ymin=InlandEco_mean-InlandEco_SD, width=0.4))+
  geom_point(aes(y=Bulk_mean,colour="LINE3"),size=3)+
  geom_errorbar(aes(ymax=Bulk_mean+Bulk_SD, ymin=Bulk_mean-Bulk_SD, width=0.4))+
  scale_color_manual(name="Type",values=cols,labels=c("Coastal ecotype","Inland ecotype","Bulk soil"))+
  ylab("Relative Abundance")+
  ggtitle("Coastal site")+
  theme_bw()+
  coord_flip()+
  xlab("")+
  theme(text = element_text(size=12),legend.position="none")+
  geom_text(data=label_coastalsite_coastaleco.df, aes(x=label_coastalsite_coastaleco.df$Taxa,y=label_coastalsite_coastaleco.df$CoastalEco_mean,label = "*"),size=6,fontface="bold",position=position_nudge(x=0.25))+
  geom_text(data=label_coastalsite_inlandeco.df, aes(x=label_coastalsite_inlandeco.df$Taxa,y=label_coastalsite_inlandeco.df$InlandEco_mean,label = "*"),size=6,fontface="bold",position=position_nudge(x=0.25))
Fig3_coastal





#INLAND DATA. 

#First, set up color palette
cols <- c("LINE1"="royalblue1","LINE2"="green", "LINE3"="gray10")

#Next, set up vector to plot askterisks for classes that differed between rhizo and bulk.
label_inlandsite_coastaleco.df<-merge(x=sigdifftaxa_inlandsite_coastaleco,y=inlandsummary[,c("names","CoastalEco_mean")],by.x='Taxa',by.y='names',all.x=FALSE,all.y=FALSE)
label_inlandsite_inlandeco.df<-merge(x=sigdifftaxa_inlandsite_inlandeco, y=inlandsummary[,c("names","InlandEco_mean")],by.x='Taxa',by.y='names',all.x=FALSE,all.y=FALSE)

#dot plot with means and standard deviation
Fig3_inland<-ggplot(inlandsummary, aes(names))+
  geom_point(aes(y=CoastalEco_mean,colour="LINE1"), size=3)+
  geom_errorbar(aes(ymax=CoastalEco_mean+CoastalEco_SD, ymin=CoastalEco_mean-CoastalEco_SD, width=0.4))+
  geom_point(aes(y=InlandEco_mean,colour="LINE2"), size=3)+
  geom_errorbar(aes(ymax=InlandEco_mean+InlandEco_SD, ymin=InlandEco_mean-InlandEco_SD, width=0.4))+
  geom_point(aes(y=Bulk_mean,colour="LINE3"),size=3)+
  geom_errorbar(aes(ymax=Bulk_mean+Bulk_SD, ymin=Bulk_mean-Bulk_SD, width=0.4))+
  scale_color_manual(name="",values=cols,labels=c("Coastal ecotype","Inland ecotype","Bulk soil"))+
  ylab("Relative Abundance")+
  ggtitle("Inland site")+
  theme_bw()+
  coord_flip()+
  xlab("")+
  theme(axis.text.y=element_blank())+
  theme(text = element_text(size=12),legend.position="right")+
  geom_text(data=label_inlandsite_coastaleco.df,aes(x=label_inlandsite_coastaleco.df$Taxa,y=label_inlandsite_coastaleco.df$CoastalEco_mean,label = "*"),size=6,fontface="bold",position=position_nudge(x=0.25))+
  geom_text(data=label_inlandsite_inlandeco.df, aes(x=label_inlandsite_inlandeco.df$Taxa,y=label_inlandsite_inlandeco.df$InlandEco_mean,label = "*"),size=6,fontface="bold",position=position_nudge(x=0.25))
Fig3_inland




#Plotting the two together
library("grid")
tiff("Figure3.tiff", height = 6, width = 14,units="in",res=300)
grid.newpage()
grid.draw(cbind(ggplotGrob(Fig3_coastal), ggplotGrob(Fig3_inland), size='last'))
dev.off()


