
#determining rare vs abundant taxa for coastal ecotype
library(vegan)
library(ggplot2)

#import otu table (remove "#Constructed from biom file " from first line before importing).
mim_otu<-read.table("zotu_table_forR.txt", row.names=1, header=T)

#subset data to include only Coastal Ecotype at inland site
mim_meta<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_meta<-rownames(subset(mim_meta, Site_planted =="Inland" & Origin=="Coastal" & Rhizo_Bulk =="Rhizo"))
mim_meta 
coastaleco_otu_inland <- subset(mim_otu, select=colnames(mim_otu) %in% mim_meta)
#remove OTUs with zero abundance
coastaleco_otu_inland <- coastaleco_otu_inland[ rowSums(coastaleco_otu_inland)!=0, ]

#make new dataframe housing the relative abundances of OTUs in the original dataset.
OTU_relabund<- as.data.frame(rowSums(coastaleco_otu_inland)/sum(rowSums(coastaleco_otu_inland)))
row.names(OTU_relabund)<- rownames(coastaleco_otu_inland)
colnames(OTU_relabund)<-c('RelAbund')

#set up loop:

#function
Rare_vs_abundant_calc<-function(threshold){
  #first subset list of OTU names to exclude those below given threshold.
  Rare_OTUs<-rownames(subset(OTU_relabund, RelAbund >= threshold))
  #then remove those OTUs from the OTU table to get the truncated dataset.
  truncated_mim_otu<-subset(coastaleco_otu_inland, rownames(coastaleco_otu_inland) %in% Rare_OTUs)
  #transpose truncated OTU table
  truncated_mim_otu_t<-t(truncated_mim_otu)
  #calculate Bray-Curtis dissimilarity matrix for truncated OTU table
  truncated_BC <- vegdist(truncated_mim_otu_t, method="bray")
  #transpose original OTU table
  coastaleco_otu_inland_t<-t(coastaleco_otu_inland)
  #calculate Bray-Curtis dissimilarity matrix for rel abund table
  original_BC <- vegdist(coastaleco_otu_inland_t, method="bray")
  #calculate correlation between original and truncated dissimilarity matrices
  mantel_output<-mantel(original_BC, truncated_BC, permutations = 999, method='pearson')
  #output the pvalue
  mantel_results <- as.data.frame(cbind(mantel_output$statistic,mantel_output$signif))
  
}


#calculate pairwise dissim matrix for truncated dataset (based on relative abundances)
Thresholds<- seq(0.000001, 0.009, 0.0001)

results<-t(as.data.frame(sapply(Thresholds,Rare_vs_abundant_calc)))
Thresholds<-as.data.frame(Thresholds)
Thresh_Pvalue<-as.data.frame(cbind(Thresholds,results))
colnames(Thresh_Pvalue)<-c("Threshold","Statistic","Pvalue")
Thresh_Pvalue$AdjP <-p.adjust(Thresh_Pvalue$Pvalue, method = "fdr")
Thresh_Fig<- ggplot(Thresh_Pvalue, aes(Threshold, AdjP))+
  geom_line()+
  geom_point()+
  ylab("Statistic")+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle("(A)")+
  xlab("Threshold")
Thresh_Fig











#determining rare vs abundant taxa for Inland ecotype
#import otu table (remove "#Constructed from biom file " from first line before importing).
mim_otu<-read.table("zotu_table_forR.txt", row.names=1, header=T)

#subset data to include only Coastal Ecotype at inland site
mim_meta<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_meta<-rownames(subset(mim_meta, Site_planted =="Inland" & Origin=="Inland" & Rhizo_Bulk =="Rhizo"))
mim_meta 
inlandeco_otu_inland <- subset(mim_otu, select=colnames(mim_otu) %in% mim_meta)
#remove OTUs with zero abundance
inlandeco_otu_inland <- inlandeco_otu_inland[ rowSums(inlandeco_otu_inland)!=0, ]

#make new dataframe housing the relative abundances of OTUs in the original dataset.
OTU_relabund<- as.data.frame(rowSums(inlandeco_otu_inland)/sum(rowSums(inlandeco_otu_inland)))
row.names(OTU_relabund)<- rownames(inlandeco_otu_inland)
colnames(OTU_relabund)<-c('RelAbund')

#set up loop:

#function
Rare_vs_abundant_calc<-function(threshold){
  #first subset list of OTU names to exclude those below given threshold.
  Rare_OTUs<-rownames(subset(OTU_relabund, RelAbund >= threshold))
  #then remove those OTUs from the OTU table to get the truncated dataset.
  truncated_mim_otu<-subset(inlandeco_otu_inland, rownames(inlandeco_otu_inland) %in% Rare_OTUs)
  #transpose truncated OTU table
  truncated_mim_otu_t<-t(truncated_mim_otu)
  #calculate Bray-Curtis dissimilarity matrix for truncated OTU table
  truncated_BC <- vegdist(truncated_mim_otu_t, method="bray")
  #transpose original OTU table
  inlandeco_otu_inland_t<-t(inlandeco_otu_inland)
  #calculate Bray-Curtis dissimilarity matrix for rel abund table
  original_BC <- vegdist(inlandeco_otu_inland_t, method="bray")
  #calculate correlation between original and truncated dissimilarity matrices
  mantel_output<-mantel(original_BC, truncated_BC, permutations = 999, method='pearson')
  #output the pvalue
  #output the pvalue
  mantel_results <- as.data.frame(cbind(mantel_output$statistic,mantel_output$signif))
  
}

#calculate pairwise dissim matrix for truncated dataset (based on relative abundances)
Thresholds<- seq(0.000001, 0.009, 0.0001)

results<-t(as.data.frame(sapply(Thresholds,Rare_vs_abundant_calc)))
Thresholds<-as.data.frame(Thresholds)
Thresh_Pvalue<-as.data.frame(cbind(Thresholds,results))
colnames(Thresh_Pvalue)<-c("Threshold","Statistic","Pvalue")
Thresh_Pvalue$AdjP <-p.adjust(Thresh_Pvalue$Pvalue, method = "fdr")


Thresh_Fig<- ggplot(Thresh_Pvalue, aes(Threshold, AdjP))+
  geom_line()+
  geom_point()+
  ylab("Statistic")+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle("(A)")+
  xlab("Threshold")
Thresh_Fig





#determining rare vs abundant taxa for Bulk soil.

#import otu table (remove "#Constructed from biom file " from first line before importing).
mim_otu<-read.table("zotu_table_forR.txt", row.names=1, header=T)

#subset data to include only Coastal Ecotype at inland site
mim_meta<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_meta<-rownames(subset(mim_meta, Site_planted =="Inland" & Rhizo_Bulk =="Bulk"))
mim_meta 
bulk_otu_inland <- subset(mim_otu, select=colnames(mim_otu) %in% mim_meta)
#remove OTUs with zero abundance
bulk_otu_inland <-bulk_otu_inland[ rowSums(bulk_otu_inland)!=0, ]

#make new dataframe housing the relative abundances of OTUs in the original dataset.
OTU_relabund<- as.data.frame(rowSums(bulk_otu_inland)/sum(rowSums(bulk_otu_inland)))
row.names(OTU_relabund)<- rownames(bulk_otu_inland)
colnames(OTU_relabund)<-c('RelAbund')

#set up loop:

#function
Rare_vs_abundant_calc<-function(threshold){
  #first subset list of OTU names to exclude those below given threshold.
  Rare_OTUs<-rownames(subset(OTU_relabund, RelAbund >= threshold))
  #then remove those OTUs from the OTU table to get the truncated dataset.
  truncated_mim_otu<-subset(bulk_otu_inland, rownames(bulk_otu_inland) %in% Rare_OTUs)
  #transpose truncated OTU table
  truncated_mim_otu_t<-t(truncated_mim_otu)
  #calculate Bray-Curtis dissimilarity matrix for truncated OTU table
  truncated_BC <- vegdist(truncated_mim_otu_t, method="bray")
  #transpose original OTU table
  bulk_otu_inland_t<-t(bulk_otu_inland)
  #calculate Bray-Curtis dissimilarity matrix for rel abund table
  original_BC <- vegdist(bulk_otu_inland_t, method="bray")
  #calculate correlation between original and truncated dissimilarity matrices
  mantel_output<-mantel(original_BC, truncated_BC, permutations = 999, method='pearson')
  #output the pvalue
  mantel_results <- as.data.frame(cbind(mantel_output$statistic,mantel_output$signif))
  
}

#calculate pairwise dissim matrix for truncated dataset (based on relative abundances)
Thresholds<- seq(0.000001, 0.009, 0.0001)

results<-t(as.data.frame(sapply(Thresholds,Rare_vs_abundant_calc)))
Thresholds<-as.data.frame(Thresholds)
Thresh_Pvalue<-as.data.frame(cbind(Thresholds,results))
colnames(Thresh_Pvalue)<-c("Threshold","Statistic","Pvalue")
Thresh_Pvalue$AdjP <-p.adjust(Thresh_Pvalue$Pvalue, method = "fdr")

Thresh_Fig<- ggplot(Thresh_Pvalue, aes(Threshold, AdjP))+
  geom_line()+
  geom_point()+
  ylab("Statistic")+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle("(A)")+
  xlab("Threshold")
Thresh_Fig



#split each group into rare vs abundant taxa.

#plot phylogenetic diversity for rare vs abundant taxa.




