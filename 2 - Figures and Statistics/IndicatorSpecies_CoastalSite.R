#import zotu table (remove "#Constructed from biom file " from first line before importing).
mim_otu<-read.table("zotu_table_forR.txt", row.names=1, header=T)

#load metadata
#remove "#" from first line before importing
#remove "BarcodeSequence" column before importing.
#remove "LinkerPrimerSequence" column before importing.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, sep='\t', row.names=1)


#subset the full zOTU table for only coastal site rhizosphere samples.
mim_otu_coastalrhizo <- subset(mim_otu, select=c(rownames(subset(mim_envir, Rhizo_Bulk =="Rhizo" & Site_planted=="Coastal"))))
#transpose matrix for compatibility with 'indicspecies'
mim_otu_coastalrhizo_t<-as.data.frame(t(mim_otu_coastalrhizo))

#subset the full metadata table for only coastal site rhizosphere samples.
mim_envir_coastalrhizo <- mim_envir[c(rownames(subset(mim_envir, Rhizo_Bulk =="Rhizo" & Site_planted=="Coastal"))),]
#subset for only the "Origin" data.
mim_envir_coastalrhizo<-subset(mim_envir_coastalrhizo, select=c("Origin"))
#label each sample: 1 for Inland ecotypes, 2 for coastal ecotypes.
mim_envir_coastalrhizo$vec <- with(mim_envir_coastalrhizo, ifelse(Origin=="Inland", 1, 2))
mim_envir_coastalrhizo$vec



#Indicspecies analysis
library(indicspecies)

set.seed(1)
result<- multipatt(mim_otu_coastalrhizo_t, mim_envir_coastalrhizo$vec, duleg=FALSE,control=how(nperm=999))
summary(result)


#BUT NEED TO CORRECT FOR MULTIPLE COMPARISONS

#first, extract OTU names and p-values from the indicspecies result file.
OTUnames<-as.data.frame(row.names(result$sign))
OTUpvalues<-as.data.frame(result$sign$p.value)
#combine those two vectors
results<- cbind(OTUnames,OTUpvalues)
colnames(results)<-c("OTU","pvalue")

#remove rows with 'na', although it doesn't make a difference in adj pvalues.
results <- results[complete.cases(results),]

#calculate adjusted p-values.
results$adjpvalue<- p.adjust(results$pvalue, method = "fdr")

#p-value is 'NA' fr species whose maximum combination is the set of all combinations




#Apparently no indicator species. 
#Now compare zOTU relative abundance for inland vs coastal ecotypes at coastal site.

#First subset the coastal metadata for inland ecotypes only
inlandeco_coastal <- subset(mim_otu, select =c(rownames(subset(mim_envir_coastalrhizo, Origin=="Inland"))))

#next subset the coastal metadata for coastal ecotypes only
coastaleco_coastal <- subset(mim_otu, select =c(rownames(subset(mim_envir_coastalrhizo, Origin=="Coastal"))))

#T-tests comparing Coastal ecotypes versus inland ecotypes at coastal Site
ttest.out <- data.frame(matrix(nrow = nrow(mim_otu_coastalrhizo), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(mim_otu_coastalrhizo)
#Perform the t-tests.
for (i in seq(nrow(inlandeco_coastal))){
  test<-c()
  test[[i]]<-t.test(inlandeco_coastal[i,], coastaleco_coastal[i,], paired=FALSE, var.equal = TRUE)
  ttest.out[i,1] <- row.names(mim_otu_coastalrhizo)[i]
  ttest.out[i,2] <-(test[[i]]$statistic)
  ttest.out[i,3] <-(test[[i]]$parameter)
  ttest.out[i,4] <-(test[[i]]$p.value)
}
#Calculate adjusted p-values and add to table.
ttest.out[,5] <-p.adjust(ttest.out$p.value, method = "fdr")
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value", "adj.p.value")

####
#Apparently no OTUs significantly differ between ecotypes at the coastal site.
#Now exploring differences in presence/absence (FIDELITY)
#FIDELITY = proportion of samples in a given treatment in which that taxa occurs
#probability of finding the species in sites belonging to the site group.


#Each ecotype has 10 samples at the coastal site.

#First, make a new dataframe. 
PresAbs.df <- data.frame(matrix(nrow = nrow(inlandeco_coastal), ncol = 2))
rownames(PresAbs.df)<-row.names(inlandeco_coastal)
#name columns
colnames(PresAbs.df)<- c("Inland","Coastal")

#calculate number of samples with at least one read
PresAbs.df$Inland <- rowSums(inlandeco_coastal > 0)
PresAbs.df$Coastal <- rowSums(coastaleco_coastal > 0)


#choose OTUs present in at least 5 Inland samples, but zero Coastal samples
InlandInd <- (subset(PresAbs.df, Inland >="5" & Coastal=="0"))

#choose OTUs present in at least 5 Coastal samples, but zero Inland samples
CoastalInd <- (subset(PresAbs.df, Coastal >="5" & Inland=="0"))


rownames(CoastalInd)  #39 obs.
coastaleco_coastal["ZOTU2522",]
coastaleco_coastal["ZOTU1851",]
coastaleco_coastal["ZOTU1146",]


rownames(InlandInd) #39 obs
inlandeco_coastal["ZOTU2343",] 
inlandeco_coastal["ZOTU10425",]  
inlandeco_coastal["ZOTU5735",]  


#Now exploring differences in specificity
#Specificity = avg # reads across samples in a given treatment 
#divided by sum of averages across treatments)

#proportion of total occurrences occuring in a specific group.

#This conditional probability is called the specificity or positive
#predictive value of the species as indicator of the site group. 

#Find average # reads for coastal ecotype in inland environment
coastaleco_avg <- as.data.frame(rowMeans(coastaleco_coastal))
rownames(coastaleco_avg) <- row.names(coastaleco_coastal)
#Find average # reads for inland ecotype in inland environment
inlandeco_avg <- as.data.frame(rowMeans(inlandeco_coastal))
rownames(inlandeco_avg) <- row.names(inlandeco_coastal)
#Calculate specificity for each OTU
CoastalEco_specificity<-coastaleco_avg/(coastaleco_avg+inlandeco_avg)
InlandEco_specificity<-inlandeco_avg/(coastaleco_avg+inlandeco_avg)

