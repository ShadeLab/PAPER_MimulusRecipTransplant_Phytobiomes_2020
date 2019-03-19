#import zotu table (remove "#Constructed from biom file " from first line before importing).
mim_otu<-read.table("zotu_table_forR.txt", row.names=1, header=T)

#load metadata
#remove "#" from first line before importing
#remove "BarcodeSequence" column before importing.
#remove "LinkerPrimerSequence" column before importing.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, sep='\t', row.names=1)

#what are the inland site rhizosphere samples?
rownames(subset(mim_envir, Rhizo_Bulk =="Rhizo" & Site_planted=="Inland"))

#subset the full zOTU table for only inland site rhizosphere samples.
mim_otu_inlandrhizo <- subset(mim_otu, select=c("Mim109", "Mim11",  "Mim12",  "Mim1",   "Mim21",  "Mim22",  "Mim2",   "Mim31",  "Mim32",  "Mim46", "Mim47",  "Mim57",  "Mim58",  "Mim66",  "Mim67",  "Mim83",  "Mim87",  "Mim93",  "Mim97",  "Mim98"))
#transpose matrix for compatibility with 'indicspecies'
mim_otu_inlandrhizo_t<-as.data.frame(t(mim_otu_inlandrhizo))

#subset the full metadata table for only inland site rhizosphere samples.
mim_envir_inlandrhizo <- mim_envir[c("Mim109", "Mim11",  "Mim12",  "Mim1",  "Mim21",  "Mim22",  "Mim2",   "Mim31",  "Mim32",  "Mim46", "Mim47",  "Mim57",  "Mim58",  "Mim66",  "Mim67",  "Mim83",  "Mim87",  "Mim93",  "Mim97",  "Mim98"),]
#subset for only the "Origin" data.
mim_envir_inlandrhizo<-subset(mim_envir_inlandrhizo, select=c("Origin"))
#label each sample: 1 for Inland ecotypes, 2 for coastal ecotypes.
mim_envir_inlandrhizo$vec <- with(mim_envir_inlandrhizo, ifelse(Origin=="Inland", 1, 2))
mim_envir_inlandrhizo$vec



#Indicspecies analysis
library(indicspecies)

set.seed(1)
result<- multipatt(mim_otu_inlandrhizo_t, mim_envir_inlandrhizo$vec, duleg=FALSE,control=how(nperm=999))
#summary(result)


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
#Now compare zOTU relative abundance for inland vs coastal ecotypes at Inland site.

#First subset the inland metadata for inland ecotypes only
inlandeco_inlandnames<-rownames(subset(mim_envir_inlandrhizo, Origin=="Inland"))
inlandeco_inlandnames
inlandeco_inland <- subset(mim_otu, select =c("Mim109", "Mim1","Mim22","Mim32","Mim66","Mim67","Mim83","Mim87","Mim93","Mim98" ))

#next subset the inland metadata for coastal ecotypes only
coastaleco_inlandnames<-rownames(subset(mim_envir_inlandrhizo, Origin=="Coastal"))
coastaleco_inlandnames
coastaleco_inland <- subset(mim_otu, select =c("Mim11","Mim12","Mim21","Mim2","Mim31","Mim46","Mim47","Mim57","Mim58","Mim97"))



#T-tests comparing Coastal ecotypes versus inland ecotypes at Inland Site
ttest.out <- data.frame(matrix(nrow = nrow(mim_otu_inlandrhizo), ncol = 4))
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value")
rownames(ttest.out) <- row.names(mim_otu_inlandrhizo)
#Perform the t-tests.
for (i in seq(nrow(inlandeco_inland))){
  test<-c()
  test[[i]]<-t.test(inlandeco_inland[i,], coastaleco_inland[i,], paired=FALSE, var.equal = TRUE)
  ttest.out[i,1] <- row.names(mim_otu_inlandrhizo)[i]
  ttest.out[i,2] <-(test[[i]]$statistic)
  ttest.out[i,3] <-(test[[i]]$parameter)
  ttest.out[i,4] <-(test[[i]]$p.value)
}
#Calculate adjusted p-values and add to table.
ttest.out[,5] <-p.adjust(ttest.out$p.value, method = "fdr")
colnames(ttest.out) <- c("Taxa", "t-statistic","df", "p.value", "adj.p.value")


####
#Apparently no OTUs significantly differ between ecotypes at the Inland site.
#Now exploring differences in presence/absence (FIDELITY)
#FIDELITY = proportion of samples in a given treatment in which that taxa occurs
#probability of finding the species in sites belonging to the site group.


#Each ecotype has 10 samples at the Inland site.

#First, make a new dataframe. 
PresAbs.df <- data.frame(matrix(nrow = nrow(inlandeco_inland), ncol = 2))
rownames(PresAbs.df)<-row.names(inlandeco_inland)
#name columns
colnames(PresAbs.df)<- c("Inland","Coastal")

#put a 1 if present at least 8 samples.
PresAbs.df$Inland <- rowSums(inlandeco_inland > 0)
PresAbs.df$Coastal <- rowSums(coastaleco_inland > 0)


#choose OTUs present in at least 6 Inland samples, but zero Coastal samples
InlandInd <- (subset(PresAbs.df, Inland >="6" & Coastal=="0"))

#choose OTUs present in at least 6 Coastal samples, but zero Inland samples
CoastalInd <- (subset(PresAbs.df, Coastal >="6" & Inland=="0"))


rownames(CoastalInd)
coastaleco_inland["ZOTU14147",]
coastaleco_inland["ZOTU16159",]
coastaleco_inland["ZOTU15124",]

rownames(InlandInd)
inlandeco_inland["ZOTU12953",] 
inlandeco_inland["ZOTU6975",]  
inlandeco_inland["ZOTU1624",]  
inlandeco_inland["ZOTU10581",] 
inlandeco_inland["ZOTU3045",]  
inlandeco_inland["ZOTU2227",]  
inlandeco_inland["ZOTU8375",] 
inlandeco_inland["ZOTU8727",]


#Now exploring differences in specificity
#Specificity = avg # reads across samples in a given treatment 
#divided by sum of averages across treatments)

#proportion of total occurrences occuring in a specific group.

#This conditional probability is called the specificity or positive
#predictive value of the species as indicator of the site group. 

#Find average # reads for coastal ecotype in inland environment
coastaleco_avg <- as.data.frame(rowMeans(coastaleco_inland))
rownames(coastaleco_avg) <- row.names(coastaleco_inland)
#Find average # reads for inland ecotype in inland environment
inlandeco_avg <- as.data.frame(rowMeans(inlandeco_inland))
rownames(inlandeco_avg) <- row.names(inlandeco_inland)
#Calculate specificity for each OTU
CoastalEco_specificity<-coastaleco_avg/(coastaleco_avg+inlandeco_avg)
InlandEco_specificity<-inlandeco_avg/(coastaleco_avg+inlandeco_avg)

