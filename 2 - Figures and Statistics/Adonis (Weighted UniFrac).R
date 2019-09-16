#Perform PERMANOVA tests to assess variables affecting community composition

library(vegan)

#load metadata
#remove "#" from first line before importing
#remove "BarcodeSequence" column before importing.
#remove "LinkerPrimerSequence" column before importing.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, sep='\t', row.names=1)

#load weighted unifrac distance matrix
mim_WU<-read.table("weighted_unifrac_single_rare.txt",header=T,sep='\t',row.names=1)

#Look at main effects of Rhizo_Bulk and Site_planted
set.seed(1)
adonis(mim_WU~Rhizo_Bulk*Site_planted, data=mim_envir)
set.seed(1)
adonis(mim_WU~Site_planted*Rhizo_Bulk, data=mim_envir)

#Interaction is significant. 


#Divide dataset into each ecotype, then test effect of site.

#Subset full metadata to only include inland ecotype samples
mim_envir_Inland<-subset(mim_envir, Origin =="Inland" & Rhizo_Bulk == "Rhizo")
#Subset full dissimilarity matrix to only include Inland ecotype samples
#first get the list of inland ecotype saamples
mim_inlandecoonly<-rownames(mim_envir_Inland)
mim_inlandecoonly
#select only those samples from dissimilarity matrix columns.
mim_WU_Inland <- subset(mim_WU, select=c(Mim109, Mim181, Mim182, Mim190, Mim191, Mim1,   Mim224, Mim228 ,Mim22,  Mim237, Mim238, Mim250 ,Mim255, Mim32  ,Mim66,  Mim67  ,Mim83,  Mim87,  Mim93,  Mim98))
#now do the same thing for the dissimilarity matrix rows.
list <- c("Mim109", "Mim181", "Mim182", "Mim190", "Mim191", "Mim1",   "Mim224", "Mim228", "Mim22"  ,"Mim237","Mim238", "Mim250", "Mim255", "Mim32",  "Mim66",  "Mim67",  "Mim83",  "Mim87",  "Mim93",  "Mim98")
mim_WU_Inland<-mim_WU_Inland[(row.names(mim_WU_Inland) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Inland~Site_planted, data=mim_envir_Inland)
#Homogeneity of dispersions test
dist <- as.dist(mim_WU_Inland)
set.seed(1)
result <- betadisper(dist, mim_envir_Inland$Site_planted, type=c("centroid"))
set.seed(1)
permutest(result)

#Subset full metadata to only include inland ecotype samples
mim_envir_Coastal<-subset(mim_envir, Origin =="Coastal" & Rhizo_Bulk == "Rhizo")
#Subset full dissimilarity matrix to only include Inland ecotype samples
#first get the list of inland ecotype saamples
mim_coastalecoonly<-rownames(mim_envir_Coastal)
mim_coastalecoonly
#select only those samples from dissimilarity matrix columns.
mim_WU_Coastal <- subset(mim_WU, select=c(Mim11,  Mim12,  Mim141, Mim143, Mim151, Mim152, Mim161, Mim162, Mim172, Mim176, Mim212, Mim21  ,Mim2 ,  Mim31,  Mim46,  Mim47,  Mim57,  Mim58,  Mim97))
#now do the same thing for the dissimilarity matrix rows.
list <- c("Mim11",  "Mim12",  "Mim141", "Mim143", "Mim151", "Mim152", "Mim161", "Mim162", "Mim172", "Mim176", "Mim212", "Mim21"  ,"Mim2" ,  "Mim31",  "Mim46",  "Mim47",  "Mim57",  "Mim58",  "Mim97")
mim_WU_Coastal<-mim_WU_Coastal[(row.names(mim_WU_Coastal) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Coastal~Site_planted, data=mim_envir_Coastal)
#Homogeneity of dispersions test
dist <- as.dist(mim_WU_Coastal)
set.seed(1)
result <- betadisper(dist, mim_envir_Coastal$Site_planted, type=c("centroid"))
set.seed(1)
permutest(result)





####Must divide dataset into "Inland" vs "Coastal" sites, then test for
#differences between rhizosphere and bulk soil samples.

#Subset full metadata to only include "Inland" site samples
mim_envir_Inland<-subset(mim_envir, Site_planted =="Inland")
#Subset full dissimilarity matrix to only include Inland site samples
#first get the list of samples collected Inland
mim_inlandonly<-rownames(mim_envir_Inland)
mim_inlandonly
#select only those samples from dissimilarity matrix columns.
mim_WU_Inland <- subset(mim_WU, select=c(Mim109, Mim116, Mim11,  Mim124, Mim129, Mim12,  Mim132, Mim139, Mim1,   Mim21, Mim22,  Mim2,   Mim31,  Mim32,  Mim46,  Mim47,  Mim57,  Mim58,  Mim66,  Mim67, Mim83,  Mim87,  Mim93,  Mim97,  Mim98))
#now do the same thing for the dissimilarity matrix rows.
list <- c("Mim109", "Mim116", "Mim11",  "Mim124", "Mim129", "Mim12",  "Mim132", "Mim139", "Mim1",   "Mim21", "Mim22","Mim2",   "Mim31",  "Mim32",  "Mim46",  "Mim47",  "Mim57",  "Mim58",  "Mim66",  "Mim67", "Mim83",  "Mim87",  "Mim93",  "Mim97",  "Mim98")
mim_WU_Inland<-mim_WU_Inland[(row.names(mim_WU_Inland) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Inland~Rhizo_Bulk, data=mim_envir_Inland)
#Homogeneity of dispersions test
dist <- as.dist(mim_WU_Inland)
set.seed(1)
result <- betadisper(dist, mim_envir_Inland$Rhizo_Bulk, type=c("centroid"))
set.seed(1)
permutest(result)



#Subset full metadata to only include "Coastal" site samples
mim_envir_Coastal<-subset(mim_envir, Site_planted =="Coastal")
#Subset full dissimilarity matrix to only include "Coastal site samples
#first get the list of samples collected Coastal
mim_Coastalonly<-rownames(mim_envir_Coastal)
mim_Coastalonly
#select only those samples from full dissimilarity matrix columns.
mim_WU_Coastal <- subset(mim_WU, select=c(Mim141, Mim143, Mim151, Mim152, Mim161, Mim162, Mim172, Mim176, Mim181, Mim182,Mim190, Mim191, Mim212, Mim224, Mim228, Mim237, Mim238, Mim250, Mim255, Mim261,Mim266, Mim271, Mim277, Mim286))
#now select only those samples from the the full dissimilarity matrix rows.
list <- c("Mim141", "Mim143", "Mim151", "Mim152", "Mim161", "Mim162", "Mim172", "Mim176", "Mim181", "Mim182","Mim190", "Mim191", "Mim212", "Mim224", "Mim228", "Mim237", "Mim238", "Mim250", "Mim255", "Mim261","Mim266", "Mim271", "Mim277", "Mim286")
mim_WU_Coastal<-mim_WU_Coastal[(row.names(mim_WU_Coastal) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Coastal~Rhizo_Bulk, data=mim_envir_Coastal)
#Homogeneity of dispersions test
dist <- as.dist(mim_WU_Coastal)
set.seed(1)
result <- betadisper(dist, mim_envir_Coastal$Rhizo_Bulk, type=c("centroid"))
set.seed(1)
permutest(result)




####
#Next, investigate effects of each ecotype (versus bulk soil) at each site.

#Subset "Inland" site samples to get only bulk soils and Inland ecotypes
mim_envir_Inland_BulkInlandeco<-subset(mim_envir_Inland, Genotype=="BulkSoil" | Origin=="Inland")
#Subset full dissimilarity matrix to only include these samples
#first get the list of samples collected Inland
mim_bulkandinlandecoonly<-rownames(mim_envir_Inland_BulkInlandeco)
mim_bulkandinlandecoonly
#select only those samples from dissimilarity matrix columns.
mim_WU_Inland_BulkInlandeco <- subset(mim_WU, select=c(Mim109,Mim116,Mim124,Mim129,Mim132,Mim139,Mim1,Mim22,Mim32,Mim66,Mim67,Mim83,Mim87,Mim93,Mim98))
#now do the same things for the dissimilarity matrix rows.
list <- c("Mim109", "Mim116", "Mim124", "Mim129", "Mim132", "Mim139", "Mim1",   "Mim22",  "Mim32",  "Mim66", "Mim67",  "Mim83",  "Mim87",  "Mim93",  "Mim98")
mim_WU_Inland_BulkInlandeco<-mim_WU_Inland_BulkInlandeco[(row.names(mim_WU_Inland_BulkInlandeco) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Inland_BulkInlandeco~Rhizo_Bulk, data=mim_envir_Inland_BulkInlandeco)
#Homogeneity of dispersions test
dist <- as.dist(mim_WU_Inland_BulkInlandeco)
set.seed(1)
result <- betadisper(dist, mim_envir_Inland_BulkInlandeco$Rhizo_Bulk, type=c("centroid"))
set.seed(1)
permutest(result)


#Subset "Inland" site samples to get only bulk soils and Coastal ecotypes
mim_envir_Inland_BulkCoastaleco<-subset(mim_envir_Inland, Genotype=="BulkSoil" | Origin=="Coastal")
#Subset full dissimilarity matrix to only include these samples
#first get the list of samples collected Inland
mim_bulkandcoastalecoonly<-rownames(mim_envir_Inland_BulkCoastaleco)
mim_bulkandcoastalecoonly
#select only those samples from dissimilarity matrix columns.
mim_WU_Inland_BulkCoastaleco <- subset(mim_WU, select=c(Mim116,Mim11,Mim124,Mim129,Mim12,Mim132,Mim139,Mim21,Mim2,Mim31,Mim46,Mim47,Mim57,Mim58,Mim97))
#now do the same things for the dissimilarity matrix rows.
list <- c("Mim116", "Mim11",  "Mim124", "Mim129", "Mim12",  "Mim132" ,"Mim139", "Mim21" , "Mim2",   "Mim31", "Mim46",  "Mim47",  "Mim57",  "Mim58",  "Mim97" )
mim_WU_Inland_BulkCoastaleco<-mim_WU_Inland_BulkCoastaleco[(row.names(mim_WU_Inland_BulkCoastaleco) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Inland_BulkCoastaleco~Rhizo_Bulk, data=mim_envir_Inland_BulkCoastaleco)
#Homogeneity of dispersions test
dist <- as.dist(mim_WU_Inland_BulkCoastaleco)
set.seed(1)
result <- betadisper(dist, mim_envir_Inland_BulkCoastaleco$Rhizo_Bulk, type=c("centroid"))
set.seed(1)
permutest(result)


#Subset "Coastal" site samples to get only bulk soils and Inland ecotypes
mim_envir_Coastal_BulkInlandeco<-subset(mim_envir_Coastal, Genotype=="BulkSoil" | Origin=="Inland")
#Subset full dissimilarity matrix to only include these samples
#first get the list of samples
mim_bulkandInlandecoonly<-rownames(mim_envir_Coastal_BulkInlandeco)
mim_bulkandInlandecoonly
#select only those samples from dissimilarity matrix columns.
mim_WU_Coastal_BulkInlandeco <- subset(mim_WU, select=c(Mim181, Mim182 ,Mim190, Mim191, Mim224 ,Mim228 ,Mim237 ,Mim238, Mim250, Mim255,Mim261, Mim266 ,Mim271 ,Mim277, Mim286))
#now do the same things for the dissimilarity matrix rows.
list <- c("Mim181", "Mim182" ,"Mim190", "Mim191", "Mim224" ,"Mim228" ,"Mim237" ,"Mim238", "Mim250", "Mim255","Mim261", "Mim266" ,"Mim271" ,"Mim277", "Mim286")
mim_WU_Coastal_BulkInlandeco<-mim_WU_Coastal_BulkInlandeco[(row.names(mim_WU_Coastal_BulkInlandeco) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Coastal_BulkInlandeco~Rhizo_Bulk, data=mim_envir_Coastal_BulkInlandeco)
#Homogeneity of dispersions test
dist <- as.dist(mim_WU_Coastal_BulkInlandeco)
set.seed(1)
result <- betadisper(dist, mim_envir_Coastal_BulkInlandeco$Rhizo_Bulk, type=c("centroid"))
set.seed(1)
permutest(result)


#Subset "Coastal" site samples to get only bulk soils and Coastal ecotypes
mim_envir_Coastal_BulkCoastaleco<-subset(mim_envir_Coastal, Genotype=="BulkSoil" | Origin=="Coastal")
#Subset full dissimilarity matrix to only include these samples
#first get the list of samples collected Coastal
mim_bulkandCoastalecoonly<-rownames(mim_envir_Coastal_BulkCoastaleco)
mim_bulkandCoastalecoonly
#select only those samples from dissimilarity matrix columns.
mim_WU_Coastal_BulkCoastaleco <- subset(mim_WU, select=c(Mim141,Mim143,Mim151,Mim152,Mim161,Mim162,Mim172,Mim176,Mim212,Mim261,Mim266,Mim271,Mim277,Mim286))
#now do the same things for the dissimilarity matrix rows.
list <- c("Mim141", "Mim143", "Mim151", "Mim152", "Mim161", "Mim162", "Mim172", "Mim176", "Mim212", "Mim261", "Mim266", "Mim271", "Mim277", "Mim286")
mim_WU_Coastal_BulkCoastaleco<-mim_WU_Coastal_BulkCoastaleco[(row.names(mim_WU_Coastal_BulkCoastaleco) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Coastal_BulkCoastaleco~Rhizo_Bulk, data=mim_envir_Coastal_BulkCoastaleco)
#Homogeneity of dispersions test
dist <- as.dist(mim_WU_Coastal_BulkCoastaleco)
set.seed(1)
result <- betadisper(dist, mim_envir_Coastal_BulkCoastaleco$Rhizo_Bulk, type=c("centroid"))
set.seed(1)
permutest(result)




####
#Finally, investigate effects of Ecotype and Genotype at each site.

#Subset Inland metadata to only include "rhizosphere" samples
mim_envir_Inland_rhizo<-subset(mim_envir_Inland, Rhizo_Bulk =="Rhizo")
mim_envir_Inland_rhizo
#Subset Inland dissimilarity matrix to only include "Inland site samples
#first get the list of samples collected Inland
mim_Inlandrhizoonly<-rownames(mim_envir_Inland_rhizo)
mim_Inlandrhizoonly
#select only those samples from Inland dissimilarity matrix columns.
mim_WU_Inland_rhizo <- subset(mim_WU_Inland, select=c(Mim109, Mim11,  Mim12,  Mim1,   Mim21,  Mim22,  Mim2,   Mim31,  Mim32,  Mim46, Mim47,  Mim57,  Mim58,  Mim66,  Mim67,  Mim83,  Mim87,  Mim93,  Mim97,  Mim98))
#now do the same things for the dissimilarity matrix rows.
list <- c("Mim109", "Mim11",  "Mim12",  "Mim1",   "Mim21",  "Mim22",  "Mim2",   "Mim31",  "Mim32",  "Mim46", "Mim47",  "Mim57",  "Mim58",  "Mim66",  "Mim67",  "Mim83",  "Mim87",  "Mim93",  "Mim97",  "Mim98")
mim_WU_Inland_rhizo<-mim_WU_Inland_rhizo[(row.names(mim_WU_Inland_rhizo) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Inland_rhizo~Origin,data=mim_envir_Inland_rhizo)
#Homogeneity of dispersions test
dist <- as.dist(mim_WU_Inland_rhizo)
set.seed(1)
result <- betadisper(dist, mim_envir_Inland_rhizo$Origin, type=c("centroid"))
set.seed(1)
permutest(result)
#With strata=Origin, you limit permutations (re-labeling) to be between only those samples of the same Ecotype.
set.seed(1)
adonis(mim_WU_Inland_rhizo~Genotype,data=mim_envir_Inland_rhizo, strata=mim_envir_Inland_rhizo$Origin)

#Subset Inland rhizosphere metadata to only include "Coastal" ecotypes
mim_envir_Inland_rhizo_coastaleco<-subset(mim_envir_Inland_rhizo, Origin =="Coastal")
mim_envir_Inland_rhizo_coastaleco
#Subset Inland rhizosphere dissimilarity matrix to only include Coastal Ecotypes
#first get the list of rhizo samples collected Inland
mim_Inlandrhizoonly_coastaleco<-rownames(mim_envir_Inland_rhizo_coastaleco)
mim_Inlandrhizoonly_coastaleco
#select only those samples from Inland rhizo dissimilarity matrix columns.
mim_WU_Inland_rhizo_coastaleco <- subset(mim_WU_Inland_rhizo, select=c(Mim11, Mim12, Mim21, Mim2, Mim31, Mim46, Mim47, Mim57, Mim58, Mim97))
#now do the same things for the dissimilarity matrix rows.
list <- c("Mim11", "Mim12", "Mim21", "Mim2",  "Mim31", "Mim46", "Mim47", "Mim57", "Mim58", "Mim97")
mim_WU_Inland_rhizo_coastaleco<-mim_WU_Inland_rhizo_coastaleco[(row.names(mim_WU_Inland_rhizo_coastaleco) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Inland_rhizo_coastaleco~Genotype,data=mim_envir_Inland_rhizo_coastaleco)

#Subset Inland rhizosphere metadata to only include "Inland" ecotypes
mim_envir_Inland_rhizo_Inlandeco<-subset(mim_envir_Inland_rhizo, Origin =="Inland")
mim_envir_Inland_rhizo_Inlandeco
#Subset Inland rhizosphere dissimilarity matrix to only include Inland Ecotypes
#first get the list of rhizo samples collected Inland
mim_Inlandrhizoonly_Inlandeco<-rownames(mim_envir_Inland_rhizo_Inlandeco)
mim_Inlandrhizoonly_Inlandeco
#select only those samples from Inland rhizo dissimilarity matrix columns.
mim_WU_Inland_rhizo_Inlandeco <- subset(mim_WU_Inland_rhizo, select=c(Mim109, Mim1,Mim22,Mim32,  Mim66,  Mim67,  Mim83,  Mim87,  Mim93,  Mim98 ))
#now do the same things for the dissimilarity matrix rows.
list <- c("Mim109", "Mim1",   "Mim22",  "Mim32",  "Mim66",  "Mim67",  "Mim83",  "Mim87",  "Mim93",  "Mim98")
mim_WU_Inland_rhizo_Inlandeco<-mim_WU_Inland_rhizo_Inlandeco[(row.names(mim_WU_Inland_rhizo_Inlandeco) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Inland_rhizo_Inlandeco~Genotype,data=mim_envir_Inland_rhizo_Inlandeco)










#Subset Coastal metadata to only include "rhizosphere" samples
mim_envir_Coastal_rhizo<-subset(mim_envir_Coastal, Rhizo_Bulk =="Rhizo")
mim_envir_Coastal_rhizo
#Subset Coastal dissimilarity matrix to only include Coastal site samples
#first get the list of samples collected Coastal
mim_Coastalrhizoonly<-rownames(mim_envir_Coastal_rhizo)
mim_Coastalrhizoonly
#select only those samples from Coastal dissimilarity matrix columns.
mim_WU_Coastal_rhizo <- subset(mim_WU_Coastal, select=c(Mim141, Mim143, Mim151, Mim152, Mim161, Mim162, Mim172, Mim176, Mim181, Mim182, Mim190, Mim191, Mim212, Mim224, Mim228, Mim237, Mim238, Mim250, Mim255))
#now do the same things for the dissimilarity matrix rows.
list <- c("Mim141", "Mim143", "Mim151", "Mim152", "Mim161", "Mim162","Mim172", "Mim176", "Mim181", "Mim182","Mim190", "Mim191", "Mim212", "Mim224", "Mim228","Mim237", "Mim238", "Mim250", "Mim255")
mim_WU_Coastal_rhizo<-mim_WU_Coastal_rhizo[(row.names(mim_WU_Coastal_rhizo) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Coastal_rhizo~Origin,data=mim_envir_Coastal_rhizo)
#Homogeneity of dispersions test
dist <- as.dist(mim_WU_Coastal_rhizo)
set.seed(1)
result <- betadisper(dist, mim_envir_Coastal_rhizo$Origin, type=c("centroid"))
set.seed(1)
permutest(result)
#With strata=Origin, you limit permutations (re-labeling) to be between only those samples of the same Ecotype.
set.seed(1)
adonis(mim_WU_Coastal_rhizo~Genotype,data=mim_envir_Coastal_rhizo, strata=mim_envir_Coastal_rhizo$Origin)

#Subset Coastal rhizosphere metadata to only include "Coastal" ecotypes
mim_envir_Coastal_rhizo_coastaleco<-subset(mim_envir_Coastal_rhizo, Origin =="Coastal")
mim_envir_Coastal_rhizo_coastaleco
#Subset Coastal rhizosphere dissimilarity matrix to only include Coastal Ecotypes
#first get the list of rhizo samples collected Coastal
mim_Coastalrhizoonly_coastaleco<-rownames(mim_envir_Coastal_rhizo_coastaleco)
mim_Coastalrhizoonly_coastaleco
#select only those samples from Coastal rhizo dissimilarity matrix columns.
mim_WU_Coastal_rhizo_coastaleco <- subset(mim_WU_Coastal_rhizo, select=c(Mim141,Mim143,Mim151,Mim152,Mim161,Mim162,Mim172,Mim176,Mim212))
#now do the same things for the dissimilarity matrix rows.
list <- c("Mim141", "Mim143", "Mim151", "Mim152", "Mim161", "Mim162", "Mim172", "Mim176", "Mim212")
mim_WU_Coastal_rhizo_coastaleco<-mim_WU_Coastal_rhizo_coastaleco[(row.names(mim_WU_Coastal_rhizo_coastaleco) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Coastal_rhizo_coastaleco~Genotype,data=mim_envir_Coastal_rhizo_coastaleco)

#Subset Coastal rhizosphere metadata to only include "Inland" ecotypes
mim_envir_Coastal_rhizo_Inlandeco<-subset(mim_envir_Coastal_rhizo, Origin =="Inland")
mim_envir_Coastal_rhizo_Inlandeco
#Subset Coastal rhizosphere dissimilarity matrix to only include Inland Ecotypes
#first get the list of rhizo samples collected Coastal
mim_Coastalrhizoonly_Inlandeco<-rownames(mim_envir_Coastal_rhizo_Inlandeco)
mim_Coastalrhizoonly_Inlandeco
#select only those samples from Coastal rhizo dissimilarity matrix columns.
mim_WU_Coastal_rhizo_Inlandeco <- subset(mim_WU_Coastal_rhizo, select=c("Mim181", "Mim182", "Mim190", "Mim191","Mim224", "Mim228", "Mim237", "Mim238", "Mim250", "Mim255"))
#now do the same things for the dissimilarity matrix rows.
list <- c("Mim181", "Mim182", "Mim190", "Mim191", "Mim224", "Mim228", "Mim237", "Mim238", "Mim250", "Mim255")
mim_WU_Coastal_rhizo_Inlandeco<-mim_WU_Coastal_rhizo_Inlandeco[(row.names(mim_WU_Coastal_rhizo_Inlandeco) %in% list), ]
#Finally, do Adonis test
set.seed(1)
adonis(mim_WU_Coastal_rhizo_Inlandeco~Genotype,data=mim_envir_Coastal_rhizo_Inlandeco)
