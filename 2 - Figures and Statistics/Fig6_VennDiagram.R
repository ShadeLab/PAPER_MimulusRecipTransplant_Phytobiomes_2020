#import otu table (remove "#Constructed from biom file " from first line before importing).
mim_otu<-read.table("zotu_table_forR.txt", row.names=1, header=T)
#make sure we get rid of Coastal site samples.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted =="Inland"))
mim_envir 
mim_otu<- subset(mim_otu, select= colnames(mim_otu) %in% mim_envir)

#First, subset the full dataframe for just inland site coastal ecotype samples.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted =="Inland" & Origin=="Coastal" & Rhizo_Bulk =="Rhizo"))
mim_envir 
inlandrhizo_coastaleco <- subset(mim_otu, select=colnames(mim_otu) %in% mim_envir)
#determine present OTUs
inlandrhizo_coastaleco$Sums <-rowSums(inlandrhizo_coastaleco)
coastaleco<-rownames(subset(inlandrhizo_coastaleco, Sums >= 1))
#this is 8931 OTUs.

#next, subset the full dataframe for just inland site inland ecotype samples.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Site_planted =="Inland" & Origin=="Inland" & Rhizo_Bulk =="Rhizo"))
mim_envir 
inlandrhizo_inlandeco <- subset(mim_otu, select=colnames(mim_otu) %in% mim_envir)
#determine present OTUs
inlandrhizo_inlandeco$Sums <-rowSums(inlandrhizo_inlandeco)
inlandeco<-rownames(subset(inlandrhizo_inlandeco, Sums >= 1))
#this is 9839 OTUs.


#lastly subset the full dataframe for just inland site bulk soils.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')
mim_envir<-rownames(subset(mim_envir, Origin =="Inland" & Rhizo_Bulk =="Bulk"))
mim_envir 
inlandbulk <- subset(mim_otu, select=colnames(mim_otu) %in% mim_envir)
#determine present OTUs
inlandbulk$Sums <-rowSums(inlandbulk)
bulk<-rownames(subset(inlandbulk, Sums >= 1))
#this is 8094 OTUs.



#First generate basic Venn Diagram.
library(gplots)
ItemsList <- venn(data = list("CoastalEco" = coastaleco,"InlandEco" = inlandeco,"Bulk" = bulk), show.plot=TRUE)
#This shows us what our final result should look like.
#It also helps us calculate proper values for each pairwise comparison
#This is important input for the 'VennDiagram' package, which is much more customizable.


####Preparing to use 'venndiagram' package.
#Verify total number of OTUs in or intersecting with CoastalEco. Should be 8931.
length(attr(ItemsList,"intersections")$CoastalEco)+length(attr(ItemsList,"intersections")$'CoastalEco:Bulk')+length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco')+length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco:Bulk')
#Verify total number of OTUs in or intersecting with InlandEco. Should be 9839
length(attr(ItemsList,"intersections")$InlandEco)+length(attr(ItemsList,"intersections")$'InlandEco:Bulk')+length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco')+length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco:Bulk')
#Verify total number of OTUs in or intersecting with Bulk. Should be 8094.
length(attr(ItemsList,"intersections")$Bulk)+length(attr(ItemsList,"intersections")$'InlandEco:Bulk')+length(attr(ItemsList,"intersections")$'CoastalEco:Bulk')+length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco:Bulk')
#Calculate number of OTUs shared by Coastal and Inland. Should be 7774.
length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco')+length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco:Bulk')
#Calculate number of OTUs shared by Inland and Bulk. Should be 7121.
length(attr(ItemsList,"intersections")$'InlandEco:Bulk')+length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco:Bulk')
#Calculate number of OTUs shared by Coastal and Bulk. Should be 6706.
length(attr(ItemsList,"intersections")$'CoastalEco:Bulk')+length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco:Bulk')
#Calculate number of OTUs shared by Coastal, Inland, and Bulk. Should be 6290.
length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco:Bulk')


#Prepare VennDiagram using the above values.
#the values should exactly match the diagram made with the 'venn' function above.
library(VennDiagram)
grid.newpage()
v <-draw.triple.venn(area1 = 8931, area2 = 9839, area3 = 8094, n12 = 7774, n23 = 7121, n13 = 6706,  n123 = 6290, category = c("Coastal", "Inland", "Bulk"), lty = "blank",     fill = c("blue", "red", "green"))



#Calculate number of reads in Coastal ecotype only, and percent of total reads
CoastalOnly_OTUs <- subset(mim_otu, rownames(mim_otu) %in% attr(ItemsList,"intersections")$CoastalEco)
sum(rowSums(CoastalOnly_OTUs)) 
sum(rowSums(CoastalOnly_OTUs))/sum(rowSums(mim_otu))*100 
#remove columns that sum to zero. Then calc avg rel abund for each OTU.
CoastalOnly_OTUs<- CoastalOnly_OTUs[,colSums(CoastalOnly_OTUs) > 0]
CoastalOnly_OTUs$mean<- rowMeans(CoastalOnly_OTUs)
#find relative abundance of each of those OTUs
CoastalOnly_OTUs$relabund <- CoastalOnly_OTUs$mean/22354
#find average percent relative abundance for each of those OTUs
mean(CoastalOnly_OTUs$relabund)*100


#Calculate number of reads in Inland ecotype only
InlandOnly_OTUs <- subset(mim_otu, rownames(mim_otu) %in% attr(ItemsList,"intersections")$InlandEco)
sum(rowSums(InlandOnly_OTUs)) 
sum(rowSums(InlandOnly_OTUs))/sum(rowSums(mim_otu))*100
#remove columns that sum to zero. Then calc avg rel abund for each OTU.
InlandOnly_OTUs<- InlandOnly_OTUs[,colSums(InlandOnly_OTUs) > 0]
InlandOnly_OTUs$mean<- rowMeans(InlandOnly_OTUs)
#find relative abundance of each of those OTUs
InlandOnly_OTUs$relabund <- InlandOnly_OTUs$mean/22354
#find average percent relative abundance for each of those OTUs
mean(InlandOnly_OTUs$relabund)*100


#Calculate number of reads in bulk soil only
BulkOnly_OTUs <- subset(mim_otu, rownames(mim_otu) %in% attr(ItemsList,"intersections")$Bulk)
sum(rowSums(BulkOnly_OTUs)) 
sum(rowSums(BulkOnly_OTUs))/sum(rowSums(mim_otu))*100
#remove columns that sum to zero. Then calc avg rel abund for each OTU.
BulkOnly_OTUs<- BulkOnly_OTUs[,colSums(BulkOnly_OTUs) > 0]
BulkOnly_OTUs$mean<- rowMeans(BulkOnly_OTUs)
#find relative abundance of each of those OTUs
BulkOnly_OTUs$relabund <- BulkOnly_OTUs$mean/22354
#find average percent relative abundance for each of those OTUs
mean(BulkOnly_OTUs$relabund)*100

#Calculate number of reads in Coastal:Inland intersection
CoastalInland_OTUs <- subset(mim_otu, rownames(mim_otu) %in% attr(ItemsList,"intersections")$'CoastalEco:InlandEco')
sum(rowSums(CoastalInland_OTUs)) 
sum(rowSums(CoastalInland_OTUs))/sum(rowSums(mim_otu))*100
#remove columns that sum to zero. Then calc avg rel abund for each OTU.
CoastalInland_OTUs<- CoastalInland_OTUs[,colSums(CoastalInland_OTUs) > 0]
CoastalInland_OTUs$mean<- rowMeans(CoastalInland_OTUs)
#find relative abundance of each of those OTUs
CoastalInland_OTUs$relabund <- CoastalInland_OTUs$mean/22354
#find average percent relative abundance for each of those OTUs
mean(CoastalInland_OTUs$relabund)*100



#Calculate number of reads in Inland:Bulk intersection
InlandBulk_OTUs <- subset(mim_otu, rownames(mim_otu) %in% attr(ItemsList,"intersections")$'InlandEco:Bulk')
sum(rowSums(InlandBulk_OTUs)) 
sum(rowSums(InlandBulk_OTUs))/sum(rowSums(mim_otu))*100
#remove columns that sum to zero. Then calc avg rel abund for each OTU.
InlandBulk_OTUs<- InlandBulk_OTUs[,colSums(InlandBulk_OTUs) > 0]
InlandBulk_OTUs$mean<- rowMeans(InlandBulk_OTUs)
#find relative abundance of each of those OTUs
InlandBulk_OTUs$relabund <- InlandBulk_OTUs$mean/22354
#find average percent relative abundance for each of those OTUs
mean(InlandBulk_OTUs$relabund)*100

#Calculate number of reads in Coastal:Bulk intersection
CoastalBulk_OTUs <- subset(mim_otu, rownames(mim_otu) %in% attr(ItemsList,"intersections")$'CoastalEco:Bulk')
sum(rowSums(CoastalBulk_OTUs)) 
sum(rowSums(CoastalBulk_OTUs))/sum(rowSums(mim_otu))*100
#remove columns that sum to zero. Then calc avg rel abund for each OTU.
CoastalBulk_OTUs<- CoastalBulk_OTUs[,colSums(CoastalBulk_OTUs) > 0]
CoastalBulk_OTUs$mean<- rowMeans(CoastalBulk_OTUs)
#find relative abundance of each of those OTUs
CoastalBulk_OTUs$relabund <- CoastalBulk_OTUs$mean/22354
#find average percent relative abundance for each of those OTUs
mean(CoastalBulk_OTUs$relabund)*100


#Calculate number of reads in Coastal:Inland:Bulk intersection
CoastalInlandBulk_OTUs <- subset(mim_otu, rownames(mim_otu) %in% attr(ItemsList,"intersections")$'CoastalEco:InlandEco:Bulk')
sum(rowSums(CoastalInlandBulk_OTUs)) 
sum(rowSums(CoastalInlandBulk_OTUs))/sum(rowSums(mim_otu))*100
#remove columns that sum to zero. Then calc avg rel abund for each OTU.
CoastalInlandBulk_OTUs<- CoastalInlandBulk_OTUs[,colSums(CoastalInlandBulk_OTUs) > 0]
CoastalInlandBulk_OTUs$mean<- rowMeans(CoastalInlandBulk_OTUs)
#find relative abundance of each of those OTUs
CoastalInlandBulk_OTUs$relabund <- CoastalInlandBulk_OTUs$mean/22354
#find average percent relative abundance for each of those OTUs
mean(CoastalInlandBulk_OTUs$relabund)*100

#make vectors of labels in order to manually change labels of Venn Diagram.
Coastal <- c('741 OTUs',paste(round((mean(CoastalOnly_OTUs$relabund)*100),digits=5),'%',sep=''))

Coastal <- c((paste((length(attr(ItemsList,"intersections")$'CoastalEco')),'OTUs',sep=' ')),paste(round((mean(CoastalOnly_OTUs$relabund)*100),digits=5),'%',sep=''))
Inland <- c((paste((length(attr(ItemsList,"intersections")$'InlandEco')),'OTUs',sep=' ')),paste(round((mean(InlandOnly_OTUs$relabund)*100),digits=5),'%',sep=''))
Bulk <- c((paste((length(attr(ItemsList,"intersections")$'Bulk')),'OTUs',sep=' ')),paste(round((mean(BulkOnly_OTUs$relabund)*100),digits=5),'%',sep=''))
CoastalInland <- c((paste((length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco')),'OTUs',sep=' ')),paste(round((mean(CoastalInland_OTUs$relabund)*100),digits=5),'%',sep=''))
InlandBulk <- c((paste((length(attr(ItemsList,"intersections")$'InlandEco:Bulk')),'OTUs',sep=' ')),paste(round((mean(InlandBulk_OTUs$relabund)*100),digits=5),'%',sep=''))
CoastalBulk <- c((paste((length(attr(ItemsList,"intersections")$'CoastalEco:Bulk')),'OTUs',sep=' ')),paste(round((mean(CoastalBulk_OTUs$relabund)*100),digits=5),'%',sep=''))
CoastalInlandBulk <- c((paste((length(attr(ItemsList,"intersections")$'CoastalEco:InlandEco:Bulk')),'OTUs',sep=' ')),paste(round((mean(CoastalInlandBulk_OTUs$relabund)*100),digits=5),'%',sep=''))


# have a look at the names in the plot object v
lapply(v,  names)
# We are interested in the labels
lapply(v, function(i) i$label)

#Manually change labels of venn diagram based on vectors created above.
v[[7]]$label  <- paste(Coastal, collapse="\n")
v[[8]]$label  <- paste(CoastalInland, collapse="\n")
v[[9]]$label  <- paste(Inland, collapse="\n")
v[[10]]$label  <- paste(CoastalBulk, collapse="\n")
v[[11]]$label  <- paste(CoastalInlandBulk, collapse="\n")
v[[12]]$label  <- paste(InlandBulk, collapse="\n")
v[[13]]$label  <- paste(Bulk, collapse="\n")

lapply(v, function(i) i$label)
library(VennDiagram)
grid.newpage()
grid.draw(v)


#https://stackoverflow.com/questions/25019794/venn-diagram-with-item-labels/45906498
