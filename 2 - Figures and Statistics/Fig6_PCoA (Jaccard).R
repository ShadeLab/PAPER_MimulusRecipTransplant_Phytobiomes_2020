library(vegan)

#import otu table (remove "#Constructed from biom file " from first line before importing).
mim_otu<-read.table("zotu_table_forR.txt", row.names=1, header=T)

#import environmental data (remove "#" from first line before importing).
#remove "BarcodeSequence" column before importing.
#remove "LinkerPrimerSequence" column before importing.
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, sep='\t', row.names=1)



#convert otu table to data matrix
mim_otu <- data.matrix(t(mim_otu))
jac_dist <- vegdist(mim_otu, method="jaccard", binary=TRUE) 
jcpcoa <- cmdscale(jac_dist, eig=TRUE)
explainvar1 <- round(jcpcoa$eig[1] / sum(jcpcoa$eig), 3) * 100
explainvar2 <- round(jcpcoa$eig[2] / sum(jcpcoa$eig), 3) * 100

#convert scores to dataframe

mim.pcoa<- data.frame(jcpcoa$points[,1], jcpcoa$points[,2])  #axis 2 scores
row.names(mim.pcoa)<-row.names(jcpcoa$points)
colnames(mim.pcoa) <- c("Jac_Axis1","Jac_Axis2")

#also, add scores to metadata for plotting
rownames(mim_envir)==rownames(mim.pcoa)
mim_envir <- cbind(mim_envir,mim.pcoa)


#reduce environmental data to variables of interest
mim_envir_trimmed <- subset(mim_envir, select=c("pH",	"Phosphorus",	"Potassium",	"Calcium",	"Magnesium",	"Copper",	"Percent_Organic_Matter",	"Sodium",	"Nitrate",	"Ammonium",	"Percent_Moisture",	"TotalN",	"Sulfur"))
#replace column headers with abbreviations
names(mim_envir_trimmed) <- c("pH",	"P",	"K",	"Ca",	"Mg",	"Cu",	"OM",	"Na",	"NO3",	"NH4",	"PercMois",	"N",	"S")

#fit environment data to PCoA axis scores
mim.envfit<-envfit(mim.pcoa, mim_envir_trimmed)
mim.envfit

#environment scores (vectors scaled by R2 values)
mim.scores<-as.data.frame(scores(mim.envfit, display="vectors"))
mim.scores <- cbind(mim.scores, Variable = rownames(mim.scores))

#This is how to calculate the multiplier
#but it doesn't make the graph look any better.
#mult <-vegan:::ordiArrowMul(result_of_envfit)

mult <-.15

#Remove variables that were not significantly correlated with axis1 or axis2
list <- c("NH4")
mim.scores<-mim.scores[!(row.names(mim.scores) %in% list), ]





#split mim_envir into coastal site and inland site
mim_envir_Inland<-subset(mim_envir, Site_planted =="Inland")
mim_envir_Coastal<-subset(mim_envir,Site_planted =="Coastal")



library(ggplot2)
Fig6<-ggplot()+
 theme_bw()+
 geom_point(aes(mim_envir_Inland$Jac_Axis1, y=mim_envir_Inland$Jac_Axis2,color=mim_envir_Inland$SampleType),size=2.5,shape=21,stroke=1.5)+
 geom_point(aes(mim_envir_Coastal$Jac_Axis1,y=mim_envir_Coastal$Jac_Axis2,fill=mim_envir_Coastal$SampleType,color=mim_envir_Coastal$SampleType),size=2.5,shape=21,stroke=1.5)+  
 xlab(paste("PC1-", explainvar1,"%", sep=""))+
 ylab(paste("PC2-",explainvar2,"%",sep=""))+
 coord_fixed()+
 scale_color_manual(name="Inland site", breaks=c("Bulk Soil","LMC (Inland ecotype)","OCC (Inland ecotype)","MRR (Coastal ecotype)","SWB (Coastal ecotype)"),values=c("gray10","green","royalblue1","green4","royalblue4"))+
 scale_fill_manual(name="Coastal site", breaks=c("Bulk Soil","LMC (Inland ecotype)","OCC (Inland ecotype)","MRR (Coastal ecotype)","SWB (Coastal ecotype)"),values=c("gray10","green","royalblue1","green4","royalblue4"))+
 geom_segment(data=mim.scores, aes(x=0, xend=mult*Jac_Axis1, y=0, yend=mult*Jac_Axis2), arrow = arrow(length = unit(0.3, "cm")), colour = "black")+
 geom_text(data = mim.scores, aes(x = mult*Jac_Axis1, y = mult*Jac_Axis2, label = Variable), size = 2.5,fontface="bold",position=position_dodge2(width=0.15))+
 theme(text = element_text(size=8),legend.title=element_text(size=7))+
 guides(fill = guide_legend(override.aes = list(colour = c("gray10","green","green4","royalblue1","royalblue4"))))
Fig6   


ggsave("Figure6.tiff",Fig6,width=6, height=4, units="in",dpi=300)

