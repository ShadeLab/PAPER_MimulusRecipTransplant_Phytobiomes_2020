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
  ylim(2750,4700)+
  scale_fill_manual(values=c("gray10","royalblue1","green"))+
  geom_boxplot()+
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
  geom_boxplot()+
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
  geom_boxplot()+
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



mim_diversity_CoastalSite <- subset(mim_diversity, Site_planted =="Coastal")

###Coastal site, comparing bulk soil and inland ecotype
mim_diversity_CoastalSite_InlandandBulk <- subset(mim_diversity_CoastalSite, Origin=="Inland"|Rhizo_Bulk=="Bulk")
t.test(Richness~Rhizo_Bulk, paired=FALSE, var.equal = TRUE, data=mim_diversity_CoastalSite_InlandandBulk)
wilcox.test(Richness~Rhizo_Bulk, data=mim_diversity_CoastalSite_InlandandBulk)
t.test(PD_whole_tree_alpha~Rhizo_Bulk,paired=FALSE, var.equal = TRUE, data=mim_diversity_CoastalSite_InlandandBulk)
wilcox.test(PD_whole_tree_alpha~Rhizo_Bulk,data=mim_diversity_CoastalSite_InlandandBulk)
t.test(shannon_alpha~Rhizo_Bulk,paired=FALSE, var.equal = TRUE, data=mim_diversity_CoastalSite_InlandandBulk)
wilcox.test(shannon_alpha~Rhizo_Bulk,data=mim_diversity_CoastalSite_InlandandBulk)


###Coastal site, comparing bulk soil and coastal ecotype
mim_diversity_CoastalSite_CoastalandBulk <- subset(mim_diversity_CoastalSite, Origin=="Coastal"|Rhizo_Bulk=="Bulk")
t.test(Richness~Rhizo_Bulk, paired=FALSE, var.equal = TRUE, data=mim_diversity_CoastalSite_CoastalandBulk)
kruskal.test(Richness~Rhizo_Bulk, data=mim_diversity_CoastalSite_CoastalandBulk)
t.test(PD_whole_tree_alpha~Rhizo_Bulk,paired=FALSE, var.equal = TRUE, data=mim_diversity_CoastalSite_CoastalandBulk)
wilcox.test(PD_whole_tree_alpha~Rhizo_Bulk,data=mim_diversity_CoastalSite_CoastalandBulk)
t.test(shannon_alpha~Rhizo_Bulk,paired=FALSE, var.equal = TRUE, data=mim_diversity_CoastalSite_CoastalandBulk)
wilcox.test(shannon_alpha~Rhizo_Bulk,data=mim_diversity_CoastalSite_CoastalandBulk)



mim_diversity_InlandSite <- subset(mim_diversity, Site_planted =="Inland")

###Inland site, comparing bulk soil and inland ecotype
mim_diversity_InlandSite_InlandandBulk <- subset(mim_diversity_InlandSite, Origin=="Inland"|Rhizo_Bulk=="Bulk")
t.test(Richness~Rhizo_Bulk, paired=FALSE, var.equal = TRUE, data=mim_diversity_InlandSite_InlandandBulk)
kruskal.test(Richness~Rhizo_Bulk, data=mim_diversity_InlandSite_InlandandBulk)
t.test(PD_whole_tree_alpha~Rhizo_Bulk,paired=FALSE, var.equal = TRUE, data=mim_diversity_InlandSite_InlandandBulk)
wilcox.test(PD_whole_tree_alpha~Rhizo_Bulk,data=mim_diversity_InlandSite_InlandandBulk)
t.test(shannon_alpha~Rhizo_Bulk,paired=FALSE, var.equal = TRUE, data=mim_diversity_InlandSite_InlandandBulk)
wilcox.test(shannon_alpha~Rhizo_Bulk,data=mim_diversity_InlandSite_InlandandBulk)


###Inland site, comparing bulk soil and coastal ecotype
mim_diversity_InlandSite_CoastalandBulk <- subset(mim_diversity_InlandSite, Origin=="Coastal"|Rhizo_Bulk=="Bulk")
t.test(Richness~Rhizo_Bulk, paired=FALSE, var.equal = TRUE, data=mim_diversity_InlandSite_CoastalandBulk)
kruskal.test(Richness~Rhizo_Bulk, data=mim_diversity_InlandSite_CoastalandBulk)
t.test(PD_whole_tree_alpha~Rhizo_Bulk,paired=FALSE, var.equal = TRUE, data=mim_diversity_InlandSite_CoastalandBulk)
wilcox.test(PD_whole_tree_alpha~Rhizo_Bulk,data=mim_diversity_InlandSite_CoastalandBulk)
t.test(shannon_alpha~Rhizo_Bulk,paired=FALSE, var.equal = TRUE, data=mim_diversity_InlandSite_CoastalandBulk)
wilcox.test(shannon_alpha~Rhizo_Bulk,data=mim_diversity_InlandSite_CoastalandBulk)





####
#Coastal site, comparing Inland and Coastal Ecotypes
mim_diversity_CoastalSite <- subset(mim_diversity, Site_planted =="Coastal")
mim_diversity_CoastalSite_Rhizosphere <- subset(mim_diversity_CoastalSite, Rhizo_Bulk=="Rhizo")

t.test(Richness~Origin, paired=FALSE, var.equal = TRUE, data=mim_diversity_CoastalSite_Rhizosphere)
kruskal.test(Richness~Origin, data=mim_diversity_CoastalSite_Rhizosphere)
t.test(PD_whole_tree_alpha~Origin,paired=FALSE, var.equal = TRUE, data=mim_diversity_CoastalSite_Rhizosphere)
wilcox.test(PD_whole_tree_alpha~Origin,data=mim_diversity_CoastalSite_Rhizosphere)
t.test(shannon_alpha~Origin,paired=FALSE, var.equal = TRUE, data=mim_diversity_CoastalSite_Rhizosphere)
wilcox.test(shannon_alpha~Origin,data=mim_diversity_CoastalSite_Rhizosphere)


#Inland site, comparing Inland and Coastal Ecotypes
mim_diversity_InlandSite <- subset(mim_diversity, Site_planted =="Inland")
mim_diversity_InlandSite_Rhizosphere <- subset(mim_diversity_InlandSite, Rhizo_Bulk=="Rhizo")

t.test(Richness~Origin, paired=FALSE, var.equal = FALSE, data=mim_diversity_InlandSite_Rhizosphere)
wilcox.test(Richness~Origin, data=mim_diversity_InlandSite_Rhizosphere)
t.test(PD_whole_tree_alpha~Origin,paired=FALSE, var.equal = TRUE, data=mim_diversity_InlandSite_Rhizosphere)
wilcox.test(PD_whole_tree_alpha~Origin,data=mim_diversity_InlandSite_Rhizosphere)
t.test(shannon_alpha~Origin,paired=FALSE, var.equal = FALSE, data=mim_diversity_InlandSite_Rhizosphere)
wilcox.test(shannon_alpha~Origin,data=mim_diversity_InlandSite_Rhizosphere)


bartlett.test(Richness ~Origin, data=mim_diversity_InlandSite_Rhizosphere)
library(car)
leveneTest(Richness ~Origin, data=mim_diversity_InlandSite_Rhizosphere)

bartlett.test(PD_whole_tree_alpha ~Origin, data=mim_diversity_InlandSite_Rhizosphere)
library(car)
leveneTest(PD_whole_tree_alpha ~Origin, data=mim_diversity_InlandSite_Rhizosphere)

bartlett.test(shannon_alpha ~Origin, data=mim_diversity_InlandSite_Rhizosphere)
library(car)
leveneTest(shannon_alpha ~Origin, data=mim_diversity_InlandSite_Rhizosphere)
