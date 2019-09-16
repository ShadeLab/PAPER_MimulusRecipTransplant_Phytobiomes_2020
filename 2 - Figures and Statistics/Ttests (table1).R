## Change in geochemistry between sites

library(ggplot2)
library(reshape2)
mim_envir<-read.table("Mimulus_metadata_FULL.txt", header=T, row.names=1, sep='\t')

#trim to variables of interest
mim_envir <- subset(mim_envir, select=c("Rhizo_Bulk", "Site_planted", "pH",	"Phosphorus",	"Potassium",	"Calcium",	"Magnesium",	"Copper",	"Percent_Organic_Matter",	"Sodium",	"Nitrate",	"Ammonium",	"Percent_Moisture",	"TotalN",	"Sulfur"))

#retain only bulk soil samples
mim_envir<-subset(mim_envir, Rhizo_Bulk =="Bulk")

mimulus_m<-melt(mim_envir)
str(mimulus_m)


#summary stats for table
library(plyr)
mim_geo_sum<-ddply(mimulus_m, c("Site_planted", "variable"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))
mim_geo_sum
write.csv(mim_geo_sum, "mim_geo_sum.csv")  


#statistical tests

#Ttests (and tests of homogeneity of variance)
bartlett.test(pH ~ Site_planted, data=mim_envir)
library(car)
leveneTest(pH ~ Site_planted, data=mim_envir)
t.test(pH ~ Site_planted, data=mim_envir)
wilcox.test(pH~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$pH, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$pH, mim_envir$Site_planted=="Inland"))


bartlett.test(Phosphorus ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Phosphorus ~ Site_planted, data=mim_envir) 
t.test(Phosphorus ~ Site_planted, data=mim_envir, var.equal=FALSE)
wilcox.test(Phosphorus ~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Phosphorus, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Phosphorus, mim_envir$Site_planted=="Inland"))

bartlett.test(Potassium ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Potassium ~ Site_planted, data=mim_envir)
t.test(Potassium ~ Site_planted, data=mim_envir, var.equal=TRUE)
wilcox.test(Potassium~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Potassium, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Potassium, mim_envir$Site_planted=="Inland"))

bartlett.test(Calcium ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Calcium ~ Site_planted, data=mim_envir)
t.test(Calcium ~ Site_planted, data=mim_envir)
wilcox.test(Calcium~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Calcium, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Calcium, mim_envir$Site_planted=="Inland"))

bartlett.test(Magnesium ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Magnesium ~ Site_planted, data=mim_envir)
t.test(Magnesium ~ Site_planted, data=mim_envir, var.equal=FALSE)
wilcox.test(Magnesium~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Magnesium, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Magnesium, mim_envir$Site_planted=="Inland"))

bartlett.test(Copper ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Copper ~ Site_planted, data=mim_envir)
t.test(Copper ~ Site_planted, data=mim_envir, var.equal=FALSE)
wilcox.test(Copper~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Copper, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Copper, mim_envir$Site_planted=="Inland"))

bartlett.test(Percent_Organic_Matter ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Percent_Organic_Matter ~ Site_planted, data=mim_envir)
t.test(Percent_Organic_Matter ~ Site_planted, data=mim_envir, var.equal=TRUE)
wilcox.test(Percent_Organic_Matter~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Percent_Organic_Matter, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Percent_Organic_Matter, mim_envir$Site_planted=="Inland"))

bartlett.test(Sodium ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Sodium ~ Site_planted, data=mim_envir)
t.test(Sodium ~ Site_planted, data=mim_envir, var.equal=FALSE)
wilcox.test(Sodium~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Sodium, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Sodium, mim_envir$Site_planted=="Inland"))

bartlett.test(Nitrate ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Nitrate ~ Site_planted, data=mim_envir)
t.test(Nitrate ~ Site_planted, data=mim_envir, var.equal=FALSE)
wilcox.test(Nitrate~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Nitrate, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Nitrate, mim_envir$Site_planted=="Inland"))


bartlett.test(Ammonium ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Ammonium ~ Site_planted, data=mim_envir)
t.test(Ammonium ~ Site_planted, data=mim_envir)
wilcox.test(Ammonium~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Ammonium, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Ammonium, mim_envir$Site_planted=="Inland"))

bartlett.test(Percent_Moisture ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Percent_Moisture ~ Site_planted, data=mim_envir)
t.test(Percent_Moisture ~ Site_planted, data=mim_envir)
wilcox.test(Percent_Moisture~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Percent_Moisture, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Percent_Moisture, mim_envir$Site_planted=="Inland"))

bartlett.test(TotalN ~ Site_planted, data=mim_envir)
library(car)
leveneTest(TotalN ~ Site_planted, data=mim_envir)
t.test(TotalN ~ Site_planted, data=mim_envir)
wilcox.test(TotalN~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$TotalN, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$TotalN, mim_envir$Site_planted=="Inland"))


bartlett.test(Sulfur ~ Site_planted, data=mim_envir)
library(car)
leveneTest(Sulfur ~ Site_planted, data=mim_envir)
t.test(Sulfur ~ Site_planted, data=mim_envir, var.equal=TRUE)
wilcox.test(Sulfur~ Site_planted, paired=FALSE, data=mim_envir)
shapiro.test(subset(mim_envir$Sulfur, mim_envir$Site_planted=="Coastal"))
shapiro.test(subset(mim_envir$Sulfur, mim_envir$Site_planted=="Inland"))
