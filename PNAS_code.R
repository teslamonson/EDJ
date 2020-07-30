###PLEASE CITE OUR MANUSCRIPT IF YOU USE THIS CODE
####T.M. MONSON, D. FECKER, M. SCHERRER. NEUTRAL EVOLUTION OF HUMAN ENAMEL-DENTINE JUNCTION MORPHOLOGY. PROCEEDINGS OF THE NATIONAL ACADEMY OF SCIENCES U.S.A. (2020)

#packages to run the analysis (multicore not needed on the laptop)
library(rgl)
library(car)
library(mgcv)
library(MASS)
library(pls)
library(abind)
library(Morpho)
library(geomorph)
library(multicore)
library(plyr)
library(devtools)
library(shapes)
library(sf)

##PROCRUSTES
# putting NTS files in a format you can run the analysis
lm_list_LLM1<-list.files()

#combine to dataframe
lm_LLM1<-ldply(lm_list_LLM1,read.table,skip=3)

##make it a matrix
lm_LLM1<-data.matrix(lm_LLM1, rownames.force = NA)

##array
lm_LLM1<-arrayspecs(lm_LLM1,104,3,sep=".")

#Run Procrustes Superimposition
symproc<-procSym(lm_LLM1,SMvector=c(1:9),outlines=list(c(1,10:15,2),c(2,16:21,3),c(3,22:32,4),c(4,33:38,5),c(5,39:44,6),c(6,45:55,7),c(8,56:104,9)),deselect=T,bending=TRUE,pairedLM=NULL, sizeshape=TRUE)

###RIJ PHENETIC MATRIX
##extract coordinates for Rij phenetic matrix
regionwise.coords=coords.subset(symproc$rotated,continent)
region_means=lapply(regionwise.coords,mshape)
avg<-mean(region_means$Australia)

###REGRESSION
qplot(GenRij, Rijmshape_100e10, data=rat)+theme_bw()+geom_point(size = 5)+geom_abline(intercept=839.5, slope=72051.6,color="black",linetype="dashed", size=0.5)+ geom_smooth(method="lm")+geom_text(aes(label=PopPair),hjust=0.5, vjust=2)
rcorr(rat$GenRij,rat$RijPC1_new,)
m<-lm(Rijmshape_100e10~GenRij, data=rat)
summary(m)

##PLOT PCS
###after extracting PCs from symproc, plot PCs with convex hulls
p<-qplot(PC1,PC2,  data=plot2, col=Continent,shape=factor(Continent), label=Specimen)+ scale_colour_manual(name="Continent", values=c("dodgerblue","red","orange","black","deeppink","turquoise1","gray","olivedrab3","brown","forestgreen","rosybrown3","darkmagenta","dodgerblue","cyan3", "yellow"))+ scale_shape_manual(name="Continent", values=c(3,13,15,16,17,18,20,21))+theme_bw()+geom_point(size = 8)
p + aes(fill = factor(Continent)) + geom_polygon(data = hull_cyl, alpha = 0.1)

###ANOVA
library(geomorph)
##FACTORS
region<-data$Final.Region
continent<-factor(data$Continent)
ID<-factor(data$ID)
gdf<-geomorph.data.frame(coords=symproc$rotated,site=continent)
anova <- procD.lm(coords~site, data = gdf, iter = 999, RRPP = TRUE)

####CAC ALLOMETRY
cac<-CAC(symproc$PCscores, symproc$size, groups=NULL, log=TRUE)
cac$CACscores
summary(cac)
plot(cac$CACscores,cac$size)#plot scores against Centroid size
cor.test(cac$CACscores,cac$size)

####CV
regionwise.coords=coords.subset(symproc$rotated,continent)
region_means=lapply(regionwise.coords,mshape)
CV_Aus=sd(region_means$Australia, na.rm=TRUE)/mean(region_means$Australia, na.rm=TRUE)*100
CV_Eur=sd(region_means$Europe, na.rm=TRUE)/mean(region_means$Europe, na.rm=TRUE)*100
CV_Asia=sd(region_means$Asia, na.rm=TRUE)/mean(region_means$Asia, na.rm=TRUE)*100
CV_Africa=sd(region_means$Africa, na.rm=TRUE)/mean(region_means$Africa, na.rm=TRUE)*100
CV_Amer=sd(region_means$Americas, na.rm=TRUE)/mean(region_means$Americas, na.rm=TRUE)*100
CV_SAf=sd(region_means$'South Africa', na.rm=TRUE)/mean(region_means$'South Africa', na.rm=TRUE)*100

###PROCRUSTES DISTANCES
procdist(region_means$Australia, region_means$'South Africa')
procdist(region_means$Europe, region_means$'South Africa')
procdist(region_means$Asia, region_means$'South Africa')
procdist(region_means$Americas, region_means$'South Africa')
procdist(region_means$Africa, region_means$'South Africa')