#set the work routine
setwd("~/Dropbox/Paper_III_MMB/Without_taxa_merging/Delete")
# This script was written for the final checking of the result
# The starting file was MMB_merge_derep_no0.csv.
# This file has been:
# 1. NXT, TSQ OTU table merged
# 2. TSQ OTU table dereplicated
# 3. The similarity lower than 97% was considered to be same species
# 4. The order for PM extraction was considered

################################ Prepare dataset ready for further organis ####################################
MMB_merge_count_pre<-read.csv("MMB_merge_derep_no0.csv")　

# merge notablis
MMB_merge_count_pre$TAXA<- gsub("k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis_L1",
                                "k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis",
                                MMB_merge_count_pre$TAXA)
MMB_merge_count_pre$TAXA<- gsub("k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis_L2",
                                "k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis",
                                MMB_merge_count_pre$TAXA)
MMB_merge_count_pre$TAXA<- gsub("k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis_L4",
                                "k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis",
                                MMB_merge_count_pre$TAXA)

MMB_merge_count_pre<- MMB_merge_count_pre[,2:106]

# Give column name to the table
sample<-replicate(104,"PB")

sample[1:24]<-replicate(24,"PB_NXT")
sample[25:44]<-replicate(20,"PM_NXT")
sample[45:52]<-replicate(8,"SP_NXT")
sample[53:76]<-replicate(24,"PB_TSQ")
sample[77:96]<-replicate(20,"PM_TSQ")
sample[97:104]<-replicate(8,"SP_TSQ")
a<-sample
sample<- replicate(105,"Taxa")
sample[2:105]<-a

names(MMB_merge_count_pre) = sample

#Merge OTU with the same taxa
MMB_merge_count<-aggregate(MMB_merge_count_pre[,2:105], by=list(Taxa=MMB_merge_count_pre$Taxa), FUN=sum)

# remove OTU total reads lower than 5
MMB_merge_count <- MMB_merge_count[rowSums(MMB_merge_count[,2:105])>5,]

# remove collembolan_unassigned
MMB_merge_count <- MMB_merge_count[-which(MMB_merge_count$Taxa=="k:Animilia,p:Arthropoda,c:Collembola"),]
# remove mite_unassigned
MMB_merge_count <- MMB_merge_count[-which(MMB_merge_count$Taxa=="k:Animilia,p:Arthropoda,c:Arachnida"),]
# remove F.candida
MMB_merge_count <- MMB_merge_count[-which(MMB_merge_count$Taxa=="k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Folsomia,s:Folsomia_candida"),]

#### calculate relative abundace for each OTUs ####
library(vegan)
MMB_merge_rltv<-decostand(MMB_merge_count[,2:105],"total",2)
# decostand function for standarization of ecology community data
# Method = total to calculate the relative abundance
# check<-colSums(MMB_taxa_rltv) # double check the calculation
# write.csv(cbind(TAXA,MMB_taxa_rltv),"MMB_taxa_relative.csv")

MMB_merge_rltv<-cbind(MMB_merge_count[,1],MMB_merge_rltv)

# Prepare dataset for further analysis ----------------------------------

#### readcounts, extract_method level ######

PB_NXT_count<-MMB_merge_count[,2:25]
PM_NXT_count<-MMB_merge_count[,30:45] # exclude 10 g extraction
SP_NXT_count<-MMB_merge_count[,46:53]
PB_TSQ_count<-MMB_merge_count[,54:77]
PM_TSQ_count<-MMB_merge_count[,82:97] # exclude 10 g extraction
SP_TSQ_count<-MMB_merge_count[,98:105]

MMB_merge_count_extract<-rbind(rowSums(PB_NXT_count),
                               rowSums(PM_NXT_count),
                               rowSums(SP_NXT_count),
                               rowSums(PB_TSQ_count),
                               rowSums(PM_TSQ_count),
                               rowSums(SP_TSQ_count))

MMB_merge_count_extract_rarefy<- decostand(MMB_merge_count_extract,"pa",2)
extract<-c("PB_NXT","PM_NXT","SP_NXT","PB_TSQ","PM_TSQ","SP_TSQ")
row.names(MMB_merge_count_extract_rarefy) <- extract

#### readcounts, soil core level ######
##### NXT
## PB_NXT_count
PB_NXT_count_derep =matrix(nrow=81,ncol=8)  #creat an empty matrix for PB
for ( i in 1:8)
{
  m<-(i-1)*3+1
  PB_NXT_count_derep[,i]<- rowSums(PB_NXT_count[,m:(m+2)])
}

## PM_NXT_count
PM_NXT_count_derep =matrix(nrow=81,ncol=8)  #creat an empty matrix for PM
for ( i in 1:8)
{
  m<-(i-1)*2+1
  PM_NXT_count_derep[,i]<- rowSums(PM_NXT_count[,m:(m+1)])
}

#### TSQ
##PB_TSQ_count
PB_TSQ_count_derep =matrix(nrow=81,ncol=8)  #creat an empty matrix for PB
for ( i in 1:8)
{
  m<-(i-1)*3+1
  PB_TSQ_count_derep[,i]<- rowSums(PB_TSQ_count[,m:(m+2)])
}

## PM_TSQ_count
PM_TSQ_count_derep =matrix(nrow=81,ncol=8)  #creat an empty matrix for PM
for ( i in 1:8)
{
  m<-(i-1)*2+1
  PM_TSQ_count_derep[,i]<- rowSums(PM_TSQ_count[,m:(m+1)])
}

PB_NXT_count_derep
PM_NXT_count_derep
SP_NXT_count
PB_TSQ_count_derep
PM_TSQ_count_derep
SP_TSQ_count

MMB_merge_count_soilcore <- cbind(MMB_merge_count$Taxa,PB_NXT_count_derep,PM_NXT_count_derep,SP_NXT_count,PB_TSQ_count_derep,PM_TSQ_count_derep,SP_TSQ_count)

write.csv(MMB_merge_count_soilcore,"check_the_list.csv")

sample_soilcore<-replicate(48,"PB")
sample_soilcore<-c(replicate(8,"PB_NXT"),replicate(8,"PM_NXT"),replicate(8,"SP_NXT"),replicate(8,"PB_TSQ"),replicate(8,"PM_TSQ"),replicate(8,"SP_TSQ"))

################################ dataset ready to use ####################################
#### Venn diagram ####

library(gplots)
#prepare dataset
nr <- nrow(MMB_merge_count)
PBL_NXT = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(PB_NXT_count[i,]) > 0) {PBL_NXT[i]= MMB_merge_count[i,1]}}
PBL_NXT <- PBL_NXT[PBL_NXT!="NA"]

PML_NXT = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(PM_NXT_count[i,]) > 0) {PML_NXT[i]= MMB_merge_count[i,1]}}
PML_NXT <- PML_NXT[PML_NXT!="NA"]

SPL_NXT = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(SP_NXT_count[i,]) > 0) {SPL_NXT[i]= MMB_merge_count[i,1]}}
SPL_NXT <- SPL_NXT[SPL_NXT!="NA"]

PBL_TSQ = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(PB_TSQ_count[i,]) > 0) {PBL_TSQ[i]= MMB_merge_count[i,1]}}
PBL_TSQ <- PBL_TSQ[PBL_TSQ!="NA"]

PML_TSQ = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(PM_TSQ_count[i,]) > 0) {PML_TSQ[i]= MMB_merge_count[i,1]}}
PML_TSQ <- PML_TSQ[PML_TSQ!="NA"]

SPL_TSQ = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(SP_TSQ_count[i,]) > 0) {SPL_TSQ[i]= MMB_merge_count[i,1]}}
SPL_TSQ <- SPL_TSQ[SPL_TSQ!="NA"]

##### draw graph #####
# After investigation again, trip scaled Venn diagrame R procedure is not existin right now, 20170209
# Use BioVenn instead
# This Venn diagramme should:
# 1: remove collembolan_unassigned mite_unassigned,
# 2: remove F.candida (Considered as a contamination from our lab)
# 3: combine 3 strains of notablis
input_1 <- list (SPL_NXT,SPL_TSQ)
venn(input_1)

input_2 <- list (PBL_TSQ,PML_TSQ,SPL_TSQ)
venn(input_2)

input_3 <-list (PBL_NXT,PML_NXT,SPL_NXT)
venn(input_3)

input_4<-list(PBL_NXT,PML_NXT,SPL_NXT,PBL_TSQ,PML_TSQ,SPL_TSQ)

#### individual based rarefaction curve #####################

# The rarefaction curve produced from this part is based the individuals (abundance data)
library(vegan)

# dataset
# subset for each method
# use counts for each method
tPB_NXT_count_derep<-t(PB_NXT_count_derep)
tPM_NXT_count_derep<-t(PM_NXT_count_derep) 
tSP_NXT_count<-t(SP_NXT_count)
tPB_TSQ_count_derep<-t(PB_TSQ_count_derep) 
tPM_TSQ_count_derep<-t(PM_TSQ_count_derep) 
tSP_TSQ_count<-t(SP_TSQ_count)

## PB_NXT accumulation 
S_PB_NXT <- specnumber(tPB_NXT_count_derep) 
raremax_PB_NXT<-min(rowSums(tPB_NXT_count_derep))
Srare_PB_NXT <- rarefy(tPB_NXT_count_derep, raremax_PB_NXT)

plot(S_PB_NXT, Srare_PB_NXT, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

abline(0, 1)

## PM_NXT accumulation 
S_PM_NXT <- specnumber(tPM_NXT_count_derep) 
raremax_PM_NXT<-min(rowSums(tPM_NXT_count_derep))
Srare_PM_NXT <- rarefy(tPM_NXT_count_derep, raremax_PM_NXT)

plot(S_PM_NXT, Srare_PM_NXT, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

abline(0, 1)

## SP_NXT accumulation
S_SP_NXT <- specnumber(tSP_NXT_count) 
raremax_SP_NXT<-min(rowSums(tSP_NXT_count))
Srare_SP_NXT <- rarefy(tSP_NXT_count, raremax_SP_NXT)

plot(S_SP_NXT, Srare_SP_NXT, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

abline(0, 1)

## PB_TSQ accumulation 
S_PB_TSQ <- specnumber(round(tPB_TSQ_count_derep)) 
raremax_PB_TSQ<-min(rowSums(round(tPB_TSQ_count_derep)))
Srare_PB_TSQ <- rarefy(round(tPB_TSQ_count_derep), raremax_PB_TSQ)

plot(S_PB_TSQ, Srare_PB_TSQ, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species",col="blue")

abline(0, 1,col="blue")

## PM_TSQ accumulation 
S_PM_TSQ <- specnumber(round(tPM_TSQ_count_derep)) 
raremax_PM_TSQ<-min(rowSums(round(tPM_TSQ_count_derep)))
Srare_PM_TSQ <- rarefy(round(tPM_TSQ_count_derep), raremax_PM_TSQ)

plot(S_PM_TSQ, Srare_PM_TSQ, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species",col="red")

abline(0, 1,col="red")

## SP_TSQ accumulation
S_SP_TSQ <- specnumber(round(tSP_TSQ_count)) 
raremax_SP_TSQ<-min(rowSums(round(tSP_TSQ_count)))
Srare_SP_TSQ <- rarefy(round(tSP_TSQ_count), raremax_SP_TSQ)

plot(S_SP_TSQ, Srare_SP_TSQ, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species",col="orange")

abline(0, 1,col="orange")

# export individual based rarefaction curve
par(cex = 0.8)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
par(tcl = -0.25)

par(mfrow=c(2,3))
rarecurve(tPB_NXT_count_derep, 
          sample = raremax_PB_NXT, 
          col = "blue", 
          xlim=c(0,40000),
          ylim=c(0,60),
          xaxt="n",
          label = FALSE,
          cex = 0.6)
mtext("Nextera", font=2,side=2,line=2.5,cex =1)
mtext("Phosphate Buffer",font=2,side=3,line=1.5,cex =1)
rarecurve(tPM_NXT_count_derep, 
          sample = raremax_PM_NXT, 
          col = "red", 
          xlim=c(0,40000),
          ylim=c(0,60),
          label = FALSE,
          xaxt="n",
          yaxt="n",
          cex = 0.6)
mtext("Powermax",font=2,side=3,line=1.5,cex =1)
rarecurve(tSP_NXT_count, 
          sample = raremax_SP_NXT, 
          col = "orange", 
          xlim=c(0,40000),
          ylim=c(0,60),
          label = FALSE,
          xaxt="n",
          yaxt="n",
          cex = 0.6)
mtext("DNA soup",font=2,side=3,line=1.5,cex =1)

rarecurve(round(tPB_TSQ_count_derep), 
          sample = raremax_PB_TSQ, 
          col = "blue", 
          xlim=c(0,40000),
          ylim=c(0,60),
          label = FALSE,
          cex = 0.6)
mtext("TruSeq", font=2,side=2,line=2.5,cex =1)
rarecurve(round(tPM_TSQ_count_derep), 
          sample = raremax_PM_TSQ, 
          col = "red", 
          xlim=c(0,40000),
          ylim=c(0,60),
          label = FALSE,
          yaxt="n",
          cex = 0.6)

rarecurve(round(tSP_TSQ_count), 
          sample = raremax_SP_TSQ, 
          col = "orange", 
          xlim=c(0,40000),
          ylim=c(0,60),
          label = FALSE,
          yaxt="n",
          cex = 0.6)

#### Sample based rarefaction curve #############
# Use vegan based package:rareNMtests to maximize the possibility of drawing rarefaction curve

library("rareNMtests")

# Ecological null model test using sample-based rarefaction curves

# for species richness (q = 0)

## sample based rarefaction curve for entire dataset #####

## sample based rarefaction with readcount in soil core level #########
#dataset
MMB_merge_count_soilcore
sample_soilcore
MMB_merge_count_soilcore_rafy<-decostand(MMB_merge_count_soilcore[,2:49],"pa",2)
sample_rare_soilcore<-cbind(sample_soilcore,t(MMB_merge_count_soilcore_rafy))

FacToCha <- function(x) {if(is.factor(x)) 
{as.character(paste(x))} else {x}} 
Library<-substr(sample_soilcore,4,6)
Extract<-substr(sample_soilcore,1,2)
#p=0.07
sbecoq_sample_rare_soilcore <- EcoTest.sample(sample_rare_soilcore[,-1], by=sample_rare_soilcore[,1], MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

# fix PB compare NXT and TSQ p=0.07
sbecoq_sample_rare_soilcore <- EcoTest.sample(rbind(sample_rare_soilcore[1:8,-1],sample_rare_soilcore[25:32,-1]), by=c(Library[1:8],Library[25:32]), MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

# fix PM compare NXT and TSQ p=0.125
sbecoq_sample_rare_soilcore <- EcoTest.sample(rbind(sample_rare_soilcore[9:16,-1],sample_rare_soilcore[33:40,-1]), by=c(Library[1:8],Library[25:32]), MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

# fix SP compare NXT and TSQ p=0.215
sbecoq_sample_rare_soilcore <- EcoTest.sample(rbind(sample_rare_soilcore[17:24,-1],sample_rare_soilcore[41:48,-1]), by=c(Library[1:8],Library[25:32]), MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

#Compare within NXT p=0.27
sbecoq_sample_rare_soilcore <- EcoTest.sample(sample_rare_soilcore[1:24,-1], by=Extract[1:24], MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

#Compare within TSQ p=0.06 
sbecoq_sample_rare_soilcore <- EcoTest.sample(sample_rare_soilcore[25:48,-1], by=Extract[25:48], MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

# library no difference regarding to the rarefaction curve
# extraction no difference regarding to the rarefaction cure as well

###### accumulation curve ######################
library(vegan)
# Accumulation curve

# Function specaccum finds species accumulation curves or the number of species for a certain number of sampled sites or individuals.

#dataset
PB_NXT_count_derep
PM_NXT_count_derep
SP_NXT_count
PB_TSQ_count_derep
PM_TSQ_count_derep
SP_TSQ_count

#alpha("blue", 0.1)
#adjustcolor( "blue", alpha.f = 0.2)
sp1 <- specaccum(t(PB_NXT_count_derep))
sp2 <- specaccum(t(PM_NXT_count_derep))
sp3 <- specaccum(t(SP_NXT_count))
sp4 <- specaccum(t(PB_TSQ_count_derep))
sp5 <- specaccum(t(PM_TSQ_count_derep))
sp6 <- specaccum(t(SP_TSQ_count))


#### draw accumulation curve -------------------------------------------------
par(mfrow=c(1,1))
par(mar = c(4, 4, 2, 0), oma = c(1, 1, 1, 0.5))
plot(sp1,  
     xlim=c(1,8),
     ylim=c(15,80),
     main="Accumulation curve",#main title of the graph
     xlab = "Number of samples",#Lable for x axis
     ylab = "Number of species",# Lable for y axis
     ci.type="poly", 
     col="blue", 
     lwd=2, 
     ci.lty=0, 
     ci.col=adjustcolor( "lightblue", alpha.f = 0.5) #adjust the transparency 
)
#boxplot(sp2, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp2, 
     xlim=c(1,8),
     ylim=c(15,80),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     axes=F,
     col="red", 
     lwd=2, 
     ci.lty=0, 
     ci.type="poly", 
     ci.col=adjustcolor( "pink", alpha.f = 0.5)
)
#boxplot(sp4, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp3, 
     xlim=c(1,8),
     ylim=c(15,80),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     axes=F,
     ci.type="poly", 
     col="orange", 
     lwd=2, 
     ci.lty=0, 
     ci.col=adjustcolor( "yellow", alpha.f = 0.5)
)
#boxplot(sp6, col="yellow", add=TRUE, pch="+")
legend(5,30, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       text.col = "black", 
       title="Nextera",
       lty = c(1, 1, 1),
       bty="n",lwd=2,cex=1)

par(mfrow=c(1,1))
plot(sp4,  
     xlim=c(1,8),
     ylim=c(15,80),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     #axes=F,
     ci.type="poly", 
     col="blue", 
     lwd=2, 
     lty=2,
     ci.lty=0, 
     ci.col=adjustcolor( "lightblue", alpha.f = 0.5) #adjust the transparency 
)
#boxplot(sp2, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp5, 
     xlim=c(1,8),
     ylim=c(15,80),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     axes=F,
     ci.type="poly", 
     col="red", 
     lwd=2, 
     lty=2,
     ci.lty=0, 
     ci.col=adjustcolor( "pink", alpha.f = 0.5)
)
#boxplot(sp4, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp6, 
     xlim=c(1,8),
     ylim=c(15,80),
     main="Accumulation curve",#main title of the graph
     xlab = "Number of samples",#Lable for x axis
     ylab = "Number of species",# Lable for y axis
     ci.type="poly", 
     col="orange", 
     lwd=2, 
     lty=2,
     ci.lty=0, 
     ci.col=adjustcolor( "yellow", alpha.f = 0.5)
)



legend(5,30, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       title="TruSeq",
       text.col = "black", 
       lty = c(2, 2, 2),
       bty="n",lwd=2,cex=1)


##### sample based rarefaction curve ############
sp1_rare <- specaccum(t(PB_NXT_count_derep),method="rarefaction")
sp2_rare <- specaccum(t(PM_NXT_count_derep),method="rarefaction")
sp3_rare <- specaccum(t(SP_NXT_count),method="rarefaction")
sp4_rare <- specaccum(t(round(PB_TSQ_count_derep)),method="rarefaction")
sp5_rare <- specaccum(t(round(PM_TSQ_count_derep)),method="rarefaction")
sp6_rare <- specaccum(t(round(SP_TSQ_count)),method="rarefaction")

# graphic
par(mfrow=c(1,1))
par(mar = c(4, 4, 2, 0), oma = c(1, 1, 1, 0.5))
plot(sp1_rare,  
     xlim=c(1,8),
     ylim=c(40,60),
     main="Rarefaction curve",#main title of the graph
     xlab = "Number of samples",#Lable for x axis
     ylab = "Number of species",# Lable for y axis
     ci.type="poly", 
     col="blue", 
     lwd=2, 
     ci.lty=0, 
     ci.col=adjustcolor( "lightblue", alpha.f = 0.5) #adjust the transparency 
)
#boxplot(sp2, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp2_rare, 
     xlim=c(1,8),
     ylim=c(40,60),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     axes=F,
     col="red", 
     lwd=2, 
     ci.lty=0, 
     ci.type="poly", 
     ci.col=adjustcolor( "pink", alpha.f = 0.5)
)
#boxplot(sp4, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp3_rare, 
     xlim=c(1,8),
     ylim=c(40,60),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     axes=F,
     ci.type="poly", 
     col="orange", 
     lwd=2, 
     ci.lty=0, 
     ci.col=adjustcolor( "yellow", alpha.f = 0.5)
)
#boxplot(sp6, col="yellow", add=TRUE, pch="+")
legend(6,45, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       text.col = "black", 
       title="Nextera",
       lty = c(1, 1, 1),
       bty="n",lwd=2,cex=1)

par(mfrow=c(1,1))
plot(sp4_rare,  
     xlim=c(1,8),
     ylim=c(40,60),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     #axes=F,
     ci.type="poly", 
     col="blue", 
     lwd=2, 
     lty=2,
     ci.lty=0, 
     ci.col=adjustcolor( "lightblue", alpha.f = 0.5) #adjust the transparency 
)
#boxplot(sp2, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp5_rare, 
     xlim=c(1,8),
     ylim=c(40,60),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     axes=F,
     ci.type="poly", 
     col="red", 
     lwd=2, 
     lty=2,
     ci.lty=0, 
     ci.col=adjustcolor( "pink", alpha.f = 0.5)
)
#boxplot(sp4, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp6_rare, 
     xlim=c(1,8),
     ylim=c(40,60),
     main="Rarefaction curve",#main title of the graph
     xlab = "Number of samples",#Lable for x axis
     ylab = "Number of species",# Lable for y axis
     ci.type="poly", 
     col="orange", 
     lwd=2, 
     lty=2,
     ci.lty=0, 
     ci.col=adjustcolor( "yellow", alpha.f = 0.5)
)

legend(6,45, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       title="TruSeq",
       text.col = "black", 
       lty = c(2, 2, 2),
       bty="n",lwd=2,cex=1)
#boxplot(sp6, col="yellow", add=TRUE, pch="+")

# done

#### Alpha diversity -- Chao2/ICE/Jackkinfe #####
# dataset
tPB_NXT_count_derep<-t(PB_NXT_count_derep)
tPM_NXT_count_derep<-t(PM_NXT_count_derep)
tSP_NXT_count<-t(SP_NXT_count)

tPB_TSQ_count_derep<-t(PB_TSQ_count_derep) 
tPM_TSQ_count_derep <-t(PM_TSQ_count_derep) 
tSP_TSQ_count <-t(SP_TSQ_count)

tPB_NXT_count_derep_rafy<-decostand(tPB_NXT_count_derep,"pa",2)
tPM_NXT_count_derep_rafy<-decostand(tPM_NXT_count_derep,"pa",2)
tSP_NXT_count_rafy<-decostand(tSP_NXT_count,"pa",2)

tPB_TSQ_count_derep_rafy<-decostand(tPB_TSQ_count_derep,"pa",2)
tPM_TSQ_count_derep_rafy<-decostand(tPM_TSQ_count_derep,"pa",2)
tSP_TSQ_count_rafy<-decostand(tSP_TSQ_count,"pa",2)


library(fossil)

Chao2_PB_TSQ<-chao2(tPB_TSQ_count_derep_rafy,taxa.row=F)
Chao2_PM_TSQ<-chao2(tPM_TSQ_count_derep_rafy,taxa.row=F)
Chao2_SP_TSQ<-chao2(tSP_TSQ_count_rafy,taxa.row=F)

alpha_fossil_PB_NXT<-spp.est(t(tPB_NXT_count_derep_rafy),abund=F)
alpha_fossil_PM_NXT<-spp.est(t(tPM_NXT_count_derep_rafy),abund=F)
alpha_fossil_SP_NXT<-spp.est(t(tSP_NXT_count_rafy),abund=F) 
alpha_fossil_PB_TSQ<-spp.est(t(tPB_TSQ_count_derep_rafy),abund=F)
alpha_fossil_PM_TSQ<-spp.est(t(tPM_TSQ_count_derep_rafy),abund=F)
alpha_fossil_SP_TSQ<-spp.est(t(tSP_TSQ_count_rafy),abund=F)

alpha_fossil_PB_NXT 
alpha_fossil_PM_NXT 
alpha_fossil_SP_NXT  
alpha_fossil_PB_TSQ 
alpha_fossil_PM_TSQ 
alpha_fossil_SP_TSQ 
# aov test
# Anova test is not propiate for incidence based diversity indexes

#### graph for CHAO2/ICE/Jackkinfe2 NXT####
par(cex = 0.8)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
par(tcl = -0.25)

par(mfrow=c(2,3))
plot(alpha_fossil_PB_NXT[,5]~alpha_fossil_PB_NXT[,1],
     ylim=c(10,100),
     xlim=c(1,8),
     xaxt="n",
     col="blue",
     type="l") #Chao2
polygon(c(alpha_fossil_PB_NXT[,1],rev(alpha_fossil_PB_NXT[,1])),
        c(alpha_fossil_PB_NXT[,6],rev(alpha_fossil_PB_NXT[,7])),border=NA,col=rgb(0, 0, 1,0.5))#blue
mtext("Nextera", font=2,side=2,line=2.5,cex =1)
mtext("Chao2",font=2,side=3,line=1.5,cex =1)

par(new=TRUE)
plot(alpha_fossil_PM_NXT[,5],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="red",
     type="l") #Chao2
polygon(c(alpha_fossil_PM_NXT[,1],rev(alpha_fossil_PM_NXT[,1])),
        c(alpha_fossil_PM_NXT[,6],rev(alpha_fossil_PM_NXT[,7])),border=NA,col=rgb(1, 0, 0,0.5))#red

par(new=TRUE)
plot(alpha_fossil_SP_NXT[,5],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="orange",
     type="l") #Chao2
polygon(c(alpha_fossil_SP_NXT[,1],rev(alpha_fossil_SP_NXT[,1])),
        c(alpha_fossil_SP_NXT[,6],rev(alpha_fossil_SP_NXT[,7])),border=NA,col=rgb(1, 1, 0,0.5))#orange

plot(alpha_fossil_PB_NXT[,8],
     ylim=c(10,100),
     xaxt="n",
     yaxt="n",
     col="blue",
     type="l") #ICE
polygon(c(alpha_fossil_PB_NXT[,1],rev(alpha_fossil_PB_NXT[,1])),
        c(alpha_fossil_PB_NXT[,9],rev(alpha_fossil_PB_NXT[,10])),border=NA,col=rgb(0, 0, 1,0.5))#red
mtext("ICE",font=2,side=3,line=1.5,cex =1)

par(new=TRUE)
plot(alpha_fossil_PM_NXT[,8],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="red",
     type="l") #ICE
polygon(c(alpha_fossil_PM_NXT[,1],rev(alpha_fossil_PM_NXT[,1])),
        c(alpha_fossil_PM_NXT[,9],rev(alpha_fossil_PM_NXT[,10])),border=NA,col=rgb(1, 0, 0,0.5))

par(new=TRUE)
plot(alpha_fossil_SP_NXT[,8],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="orange",
     type="l") #ICE
polygon(c(alpha_fossil_SP_NXT[,1],rev(alpha_fossil_SP_NXT[,1])),
        c(alpha_fossil_SP_NXT[,9],rev(alpha_fossil_SP_NXT[,10])),border=NA,col=rgb(1, 1, 0,0.5))

plot(alpha_fossil_PB_NXT[,11],
     ylim=c(10,100),
     col="blue",
     xaxt="n",
     yaxt="n",
     type="l") #Jacknife
polygon(c(alpha_fossil_PB_NXT[,1],rev(alpha_fossil_PB_NXT[,1])),
        c(alpha_fossil_PB_NXT[,12],rev(alpha_fossil_PB_NXT[,13])),border=NA,col=rgb(0, 0, 1,0.5))
mtext("Jackknife1",font=2,side=3,line=1.5,cex =1)

par(new=TRUE)  
plot(alpha_fossil_PM_NXT[,11],
     ylim=c(10,100),
     axes=F,
     xlab="",
     xaxt="n",
     yaxt="n",
     col="red",
     type="l") #Jacknife
polygon(c(alpha_fossil_PM_NXT[,1],rev(alpha_fossil_PM_NXT[,1])),
        c(alpha_fossil_PM_NXT[,12],rev(alpha_fossil_PM_NXT[,13])),border=NA,col=rgb(1, 0, 0,0.5))

par(new=TRUE)
plot(alpha_fossil_SP_NXT[,11],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     yaxt="n",
     xlab="",
     col="orange",
     type="l") #Jacknife
polygon(c(alpha_fossil_SP_NXT[,1],rev(alpha_fossil_SP_NXT[,1])),
        c(alpha_fossil_SP_NXT[,12],rev(alpha_fossil_SP_NXT[,13])),border=NA,col=rgb(1, 1, 0,0.5))

legend(3,20, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       text.col = "black",
       lty=c(1, 1, 1),
       bty="n",cex=1)

####graph for CHAO2/ICE/Jackkinfe2 TSQ #################

plot(alpha_fossil_PB_TSQ[,5]~alpha_fossil_PB_TSQ[,1],
     ylim=c(10,100),
     xlim=c(1,8),
     xlab="Chao2",
     col="blue",
     lty=2,
     type="l") #Chao2
polygon(c(alpha_fossil_PB_TSQ[,1],rev(alpha_fossil_PB_TSQ[,1])),
        c(alpha_fossil_PB_TSQ[,6],rev(alpha_fossil_PB_TSQ[,7])),border=NA,col=rgb(0, 0, 1,0.5))#blue
mtext("TruSeq", font=2,side=2,line=2.5,cex =1)

par(new=TRUE)
plot(alpha_fossil_PM_TSQ[,5],
     ylim=c(10,100),
     axes=F,
     xlab="",
     col="red",
     lty=2,
     type="l") #Chao2
polygon(c(alpha_fossil_PM_TSQ[,1],rev(alpha_fossil_PM_TSQ[,1])),
        c(alpha_fossil_PM_TSQ[,6],rev(alpha_fossil_PM_TSQ[,7])),border=NA,col=rgb(1, 0, 0,0.5))#red

par(new=TRUE)
plot(alpha_fossil_SP_TSQ[,5],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="orange",
     lty=2,
     type="l") #Chao2
polygon(c(alpha_fossil_SP_TSQ[,1],rev(alpha_fossil_SP_TSQ[,1])),
        c(alpha_fossil_SP_TSQ[,6],rev(alpha_fossil_SP_TSQ[,7])),border=NA,col=rgb(1, 1, 0,0.5))#orange

plot(alpha_fossil_PB_TSQ[,8],
     ylim=c(10,100),
     xlab="ICE",
     yaxt="n",
     col="blue",
     lty=2,
     type="l") #ICE
polygon(c(alpha_fossil_PB_TSQ[,1],rev(alpha_fossil_PB_TSQ[,1])),
        c(alpha_fossil_PB_TSQ[,9],rev(alpha_fossil_PB_TSQ[,10])),border=NA,col=rgb(0, 0, 1,0.5))#red

par(new=TRUE)
plot(alpha_fossil_PM_TSQ[,8],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="red",
     lty=2,
     type="l") #ICE
polygon(c(alpha_fossil_PM_TSQ[,1],rev(alpha_fossil_PM_TSQ[,1])),
        c(alpha_fossil_PM_TSQ[,9],rev(alpha_fossil_PM_TSQ[,10])),border=NA,col=rgb(1, 0, 0,0.5))

par(new=TRUE)
plot(alpha_fossil_SP_TSQ[,8],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="orange",
     lty=2,
     type="l") #ICE
polygon(c(alpha_fossil_SP_TSQ[,1],rev(alpha_fossil_SP_TSQ[,1])),
        c(alpha_fossil_SP_TSQ[,9],rev(alpha_fossil_SP_TSQ[,10])),border=NA,col=rgb(1, 1, 0,0.5))

plot(alpha_fossil_PB_TSQ[,11],
     ylim=c(10,100),
     col="blue",
     yaxt="n",
     lty=2,
     type="l") #Jacknife
polygon(c(alpha_fossil_PB_TSQ[,1],rev(alpha_fossil_PB_TSQ[,1])),
        c(alpha_fossil_PB_TSQ[,12],rev(alpha_fossil_PB_TSQ[,13])),border=NA,col=rgb(0, 0, 1,0.5))

par(new=TRUE)  
plot(alpha_fossil_PM_TSQ[,11],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="red",
     lty=2,
     type="l") #Jacknife
polygon(c(alpha_fossil_PM_TSQ[,1],rev(alpha_fossil_PM_TSQ[,1])),
        c(alpha_fossil_PM_TSQ[,12],rev(alpha_fossil_PM_TSQ[,13])),border=NA,col=rgb(1, 0, 0,0.5))

par(new=TRUE)
plot(alpha_fossil_SP_TSQ[,11],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     lty=2,
     col="orange",
     type="l") #Jacknife
polygon(c(alpha_fossil_SP_TSQ[,1],rev(alpha_fossil_SP_TSQ[,1])),
        c(alpha_fossil_SP_TSQ[,12],rev(alpha_fossil_SP_TSQ[,13])),border=NA,col=rgb(1, 1, 0,0.5))

legend(3,20, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       text.col = "black",
       lty=c(2, 2, 2),
       bty="n",cex=1)


########### beta diversity ###################
MMB_merge_count_soilcore_rafy<-decostand(MMB_merge_count_soilcore[,2:49],"pa",2)

NumToChr <- function(x) {if(is.numeric(x)) 
{as.character(paste(x))} else {x}} 

Soilcore<- NumToChr(Soilcore)

betad_rafy<-betadiver(t(MMB_merge_count_soilcore_rafy),"z")

beta.dist<-vegdist(MMB_merge_count_soilcore_rafy,"binomial",binary=TRUE) # define the dissimilarity indexes
adonis(betad_rafy~Extract, method="Jacaard", strata=Library,perm=200)
adonis(betad_rafy~Library, method="Jacaard", strata=Extract,perm=200)
#The function also finds indices for presence/ absence data by setting binary = TRUE

Extract_Lib <- paste(Extract,Library)
par(mfrow=c(1,1))
par(mar = c(4, 4, 2, 0), oma = c(1, 1, 1, 0.5))
mod_rafy<-betadisper(betad_rafy,Extract_Lib)
plot(mod_rafy,
     col= c("blue","blue","red","red","orange","orange"),
     main="Multivariate Dispersion",
     pch = c(15,22,17,24,18,23),
     segments = FALSE, # control grey line
     hull = FALSE, ellipse = FALSE,label = F,lwd=0.5)
ordiellipse(mod_rafy, Extract_Lib, 
            kind="se", 
            conf=0.95, 
            lwd=1, 
            lty = c(1,2,1,2,1,2),
            col= c("blue","blue","red","red","orange","orange")
)

legend(-0.5,0.3, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       text.col = "black", 
       title="Nextera",
       pch = c(15, 17, 18),
       lty=0,
       bty="n",lwd=1,cex=0.75)
legend(-0.5,0.19, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       title="TruSeq",
       text.col = "black", 
       pch = c(22, 24, 23),
       lty=0,
       bty="n",lwd=1,cex=0.75)




