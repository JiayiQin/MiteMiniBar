#set the work routine
setwd("/Users/JiayiQin/Dropbox/Paper_III_MMB/FINAL/")
#Import the dataset

MMB_combine <- read.csv("MMB.csv",header = T, sep =",",check.names=FALSE)

MMB_combine <- MMB_combine[order(MMB_combine$Taxa_result_CGG),]

TREAT<-read.csv("Treat.csv",header=TRUE,check.names=FALSE)

###### Organise taxa annotation  -- ROW ----------------------------------------
# MMB_combine is the table for further analysis
# 52 variables for RISO
# 156 variables for CGG

# Dim function to convert data from factor to character, *use paste(X)
FacToCha <- function(x) {if(is.factor(x)) 
{as.character(paste(x))} else {x}} 

MMB_combine$Taxa_result_CGG <- FacToCha(MMB_combine$Taxa_result_CGG)
MMB_combine$Taxa_result_riso <- FacToCha(MMB_combine$Taxa_result_riso)

# Manually check if the same amplicon got same taxonomy annotation in Excel
# Based on Taxa_result_CGG, copy taxa result 
TAXA <- MMB_combine$Taxa_result_CGG
for ( i in 1:672) {
  if (MMB_combine$Taxa_result_CGG[i] != "NA") 
  {TAXA[i] = MMB_combine$Taxa_result_CGG[i]}
  else
  {TAXA[i] = MMB_combine$Taxa_result_riso[i]}
}

MMB_taxa = cbind(TAXA,MMB_combine)


# change the NA to value 0 for further analysis
MMB_taxa[is.na(MMB_taxa)]<- 0
# This command will lead to several warning message, 
# because there are several NAs as factor. No problem going further.
# This command will also change integer into numeric automatically.

# Minimize the table, remains the TAXA and OTU reads table
MMB_taxa_clean<- MMB_taxa[,-2:-7]
MMB_taxa_clean<- MMB_taxa_clean[,-54:-57] # extra column as character at the end of the table, need to delete

# merge notablis
MMB_taxa_clean$TAXA<- gsub("k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis_L1",
                                "k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis",
                                MMB_taxa_clean$TAXA)
MMB_taxa_clean$TAXA<- gsub("k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis_L2",
                                "k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis",
                                MMB_taxa_clean$TAXA)
MMB_taxa_clean$TAXA<- gsub("k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis_L4",
                                "k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Parisotoma,s:Parisotoma_notabilis",
                                MMB_taxa_clean$TAXA)

# Remove ? mark
MMB_taxa_clean$TAXA<- gsub("(\\?)",
                           "",
                           MMB_taxa_clean$TAXA)

#Merge OTU with the same taxa
MMB_merge_count<-aggregate(MMB_taxa_clean[,2:209], by=list(Taxa=MMB_taxa_clean$TAXA), FUN=sum)

# remove OTU total reads lower than 5
MMB_merge_count <- MMB_merge_count[rowSums(MMB_merge_count[,2:209])>5,]

# remove collembolan_unassigned
MMB_merge_count <- MMB_merge_count[-which(MMB_merge_count$Taxa=="k:Animilia,p:Arthropoda,c:Collembola"),]
# remove mite_unassigned
MMB_merge_count <- MMB_merge_count[-which(MMB_merge_count$Taxa=="k:Animilia,p:Arthropoda,c:Arachnida"),]
# remove F.candida
MMB_merge_count <- MMB_merge_count[-which(MMB_merge_count$Taxa=="k:Animilia,p:Arthropoda,c:Collembola,o:Entomobryomorpha,f:Isotomidae,g:Folsomia,s:Folsomia_candida"),]

write.csv(MMB_merge_count,"MMB_check.csv") # output the taxa table to check if they are reasonable

# colomn names
row.names(MMB_merge_count)<-MMB_merge_count$Taxa
MMB_merge_count <- MMB_merge_count[,-1]

Alpha<-spp.est(MMB_merge_count)
row.names(Alpha)<-colnames(MMB_merge_count)













tMMB_merge_count <- t(MMB_merge_count)
colnames(tMMB_merge_count) <- tMMB_merge_count [1,]
tMMB_merge_count <- tMMB_merge_count [-1,]
tMMB_merge_count<-as.data.frame(tMMB_merge_count)

#row names
tMMB_merge_count$Index<-row.names(tMMB_merge_count)

# make complete table with OTUs and metadata
TREAT$Index <-as.character(paste(TREAT$Index ))
Comeplete_MMB <- merge(tMMB_merge_count,TREAT,by="Index",all.y=TRUE,all.x = TRUE)

######### compare internally ##########

#### Chao2/ICE/Jackkinfe #####

library(fossil)

#occurrence example with sample data set
data(fdata.mat)
spp.est(fdata.mat, abund = FALSE, counter = FALSE)

