#clear work space
rm(list=ls())


#Loading E-MTAB-4890 Limma result
#Reading PPMS
PPMS_FC <- read.delim("C:/Users/Patrick/Documents/QC/raw/Limma Thesis/Limma 4890/table_PPMS_filt.tab")
#Reading the RRMS
RRMS_FC <- read.delim("C:/Users/Patrick/Documents/QC/raw/Limma Thesis/Limma 4890/table_RRMS_filt.tab")
#Reading the SPMS
SPMS_FC <- read.delim("C:/Users/Patrick/Documents/QC/raw/Limma Thesis/Limma 4890/table_SPMS_filt.tab")

#Creating an empty data set
PPMS_interesting <- PPMS_FC[1,]
#Creating an indicator value
V <- 0
#Adding data points that have are significant 
for(Probe in PPMS_FC$X) {
  V <- V + 1
  if(PPMS_FC$adj.P.Val[V] <0.05 & abs(PPMS_FC$Fold.Change[V]) > 1.3) {
    PPMS_interesting[nrow(PPMS_interesting)+1,] = PPMS_FC[V,]
  }
}
#Removing the placeholder 
PPMS_interesting <- PPMS_interesting[2:nrow(PPMS_interesting),1:7]

#Creating an empty data set
RRMS_interesting <- RRMS_FC[1,]
#Creating an indicator value
V <- 0
#Adding data points that are significant 
for(Probe in RRMS_FC$X) {
  V <- V + 1
  if(RRMS_FC$adj.P.Val[V] <0.05 & abs(RRMS_FC$Fold.Change[V]) > 1) {
    RRMS_interesting[nrow(RRMS_interesting)+1,] = RRMS_FC[V,]
  }
}
#Removing the placeholder 
RRMS_interesting <- RRMS_interesting[2:nrow(RRMS_interesting),1:7]

#Creating an empty data set
SPMS_interesting <- SPMS_FC[1,]
#Creating an indicator value
V <- 0
#Adding data points that have are significant 
for(Probe in SPMS_FC$X) {
  V <- V + 1
  if(SPMS_FC$adj.P.Val[V] <0.05 & abs(SPMS_FC$Fold.Change[V]) > 1) {
    SPMS_interesting[nrow(SPMS_interesting)+1,] = SPMS_FC[V,]
  }
}
#Removing the placeholder 
SPMS_interesting <- SPMS_interesting[2:nrow(SPMS_interesting),1:7]

#Reading the Probes and its labels (THIS IS ILLUMNA CHIP)
Probes <- read.delim("C:/Users/Patrick/Documents/DataSets/E-MTAB-4890/A-MEXP-931.adf_edited.txt")
#Matching the probes
Probes_PPMS <- Probes[match(PPMS_interesting$X, Probes$Reporter.Name),c(1,3,4)]
#Checking
data.frame(table(PPMS_interesting$X == Probes_PPMS$Reporter.Name)) #True

#Matching the probes
Probes_RRMS <- Probes[match(RRMS_interesting$X, Probes$Reporter.Name),c(1,3,4)]
#Checking
data.frame(table(RRMS_interesting$X == Probes_RRMS$Reporter.Name)) #True

#Matching the probes
Probes_SPMS <- Probes[match(SPMS_interesting$X, Probes$Reporter.Name),c(1,3,4)]
#Checking
data.frame(table(SPMS_interesting$X == Probes_SPMS$Reporter.Name)) #True


#Setting the Probes_X appart
PPMS_4890 <- Probes_PPMS
RRMS_4890 <- Probes_RRMS
SPMS_4890 <- Probes_SPMS

#Loading GEO17048
#Reading PPMS
PPMS_FC <- read.delim("C:/Users/Patrick/Documents/QC/raw/Limma Thesis/Limma GSE17048/table_PPMS_filt.tab")
#Reading the RRMS
RRMS_FC <- read.delim("C:/Users/Patrick/Documents/QC/raw/Limma Thesis/Limma GSE17048/table_RRMS_filt.tab")
#Reading the SPMS
SPMS_FC <- read.delim("C:/Users/Patrick/Documents/QC/raw/Limma Thesis/Limma GSE17048/table_SPMS_filt.tab")

#Creating an empty data set
PPMS_interesting <- PPMS_FC[1,]
#Creating an indicator value
V <- 0
#Adding data points that have are significant 
for(Probe in PPMS_FC$X) {
  V <- V + 1
  if(PPMS_FC$P.Value[V] <0.05 & abs(PPMS_FC$Fold.Change[V]) > 1.3 ) { #Removing & abs(PPMS_FC$Fold.Change[V]) > 
    PPMS_interesting[nrow(PPMS_interesting)+1,] = PPMS_FC[V,]
  }
}
#Removing the placeholder 
PPMS_interesting <- PPMS_interesting[2:nrow(PPMS_interesting),1:7]

#Creating an empty data set
RRMS_interesting <- RRMS_FC[1,]
#Creating an indicator value
V <- 0
#Adding data points that are significant 
for(Probe in RRMS_FC$X) {
  V <- V + 1
  if(RRMS_FC$P.Value[V] <0.05 ) {
    RRMS_interesting[nrow(RRMS_interesting)+1,] = RRMS_FC[V,]
  }
}
#Removing the placeholder 
RRMS_interesting <- RRMS_interesting[2:nrow(RRMS_interesting),1:7]

#Creating an empty data set
SPMS_interesting <- SPMS_FC[1,]
#Creating an indicator value
V <- 0
#Adding data points that have are significant 
for(Probe in SPMS_FC$X) {
  V <- V + 1
  if(SPMS_FC$P.Value[V] <0.05 & abs(PPMS_FC$Fold.Change[V]) > 1.3) {
    SPMS_interesting[nrow(SPMS_interesting)+1,] = SPMS_FC[V,]
  }
}
#Removing the placeholder 
SPMS_interesting <- SPMS_interesting[2:nrow(SPMS_interesting),1:7]

#Reading probe labels
Probe_info <- read.delim("C:/Users/Patrick/Documents/DataSets/GSE17048/A-MEXP-1171.adf_edited.txt")

#Matching the probes
Probes_PPMS <- Probe_info[match(PPMS_interesting$X, Probe_info$Reporter.Name),c(1,5,4)] #Few 2 no refseq or hugo
#Checking
data.frame(table(PPMS_interesting$X == Probes_PPMS$Reporter.Name)) #True

#Matching the probes
Probes_RRMS <- Probe_info[match(RRMS_interesting$X, Probe_info$Reporter.Name),c(1,5,4)]
#Checking
data.frame(table(RRMS_interesting$X == Probes_RRMS$Reporter.Name)) #True

#Matching the probes
Probes_SPMS <- Probe_info[match(SPMS_interesting$X, Probe_info$Reporter.Name),c(1,5,4)]
#Checking
data.frame(table(SPMS_interesting$X == Probes_SPMS$Reporter.Name)) #True

#Setting the Probes_X appart
PPMS_17048 <- Probes_PPMS
RRMS_17048 <- Probes_RRMS
SPMS_17048 <- Probes_SPMS



#Cleaning workspace
rm(Probes,Probe_info, PPMS_interesting,PPMS_FC, Probes_PPMS,Probes_RRMS,Probes_SPMS, RRMS_FC, RRMS_interesting, SPMS_FC, SPMS_interesting)

#Determing overlapp
PPMS_17048vs4890 <- PPMS_4890[match(PPMS_17048$Reporter.Database.Entry.hugo., PPMS_4890$Reporter.Database.Entry.hugo.),]
#Removing NAs
PPMS_17048vs4890 <- PPMS_17048vs4890[!is.na(PPMS_17048vs4890$Reporter.Database.Entry.hugo.),]

#The values were checked with manual edditing the X to 10 every one was true
X <- 10
data.frame(table(PPMS_17048vs4890$Reporter.Database.Entry.hugo.[X] == PPMS_17048$Reporter.Database.Entry.hugo. )) #TRue
data.frame(table(PPMS_17048vs4890$Reporter.Database.Entry.hugo.[X] == PPMS_4890$Reporter.Database.Entry.hugo. )) #True


#Determing overlapp
RRMS_4890vs17048 <- RRMS_17048[match(RRMS_4890$Reporter.Database.Entry.refseq., RRMS_17048$Reporter.Database.Entry.refseq.), 2:3]
RRMS_4890vs17048 <- RRMS_4890[match(RRMS_4890$Reporter.Database.Entry.refseq., RRMS_17048$Reporter.Database.Entry.refseq.), 2:3] #No mathces both

#Removing NAs
RRMS_4890vs17048 <- RRMS_4890vs17048[!is.na(RRMS_4890vs17048$Reporter.Database.Entry.hugo.),]

#Determing overlapp
SPMS_4890vs17048 <- SPMS_17048[match(SPMS_4890$Reporter.Database.Entry.hugo., SPMS_17048$Reporter.Database.Entry.hugo.), 2:3]

#Removing NAs
SPMS_4890vs17048 <- SPMS_4890vs17048[!is.na(SPMS_4890vs17048$Reporter.Database.Entry.hugo.),]

#Determing values Manual
X <- 8
  data.frame(table(SPMS_4890vs17048$Reporter.Database.Entry.hugo.[X]  == SPMS_17048$Reporter.Database.Entry.hugo. ))
  data.frame(table(SPMS_4890vs17048$Reporter.Database.Entry.hugo.[X] ==  SPMS_4890$Reporter.Database.Entry.hugo. ))

#Selecting only Hugo
  PPMS <- PPMS_17048vs4890$Reporter.Database.Entry.hugo.
  SPMS <- SPMS_4890vs17048$Reporter.Database.Entry.hugo.

#Saving the gene names
  #Saving the dataset
  write.table(PPMS, file="C:/Users/Patrick/Documents/PPMS_Gene_17048vs4890.txt", sep="\t",quote=FALSE,col.names=NA)  
  write.table(SPMS, file ="C:/Users/Patrick/Documents/SPMS_Gene_17048vs4890.txt", sep="\t",quote=FALSE,col.names=NA)

  
#Important#
  #The PPMS SPMS were saved and puto into excel to use the HUGO database itself to determine the ensemble lables, these labels will be used to chearch the 5151 database when the selection happend
  #This is done manual
  
#Reading E-MTAB-5151
#Read the made PPMS file
PPMS_FC <- read.delim("C:/Users/Patrick/Documents/QC/raw/Limma Thesis/Limma 5151/table_PPMS_filt.tab")
#Read the made RRMS file
RRMS_FC <- read.delim("C:/Users/Patrick/Documents/QC/raw/Limma Thesis/Limma 5151/table_RRMS_filt.tab")
#Reading the made SPMS file
SPMS_FC <- read.delim("C:/Users/Patrick/Documents/QC/raw/Limma Thesis/Limma 5151/table_SPMS_filt.tab")

#Creating an empty data set
PPMS_interesting <- PPMS_FC[1,]
#Creating an indicator value
V <- 0
#Adding data points that have are significant 
for(Probe in PPMS_FC$X) {
  V <- V + 1
  if(PPMS_FC$adj.P.Val[V] <0.05 & abs(PPMS_FC$Fold.Change[V]) > 1.3) {
    PPMS_interesting[nrow(PPMS_interesting)+1,] = PPMS_FC[V,]
  }
}
#Removing the placeholder 
PPMS_interesting <- PPMS_interesting[2:nrow(PPMS_interesting),1:7]

#Creating an empty data set
RRMS_interesting <- RRMS_FC[1,]
#Creating an indicator value
V <- 0
#Adding data points that are significant 
for(Probe in RRMS_FC$X) {
  V <- V + 1
  if(RRMS_FC$adj.P.Val[V] <0.05 & abs(RRMS_FC$Fold.Change[V]) > 1.3) {
    RRMS_interesting[nrow(RRMS_interesting)+1,] = RRMS_FC[V,]
  }
}
#Removing the placeholder 
RRMS_interesting <- RRMS_interesting[2:nrow(RRMS_interesting),1:7]

#Creating an empty data set
SPMS_interesting <- SPMS_FC[1,]
#Creating an indicator value
V <- 0
#Adding data points that have are significant 
for(Probe in SPMS_FC$X) {
  V <- V + 1
  if(SPMS_FC$adj.P.Val[V] <0.05 & abs(SPMS_FC$Fold.Change[V]) > 1.3) {
    SPMS_interesting[nrow(SPMS_interesting)+1,] = SPMS_FC[V,]
  }
}
#Removing the placeholder 
SPMS_interesting <- SPMS_interesting[2:nrow(SPMS_interesting),1:7]

#Reading data set with labels to find the genes that match
dat_5151 <- read.delim("C:/Users/Patrick/Documents/DataSets/E-MTAB-5151/E-MTAB-5151-complete with labels.txt")

#Matching labels is done later
#Saving the interest
PPMS_5151 <- PPMS_interesting
RRMS_5151 <- RRMS_interesting
SPMS_5151 <- SPMS_interesting


#Laebling if it is present in 5151
PPMS_5151_Labels <- dat_5151[match(PPMS_5151$X, dat_5151$probe_id), 2:3] #Manual checked if the ensemble was present in the coresponding table
RRMS_5151_Labels <- dat_5151[match(RRMS_5151$X, dat_5151$probe_id), 2:3] #RRMS had no overlap
SPMS_5151_Labels <- dat_5151[match(SPMS_5151$X, dat_5151$probe_id), 2:3]

