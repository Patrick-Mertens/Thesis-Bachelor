#clear work space
rm(list=ls())

#set the working directory
setwd("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-5151")

#load needed libraries
library(affy)
library(limma)
library(bioDist)
library(gplots)
library(gcrma)

##include some functions from ArrayAnalysis.org
source ("C:/Users/Patrick Mertens/Documents/Scripts/functions_ArrayAnalysis_v2.R")




#CTRL_01 Sample file is loaded 
CTRL_01 <- read.delim("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-5151/Data/CTRL_01.txt", as.is = TRUE)
#Selecting only the rows with data into
CTRL_01 <- CTRL_01[CTRL_01$Probe_ID!="",]
rownames(CTRL_01) <- CTRL_01$Probe_ID
#Removing probe_ID
CTRL_01 <- CTRL_01[,-c(1,3), drop=FALSE]
dataset <- CTRL_01



#Making an list of the data files
DataFiles <- list.files(file.path("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-5151/Data"))

#For loop that reads the files 1 for 1. 
for(individualFile in DataFiles[2:76]) {
  Data <- read.delim(paste("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-5151/Data",individualFile, sep = "/"))
  #removing empty rows and columns
  Data <- Data[Data$Probe_ID!="", 1:3]
  rownames(Data) <- Data$Probe_ID
  Data <- Data[,-c(1,3), drop=FALSE]
  dataset <- cbind(dataset,Data)
  }

#### SAVING THE TABLE ####
setwd("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-5151")
write.table(dataset, file="E-MTAB-5151-data+detection.txt", sep="\t",quote=FALSE,col.names=NA)

#### ####

#cleaning work space
rm(CTRL_01, Data, Data2, DataFiles, individualFile)



#Reads meta data
desc <- read.delim("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-5151/E-MTAB-5151.sdrf_2.txt", as.is = TRUE)
desc_edited <- cbind(desc["Characteristics.disease."],desc["Characteristics.organism.part."],desc["Assay.Name"])

#Creating a dataset backup
dataset_b <- dataset
#naming the colnames as the desc_edited$Assay.name, (They are are in the same order)
colnames(dataset_b) <- desc_edited$Assay.Name
#Checking if the names are right
data.frame(table(colnames(dataset_b)== desc_edited$Assay.Name)) #TRUE
#Check if the data was not changed
data.frame(table(dataset_b == dataset)) #TRUE

#This dataset is already in logscale this can be seen by the range
#E-MTAB-69 is done with the same chip and has a lot of labels that can be used for this dataset
dat_69 <- read.delim("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-69/Normalized/E-MTAB-69-A-AFFY-44-normalized-expressions.tsv")

#Ordering dat_69 so it is in the same order as dataset
dat_69_b <- dat_69[match(row.names(dataset_b), dat_69$DesignElementAccession),]

#value to select specific  value 
V <- 0
#A for loop that replace NA values with an ""
for (na in dat_69_b$DesignElementAccession) {
  V <- V + 1
  if(is.na(dat_69_b$DesignElementAccession[V])){
    dat_69_b$DesignElementAccession[V] <- ""
  }
}

Probe_lables_69 <- dat_69_b[2,c(3,1)]
#value to select specific  value 
V <- 0
#A for loop that replace NA values with an ""
for (x in dat_69_b$DesignElementAccession) {
  V <- V + 1
  if(x != ""){
    Probe_lables_69[nrow(Probe_lables_69)+1,] = c(dat_69_b$DesignElementAccession[V], dat_69_b$Gene.ID[V])
  }
}

#The  first row of Probe_Knowns_2 is the same as the second row so removing it
Probe_lables_69 <- Probe_lables_69[2:21057,1:2]

#It is just 21056 labels let get some more
A_AFFY_44 <- read.delim("C:/Users/Patrick Mertens/Documents/Made code documents/A-AFFY-44_Labeling_Product.txt")
#It imports the old rows as column 1 lets remove that
A_AFFY_44 <- A_AFFY_44[,-1]

#Lets check if the probes from A_AFFY_44 is in the same order as dataset_b
data.frame(table(row.names(dataset_b)== A_AFFY_44$probe_id)) #TRUE

#Lets make a dataset that has the probes, labels, database (indicator), and the data
dataset_complete <- A_AFFY_44[1:3] #Remember the probes from AFFY and dataset are in the same order!
#Creating an indicator value
V <- 0
T <-1
#Let be sure that the dat_69 labels are present in dataset
for(probe in dataset_complete$probe_id) {
  V <- V +1
  if(probe == Probe_lables_69$DesignElementAccession[T]) {
    dataset_complete$ensembl_id[V] <- Probe_lables_69$Gene.ID[T]
    dataset_complete$Database[V] <- "En"
    T<- T +1 
  }
} #This will give eventual an Error (missing value), because dataset_complete is longer than Probe_labels_69

#adding the data to dataset_complete
dataset_complete <- cbind(dataset_complete, dataset_b[1:76])


#Save table
setwd("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-5151")
write.table(dataset_complete, file="E-MTAB-5151-complete with labels.txt", sep="\t",quote=FALSE,col.names=NA)  



#Creating a backup
desc_edited_T <- desc_edited
#indicator value
V<- 0
#For loop to replace long disease names with abbreviation
for(x in desc_edited_T$Characteristics.disease.){
  V <- V +1
  if(x == "normal") {
    desc_edited_T$Characteristics.disease.[V] <- "CTRL"
  } else if ( x == "primary progressive multiple sclerosis ") {
    desc_edited_T$Characteristics.disease.[V] <- "PPMS"
  } else if (x == "relapsing-remitting multiple sclerosis ") {
    desc_edited_T$Characteristics.disease.[V] <- "RRMS"
  } else if (x == "secondary progressive multiple sclerosis") {
    desc_edited_T$Characteristics.disease.[V] <- "SPMS"
  }
}

#Factors is needed for plots
desc_edited_T$Characteristics.disease. <- factor(desc_edited_T$Characteristics.disease., levels = c("CTRL","PPMS", "RRMS", "SPMS" ))

#create QC plots for raw data without extreme values now
factors <- c("Characteristics.disease.")
setwd("C:/Users/Patrick Mertens/Documents/QC/raw")
createQCPlots(dataset_complete, factors, Table=desc_edited_T, normMeth="", postfix="")

#Save table
setwd("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-5151")
write.table(desc_edited_T, file="E-MTAB-5151_edited_Metadata.txt", sep="\t",quote=FALSE,col.names=NA)  

