#E-GEOD-17048 script
#clear workspace
rm(list=ls())

#set the working directory
setwd("C:/Users/Patrick/Documents/DataSets/GSE17048")

#load needed libraries
library(affy)
library(limma)
library(bioDist)
library(gplots)
library(gcrma)

##include some functions from ArrayAnalysis.org
source ("C:/Users/Patrick/Documents/Scripts/Given scripts/functions_ArrayAnalysis_v2.R")


#Reading the data file
dat <- read.delim("C:/Users/Patrick/Documents/DataSets/GSE17048/Processed Data/E-GEOD-17048-processed-data-2460642653.txt")
#First row contains "Unknown:VALUE" so the row is removed
dat_b <- dat
dat <- dat[-1,]

#Reading Meta Data
desc <- read.delim("C:/Users/Patrick/Documents/DataSets/GSE17048/Processed Data/E-GEOD-17048.sdrf_edited.txt")
desc_edited <- cbind(desc["Source.Name"],desc["Characteristics..ms.subtype."],desc["Characteristics..tissue."])

#Creating a back up
dat_b <- dat
#Creating a back up
desc_edited_b <- desc_edited

#The data names does not include to what subtype it belongs to
#Making a list of the colnames from data
dat_colnames <- colnames(dat_b[2:145])

#The healthy controles is not mentioned in the meta data
##With not mentioned is that they are blank so first a check
###A later column of meta data contains "blood, HC" in desc$Comment..Sample_source_name.

#Selecting the HC ones
desc_HC <- desc[desc$Comment..Sample_source_name.[] == "blood, HC",] #45 samples
#Checking if the HC ones are the same in desc_edited
data.frame(table(desc_edited$Characteristics..ms.subtype. == "  ")) #45 true so it is the case

#Making a list of sample names from desc_edited_b
desc_names <- desc_edited_b$Source.Name
#Making a list of the string length of each sample name
desc_Lenght <- nchar(desc_names)
#Creating DN value that will be used to select the right string
DN <- 0
#Creating Desc_names_edited that will be used to place the new strings to it
Desc_names_edited <- "Placeholder"
#Creating DE that will be used to extend the list
DE <- 0

#for loop that shorten the string #this for loop removes the spacebar after the sample name
for(DL in desc_Lenght) {
  DN <- DN + 1
  DE <- DE + 1
  if(DL == 10) {
    #select until the 9th character of the string and add it to Desc_names_edited
    Desc_names_edited[DE] <- substring(desc_names[DN], first = 0, last = 9)
  } 
}
#Control
nchar(Desc_names_edited)
#Replacing the old sample names with the new ones
desc_edited_b$Source.Name <- Desc_names_edited
desc_edited$Source.Name <- Desc_names_edited


#Ordering it that it match
desc_edited_b <- desc_edited_b[match(dat_colnames, desc_edited_b$Source.Name),]
#Checking if it went right
data.frame(table(desc_edited_b$Source.Name == dat_colnames)) #True


#Creating a value to index a select a specific data
I <- 0
#A for loop to replace the empty space in sybtype column with unkown and rewrite the Source name if it is known
for (x in desc_edited_b$Characteristics..ms.subtype.) {
  I <- I + 1
  if(x == "  ") {
    desc_edited_b$Characteristics..ms.subtype.[I] <- "HC"
  } 
}

#The data in the data frame is characters this makes it not possible to run certain functions
##e.g. range etc so the data is converted to numeric

#selecting only the data
dat_dat <- dat[,2:145]

#setting the data points to numeric
dat_dat[,1] <- as.numeric(dat_dat[,1])
#Above workes so now in for loop
for(x in 2:144){
  dat_dat[,x] <- as.numeric(dat_dat[,x])
}

#Range
range(dat_dat) #-32.73534 31607.63000

#the data will be increased with 33, because log does not like negative values
dat_N <- dat_dat + 33
#Checking if the addition is true 
data.frame(table(dat_N ==  (dat_dat + 33))) #Yes

#putting the dataset in log scale 
dat_L <- cbind(dat_b[,1],log2(dat_N[,1:144]))                                                     


#Creating factors
desc_edited_b$Characteristics..ms.subtype. <- factor(desc_edited_b$Characteristics..ms.subtype., levels = c("HC","PPMS", "RRMS", "SPMS" ))

                                                                                                             
#create QC plots for raw data
factors <- c("Characteristics..ms.subtype.")
setwd("C:/Users/Patrick/Documents/QC/raw")
createQCPlots(dat_L[,2:145], factors, Table=desc_edited_b, normMeth="", postfix="")

#Saving the dataset
write.table(dat_L, file="C:/Users/Patrick/Documents/DataSets/GSE17048/Processed Data/E-GEOD-17048-dat_L_V2.txt", sep="\t",quote=FALSE,col.names=NA)  
write.table(desc_edited_b, file ="C:/Users/Patrick/Documents/DataSets/GSE17048/Processed Data/MetaData_17048_V2.txt", sep="\t" )
