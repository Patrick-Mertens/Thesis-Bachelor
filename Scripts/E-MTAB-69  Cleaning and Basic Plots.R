#clear workspace
rm(list=ls())

#set the working directory
setwd("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-69")

#load needed libraries
library(affy)
library(limma)
library(bioDist)
library(gplots)
library(gcrma)

#include some functions from ArrayAnalysis.org
source ("C:/Users/Patrick Mertens/Documents/Scripts/functions_ArrayAnalysis_v2.R")

#reading data files
dat <- read.delim("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-69/Normalized/E-MTAB-69-A-AFFY-44-normalized-expressions.tsv")

#Opening this if it is a new session
desc <- read.delim("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-69/E-MTAB-69.sdrf_edited.txt",as.is=TRUE)
#Naming the rows of desc the same as data columns of dat
rownames(desc) <- desc$Source.Name
#Matching the desc rows with the data columns of dat
desc <- desc[match(colnames(dat[4:91]), row.names(desc)),] #First 3 columns of dat are Gene.ID, Gene.Name, DesignElementAccession

#Checking if it went right
data.frame(table(row.names(desc)==colnames(dat[4:91]))) #is True

#Making a table of control values
table(desc$Staging) # NA = 36, relapse = 24, remission = 28

#Is needed for the for loop below
Z <- 0
#Replacing Staging blanks with unknown
for(x in desc$Staging) {
  Z <- Z + 1
  if(x == "  ") {
    desc$Staging[Z] <- "unknown"
  }
}

#Making a table of control values
table(desc$Staging) # unknown = 36, relapse = 24, remission = 28

#Giving the columns factors needed
desc$Staging <- factor(desc$Staging,levels=c("remission","relapse", "unknown"))
desc$Source <- factor(desc$Source,levels=c("CSF","blood"))
desc$Disease <- factor(desc$Disease,levels=c("non.Infl","MS")) #non.infill is not healthy control

#Creating Sub sets based on CSF 
desc_CSF <- desc[desc$Source!="blood",]
desc_CSF_I <- desc_CSF[desc_CSF$Disease!="MS",]
desc_CSF_M <- desc_CSF[desc_CSF$Disease!="non.Infl",] 

#Makes a data set of only the data samples that are CSF
dat_CSF <- cbind(dat[,1:3],dat[,row.names(desc_CSF)])

#Creating Sub set based on blood
desc_blood <- desc[desc$Source!="CSF",]
desc_blood_I <- desc_blood[desc_blood$Disease!="MS",]
desc_blood_M <- desc_blood[desc_blood$Disease!="non.Infl",]
desc_blood_Remission  <- desc_blood[desc_blood$Staging!="relapse",]
desc_blood_Relapse    <- desc_blood[desc_blood$Staging!="remission",]

#Makes a data set of only the data samples that are Blood
dat_blood <- cbind(dat[,1:3], dat[,row.names(desc_blood)])
dat_blood_M <- cbind(dat[,1:3], dat[,row.names(desc_blood_M)])
dat_blood_Remission <- cbind(dat[,1:3], dat[,row.names(desc_blood_Remission)])
desc_blood_Relapse <- cbind(dat[,1:3], dat[,row.names(desc_blood_Relapse)])

#Checking if the order is done correctly
row.names(desc_CSF) == colnames(dat_CSF[4:47])
row.names(desc_blood) == colnames(dat_blood[4:47])

#create QC plots for raw data
factors <- c("Disease", "Staging", "Source")
setwd("C:/Users/Patrick Mertens/Documents/QC/raw")
createQCPlots(dat[,-c(1:3)], factors, Table=desc, normMeth="", postfix="")


#In the plot looks BB61, BB57, BB44 extreme values So they will be removed
#removing BB61 (56), BB57 (35), and BB44 (3) from desc
desc <- desc[-c(56,35,3),]
#Checking if BB61, BB57, BB44 is removed.
for(File in row.names(desc)) {
  if(File == "BB61") {
    print("Error")
  } else if(File == "BB57") {
    print("Error")
  } else if(File == "BB44") {
    print("Error")
  } 
} #this gave no Error so the values are removed


#Removing BB61 (59), BB57 (38), and BB44 (6) from dat
dat <- dat[-c(59,38,6)]
#Checking if BB61, BB57, BB44 is removed
for(File2 in colnames(dat)) {
  if(File2 == "BB61") {
    print("Error")
  } else if(File2 == "BB57") {
    print("Error")
  } else if(File2 == "BB44") {
    print("Error")
  } 
} #this gave no Error so the values are removed

#Checking if it is the same order still
data.frame(table(row.names(desc)==colnames(dat[4:88]))) #is True



#create QC plots for raw data without extreme values now
factors <- c("Disease", "Staging", "Source")
setwd("C:/Users/Patrick Mertens/Documents/QC/raw")
createQCPlots(dat[,-c(1:3)], factors, Table=desc, normMeth="", postfix="")



