#This script is directly executed after GSE17048 script

#Cleaning up workspace (limited because some variables of the earlier mentioned script will be used)
rm(dat_b, dat_dat, dat, dat_N, desc_edited, desc, desc_HC)
#Renaming the variables
dat   <- dat_L
desc  <- desc_edited_b
#Removing the old names
rm(desc_edited_b, dat_L)

#Setting workplace
setwd("C:/Users/Patrick/Documents/QC/raw/Limma Thesis/Limma GSE17048")


### Limma part ###

#Creating a back up for the complete data ste
datasetcomplete <- dat
#Selecting the only the data of the table
dat <- dat[,2:145]
#Creating a dat back up
dat_b <- dat

#The row names of dat are not probes so let change that
row.names(dat) <- datasetcomplete[,1]

#create a filtered object of expressed probes
#by selecting for an average expression over all samples of at least 3
datF <- dat[rowMeans(dat)>=3,]




#########################
##statistical modelling##
#########################

#load limma library
library(limma)

#model design: take group_ as a fixed effect for an intercept model
design <- model.matrix(~Characteristics..ms.subtype.,data=desc)
colnames(design) <- gsub("[()]","",colnames(design))


#fit model
fit <- lmFit(datF,design)
#Compute LogFC and P-values
fit <- eBayes(fit)



#build the (contrast) matrix to compute some group differences of interest based on the model parameters
cont.matrix <- makeContrasts(
  PPMS = Characteristics..ms.subtype.PPMS,
  RRMS = Characteristics..ms.subtype.RRMS,
  SPMS = Characteristics..ms.subtype.SPMS,
  #SPMSvsRRMS = Characteristics..ms.subtype.SPMS - Characteristics..ms.subtype.RRMS , #contrast between RRMS en SPMS excluded of this research
  levels = colnames(design)
)
#compute the contrast fits
contrast.fit <- contrasts.fit(fit, cont.matrix)
contrast.fit <- eBayes(contrast.fit)

#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
files.c <- saveStatOutput(cont.matrix,contrast.fit,postfix="filt",annotation=desc)

#create summary table of the contrast results
createPvalTab(files.c,postfix="filt_c",html=TRUE)


