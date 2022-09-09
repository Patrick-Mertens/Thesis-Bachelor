#Selecting the only the data of the table
dat <- dataset_complete[,4:79]
#Creating a dat back up
dat_b <- dat

#create a filtered object of expressed probes
#by selecting for an average expression over all samples of at least 3
datF <- dat[rowMeans(dat)>=3,]






#########################
##statistical modelling##
#########################



#load limma library
library(limma)

#model design: take group_ as a fixed effect for an intercept model
design <- model.matrix(~Characteristics.disease.,data=desc_edited_T)
colnames(design) <- gsub("[()]","",colnames(design))


#fit model
fit <- lmFit(datF,design)
#Compute LogFC and P-values
fit <- eBayes(fit)



#build the (contrast) matrix to compute some group differences of interest based on the model parameters
cont.matrix <- makeContrasts(
  PPMS = Characteristics.disease.PPMS,
  RRMS = Characteristics.disease.RRMS,
  SPMS = Characteristics.disease.SPMS,
  levels = colnames(design)
)
#compute the contrast fits
contrast.fit <- contrasts.fit(fit, cont.matrix)
contrast.fit <- eBayes(contrast.fit)

#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
files.c <- saveStatOutput(cont.matrix,contrast.fit,postfix="filt",annotation=desc_edited_T)

#create summary table of the contrast results
createPvalTab(files.c,postfix="filt_c",html=TRUE)


