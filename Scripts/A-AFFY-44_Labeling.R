#clear workspace
rm(list=ls())

#set the working directory
setwd("C:/Users/Patrick Mertens/Documents")

#Installing the BiocManager library to get the affymatrix chip
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")  #Affymetrix HG-U133_Plus_2 Array annotation data (chip hgu133plus2)
#It is the chipped use for E-MTAB-69

#load needed libraries
library(affy)
library(limma)
library(bioDist)
library(gplots)
library(gcrma)
library("hgu133plus2.db")

#Opening the other probe information file
A_AFFY_44 <- read.delim("C:/Users/Patrick Mertens/Documents/AFFY probe labes Edited Genebank Refeseq embel for R studio_2.txt", as.is = TRUE)


x <- hgu133plus2ENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes][1:41580])
if(length(xx) > 0) {
  # Get the Ensembl gene IDs for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}

#Changing the Library given information to data.frame
Probes_Ensembl <- as.data.frame(x)



#Making a data set of the probes used in the data set (E-MTAB-69)
dat_69 <- read.delim("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-69/Normalized/E-MTAB-69_norm_gcrma_3_5.txt",as.is=TRUE)
dat_Probes <- dat_69$X
dat_Probes <- data.frame(dat_Probes)

#Making a backup of AFFY
A_AFFY_44_b <- A_AFFY_44
#ordering the probes of AFFY to probes
A_AFFY_44_b <- A_AFFY_44_b[match(dat_Probes$dat_Probes,A_AFFY_44_b$Comment.AECompositeName.),]
row.names(A_AFFY_44_b) <- row.names(dat_Probes)
#Checking if ordering went right
data.frame(table(A_AFFY_44_b$Comment.AECompositeName. == dat_Probes$dat_Probes))#TRUE

#value to select specific genbank value (Is needed by the for loop that replace NA values in genbank column)
V <- 0
#A for loop that replace NA values with an ""
for (na in A_AFFY_44_b$Composite.Element.Database.Entry.genbank.) {
  V <- V + 1
  if(is.na(A_AFFY_44_b$Composite.Element.Database.Entry.genbank.[V])){
    A_AFFY_44_b$Composite.Element.Database.Entry.genbank.[V] <- ""
  }
}


#value to select the row of Probes_Ensembl
T <- 1
#value to select the row of A_AFFY_44_b
J <- 1
#Making an same column sized data frame, because I cannot overwrite/add between in a for loop do I need to append a empty one
Probe_Labels <- Probes_Ensembl[1,1:2]
#Creating a column for Probe_Labels to serve as database indicator
Database <- c("Placeholder")
#Combining it
Probe_Labels <- cbind(Probe_Labels, Database)

#Taking a Probe and giving it a correct label and adding it to Probe_Labels
for (x in dat_Probes$dat_Probes) {
  if(x == Probes_Ensembl$probe_id[T]) {
    Probe_Labels[nrow(Probe_Labels)+1,] = c(Probes_Ensembl[T,1], Probes_Ensembl[T,2], "En")
    T <- T+1
    J <- J+1
  } else { 
      if(A_AFFY_44_b$Composite.Element.Database.Entry.ensembl.[J] != ""){
      #Execution of the if statement
      Probe_Labels[nrow(Probe_Labels)+1,] = c(A_AFFY_44_b[J,1], A_AFFY_44_b[J,2], "En")
      J <- J+1
  } else if(A_AFFY_44_b$Composite.Element.Database.Entry.refseq.[J] != ""){
      #Execution of the else if statement
      Probe_Labels[nrow(Probe_Labels)+1,] = c(x, A_AFFY_44_b$Composite.Element.Database.Entry.refseq.[J], "Q")
      J <- J+1
  } else if(A_AFFY_44_b$Composite.Element.Database.Entry.genbank.[J] != ""){
      #Execution of the else if statement
      Probe_Labels[nrow(Probe_Labels)+1,] = c(x, A_AFFY_44_b$Composite.Element.Database.Entry.genbank.[J], "G")
      J <- J+1
  } else if(A_AFFY_44_b$Composite.Element.Database.Entry.interpro.[J]!= "") {
      #Execution of the else if statement
      Probe_Labels[nrow(Probe_Labels)+1,] = c(x, A_AFFY_44_b$Composite.Element.Database.Entry.interpro.[J], "I")
      J <- J+1
  } else if(A_AFFY_44_b$Composite.Element.Database.Entry.unigene.[J]!= "") {
      #Execution of the else if statement
      Probe_Labels[nrow(Probe_Labels)+1,] = c(x,A_AFFY_44_b$Composite.Element.Database.Entry.unigene.[J], "U")
      J <- J+1
  } else if(A_AFFY_44_b$Composite.Element.Database.Entry.unigene.[J]== "") {
      Probe_Labels[nrow(Probe_Labels)+1,] = c(x,"unknown","unknown")
      J <- J+1
  }
    } #There are still unknowns
}
#The  first row of Probe_Labes is the same as the second row so removing it
Probe_Labels <- Probe_Labels[2:54676, 1:3]
#Creating a back up
Probe_Labels_b <- Probe_Labels

#Creating a indicator value
V <- 0
#Replacing creating a dataset with the same structure as the fist one
Probe_unknown <- Probe_Labels[1,1:3]
#Making one list of all probes with unknown labels
for (I3 in Probe_Labels$ensembl_id) {
  V <- V + 1
  if(I3 == "unknown") {
    Probe_unknown[nrow(Probe_unknown)+1,] = c(Probe_Labels$probe_id[V],I3, "unknown" )
  }
}
#The first row of Probe_unknown was a placeholder lets remove it
Probe_unknown <- Probe_unknown[2:9119,1:3]

#Matching the unknowns to the Probes_Ensemble this pulls also the ensemble ID
Probe_unknown_b <- Probes_Ensembl[match(Probe_unknown$probe_id,Probes_Ensembl$probe_id),]

#value to select specific  value
V <- 0
#A for loop that replace NA values with an ""
for (na in Probe_unknown_b$probe_id) {
  V <- V + 1
  if(is.na(Probe_unknown_b$probe_id[V])){
    Probe_unknown_b[V,] <- ""
  }
}

#Creating an empty dataset with the same structure
Probe_Knowns <- Probe_unknown_b[1,1:2]

#value to select specific  value 
V <- 0
#A for loop that replace NA values with an ""
for (x in Probe_unknown_b$probe_id) {
  V <- V + 1
  if(x != ""){
    Probe_Knowns[nrow(Probe_Knowns)+1,] = c(Probe_unknown_b$probe_id[V], Probe_unknown_b$ensembl_id[V])
  }
}

#The  first row of Probe_Knowns is the same as the second row so removing it
Probe_Knowns <- Probe_Knowns[2:2637, 1:2]
#checking how many unknowns are before
DataFrame(table(Probe_Labels_b$ensembl_id == "unknown")) #9118 are TRUE (unknowns), and 2636 of them are known in Probe_Knowns

#value to select specific  value 
V <- 0
#Value to select specific value
T <- 1

#For loop to replace the unknowns
for(x in Probe_Labels_b$probe_id) {
  V <- V + 1
  if (x == Probe_Knowns$probe_id[T]) {
    Probe_Labels_b$ensembl_id[V] <- Probe_Knowns$ensembl_id[T]
    Probe_Labels_b$Database[V] <- "En"
    T <- T + 1
  }
} #This will give an error eventually because T will grow bigger than Probe_Knowns

#After check
DataFrame(table(Probe_Labels_b$ensembl_id == "unknown")) #6482 are TRUE (unknowns), so it went right

#Cleaning up workspace
rm(A_AFFY_44, A_AFFY_44_b, dat_69, dat_Probes, Probe_Knowns,Probe_Labels,Probe_unknown,Probe_unknown_b, xx, Probes_Ensembl)

#Creating a backup
Probe_Labels <- Probe_Labels_b

#The proccessed data file from E-MTAB-69 also have some IDs there could be ID present that were still unknown
dat <- read.delim("C:/Users/Patrick Mertens/Documents/DataSets/E-MTAB-69/Normalized/E-MTAB-69-A-AFFY-44-normalized-expressions.tsv")

#Creating a value
V <- 0
#Creating an dataset with the same structure
Probe_unknown <- Probe_Labels_b[1,1:3]
#Making one list of all probes with unknown labels
for (I4 in Probe_Labels_b$ensembl_id) {
  V <- V + 1
  if(I4 == "unknown") {
    Probe_unknown[nrow(Probe_unknown)+1,] = c(Probe_Labels_b$probe_id[V],I4, "unknown" )
  }
}

#The first row of Probe_unknown_2 was a place holder lets remove it
Probe_unknown <- Probe_unknown[2:6483,1:3]

#Matching the unknowns to the Probes_Ensemble
Probe_unknown_b <- dat[match(Probe_unknown$probe_id,dat$DesignElementAccession),]

#value to select specific  value 
V <- 0
#A for loop that replace NA values with an ""
for (na in Probe_unknown_b$DesignElementAccession) {
  V <- V + 1
  if(is.na(Probe_unknown_b$DesignElementAccession[V])){
    Probe_unknown_b[V,] <- ""
  }
}

#Creating an dataset with the same structure
Probe_Knowns <- Probe_unknown[1,1:2]
#value to select specific  value 
V <- 0
#A for loop that replace NA values with an ""
for (x in Probe_unknown_b$DesignElementAccession) {
  V <- V + 1
  if(x != ""){
    Probe_Knowns[nrow(Probe_Knowns)+1,] = c(Probe_unknown_b$DesignElementAccession[V], Probe_unknown_b$Gene.ID[V])
  }
}

#The  first row of Probe_Knowns_2 is the same as the second row so removing it
Probe_Knowns <- Probe_Knowns[2:993,1:2]




#checking how many unknowns are before
DataFrame(table(Probe_Labels_b$ensembl_id == "unknown")) #6482 are TRUE (unknown)

#Indicator values
V <- 0
#Value to select specific value
T <- 1

#For loop to replace the unknowns
for(x in Probe_Labels_b$probe_id) {
  V <- V + 1
  if (x == Probe_Knowns$probe_id[T]) {
    Probe_Labels_b$ensembl_id[V] <- Probe_Knowns$ensembl_id[T]
    Probe_Labels_b$Database[V] <- "En"
    T <- T + 1
  }
} ## This for loop will eventually give an error because T get bigger than Probe_Knowns

#After check
DataFrame(table(Probe_Labels_b$ensembl_id == "unknown")) #5490 are TRUE (unknowns), 

#Overwriting the old one
Probe_Labels <- Probe_Labels_b
#Cleaning workspace
rm(dat,Probe_Knowns,Probe_unknown,Probe_unknown_b)

#Loading data from emsemble database itself
Probe_Ensembl <- read.csv("C:/Users/Patrick Mertens/Documents/mart_export (1).txt", as.is = TRUE)

#Creating a value
V <- 0
#Creating an dataset with the same structure
Probe_unknown <- Probe_Labels_b[1,1:3]
#Making one list of all probes with unknown labels
for (I4 in Probe_Labels_b$ensembl_id) {
  V <- V + 1
  if(I4 == "unknown") {
    Probe_unknown[nrow(Probe_unknown)+1,] = c(Probe_Labels_b$probe_id[V],I4, "unknown" )
  }
}
#Removing the Placeholder
Probe_unknown <- Probe_unknown[2:5491,]


#Matching and ordering it against the unknown ones
Probe_Ensembl_b <- Probe_Ensembl[match(Probe_unknown$probe_id,Probe_Ensembl$AFFY.HG.U133.Plus.2.probe),]

#value to select specific  value 
V <- 0
#A for loop that replace NA values with an ""
for (na in Probe_Ensembl_b$AFFY.HG.U133.Plus.2.probe) {
  V <- V + 1
  if(is.na(Probe_Ensembl_b$AFFY.HG.U133.Plus.2.probe[V])){
    Probe_Ensembl_b$AFFY.HG.U133.Plus.2.probe[V] <- ""
  }
}

#Creating an dataset with the same structure
Probe_Knowns <- Probe_unknown[1,1:2]

#value to select specific  value 
V <- 0
#A for loop that replace NA values with an ""
for (x in Probe_Ensembl_b$AFFY.HG.U133.Plus.2.probe) {
  V <- V + 1
  if(x != ""){
    Probe_Knowns[nrow(Probe_Knowns)+1,] = c(Probe_Ensembl_b$AFFY.HG.U133.Plus.2.probe[V], Probe_Ensembl_b$Gene.stable.ID[V])
  }
}

#The  first row of Probe_Knowns is the same as the second row so removing it
Probe_Knowns <- Probe_Knowns[2:227, 1:2]


#The probes are not in the same order as the Probe_labels (or the unknown ones)
Probe_Knowns_order <- Probe_Knowns[match(Probe_unknown$probe_id,Probe_Knowns$probe_id),]

#value to select specific  value 
V <- 0
#A for loop that replace NA values with an ""
for (na in Probe_Knowns_order$probe_id) {
  V <- V + 1
  if(is.na(Probe_Knowns_order$probe_id[V])){
    Probe_Knowns_order[V,] <- ""
  }
}

#Creating an dataset with the same structure
Probe_Knowns_order_cleaned <- Probe_unknown[1,1:2]

#value to select specific  value 
V <- 0
#A for loop that replace NA values with an ""
for (x in Probe_Knowns_order$probe_id) {
  V <- V + 1
  if(x != ""){
    Probe_Knowns_order_cleaned[nrow(Probe_Knowns_order_cleaned)+1,] = c(Probe_Knowns_order$probe_id[V], Probe_Knowns_order$ensembl_id[V])
  }
}
#Removing Placeholder
Probe_Knowns_order_cleaned <- Probe_Knowns_order_cleaned[2:227,]

#Before check
DataFrame(table(Probe_Labels_b$ensembl_id == "unknown")) #5490 are TRUE (unknowns),

#value to select specific  value 
V <- 0
#Value to select specific value
T <- 1

#For loop to replace the unknowns
for(x in Probe_Labels_b$probe_id) {
  V <- V + 1
  if (x == Probe_Knowns_order_cleaned$probe_id[T]) {
    Probe_Labels_b$ensembl_id[V] <- Probe_Knowns_order_cleaned$ensembl_id[T]
    Probe_Labels_b$Database[V] <- "En"
    T <- T + 1
  }
} # This give eventually an error because T gets bigger than the table

#after check
DataFrame(table(Probe_Labels_b$ensembl_id == "unknown")) #5264 are TRUE (unknowns),



#Saving it
write.table(Probe_Labels_b,file="A-AFFY-44_Labeling_Product.txt",sep="\t",quote=FALSE,col.names=NA)
