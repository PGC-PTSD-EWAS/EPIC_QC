### Use with location.based.pc() function

### This code can be used as an example or can be directly modified by the user to compute "location-based"
### principal components as described in Barfield et al. (2014) Gen Epid, in press.

### file.pathname should be a character string containing the local pathname for the annotation file: cpgs_within_0bp_of_TGP_SNP.Rdata

### Modify this line to indicate correct pathname, if the annotation file that we've shared is in the working dir, the code will be:
file.pathname = "cpgs_within_0bp_of_TGP_SNP.Rdata"

### beta.obj is the users data; should be a data.frame or matrix of beta values or M-values with:
###### one row per CpG site, one column per sample
###### row names = CpG names, column names = sample names

### Modify this line to incorporate user methylation data
## Load beta values (not ComBAT adjusted beta values, just normalized beta values)
library(data.table)
beta.obj <- fread("noob_qcd_crossReactiveProbesRemoved.csv", data.table = F)
rownames(beta.obj)<-beta.obj$V1
beta.obj<-beta.obj[,-1]

### Load and call location.based.pc() function that we've shared to compute principal components
### This function will select only CpG sites in the location-based annotation file,
### set missing values to the CpG-site-average, and compute principal components

source("location.based.pc.R")
pc <- location.based.pc(beta.obj,file.pathname)

### The returned result will be a princomp object.  The principal components will be available as pc$loadings

### Principal components can be incorporated as covariates in regression analysis
### Samples will be sorted in the same order as beta.obj, but if using with other data, make sure IDs match up!

top10pc <- pc$loadings[,1:10]

## You can merge the mPCs with phenotype file
# Load phenotype file: Samples as rows, the first column should be methylation IDs (same as column names of beta matrix)
phen <- read.csv("Pheno_QC.csv",stringsAsFactors=FALSE,header=TRUE)
phen <- merge(phen,top10pc,by.x = 1,by.y = "row.names",all.x = T)
write.csv(phen, file = "Pheno_QC.csv",row.names = F) # Save the phenotype file with mPCs


