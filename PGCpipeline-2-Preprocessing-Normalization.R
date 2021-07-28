#########################################################################
### PGC Pipeline - Script 2 - Preprocessing and Normalization
#########################################################################

# path to source file required to install the packages
source("/Volumes/Seyma/Methyl/CrossSectionalEPIC/install_needed_packages.R") # PATH TO SOURCE FILE, this is an example path, change accordingly

#' Lets install the required packages
bioc_packages <- c("BiocManager", "minfi", "IlluminaHumanMethylationEPICmanifest",
                   "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                   "sva", "limma", "impute", "BiocParallel")

cran_packages <- c( "data.table", "pamr", "feather", "tibble", "CpGassoc")

#' Call functions
install_bioconductor_pkgs(pkgs = bioc_packages)
install_cran_pkgs(pkgs = cran_packages)

#' Run test
check_installed(pkgs = c(bioc_packages, cran_packages))

#' Load all packages, if needed
lapply(c(bioc_packages, cran_packages), require, character.only = TRUE)

########################################################################################################
######################## 1) Preprocessing ##############################################################
########################################################################################################

#' Function to read the idat files
#' Input: path to the idats, Phenotype file
#' Output: RGset, Mset, GMset in a list
read_idat_files <- function(path, pheno_file){
  RGset <- read.metharray.exp(path, recursive = TRUE, verbose = TRUE, force = TRUE)
  RGset <- RGset[, colnames(RGset) %in% pheno$SampleID]

  if(!all(colnames(RGset) %in% pheno$SampleID)){
    stop("samples in idats and sample sheet are not same")
  }

  Mset <- preprocessNoob(RGset, offset = 15, dyeCorr = TRUE, verbose = TRUE, dyeMethod="single")
  GMset <- mapToGenome(Mset, mergeManifest = FALSE)

  return(list(RGset = RGset, Mset =  Mset, GMset = GMset))
}

# Load the new phenotype file with QC info
# Remove failed samples
pheno <- read.csv("Pheno_QC.csv",as.is = T) # WRITE PHENOTYPE FILE HERE
pheno <- subset(pheno, failed == "FALSE")

#' This is the path to the folders having idats
main_dir <- "G:/methylRawData2020/"  # PATH TO IDATS FOLDER HERE, this is an example path, change accordingly  

#' Now call the function
#' you will get RGset, Mset and GMset in a list
output <- read_idat_files(path = main_dir, pheno_file = pheno)

#' Predict the sex
#' This is the function to match sex with predicted sex
#' And remove the sex mismatching samples from phenotype file, RGset, Mset and GMset
#' Input: list of (RGset, Mset and GMset), phenotype file, name of sex column in your phenotype file
#' Sex should be coded as follows: 0 or M for males and 1 or F for females  
#' Output : predicted sex and sex mismatches

check_sex_info <- function(inputdata , pheno_file, pheno_sex_col){
  message("Processing, please wait ...")
  pheno_file[,pheno_sex_col][pheno_file[,pheno_sex_col] == 0] <- "M"
  pheno_file[,pheno_sex_col][pheno_file[,pheno_sex_col] == 1] <- "F"
  predicted_sex <- getSex(object = inputdata$GMset, cutoff = -2)

  rownames(pheno_file) <- pheno_file$SampleID
  pheno <- pheno_file[order(rownames(pheno_file)), ]
  sex <- predicted_sex[order(rownames(predicted_sex)), ]

  if(!all(rownames(pheno) == rownames(sex))){
    stop("Phenotype and predicted sex samples doesn't match")
  }
  RGset <- inputdata$RGset
  Mset <- inputdata$Mset
  GMset <- inputdata$GMset
  mm <- rownames(pheno[pheno[[pheno_sex_col]] != sex$predictedSex, ])
  message("Samples with sex mismatch: ", length(mm))
  if(length(mm)){
    RGset <- RGset[, which(!colnames(RGset) %in% mm)]
    Mset <- Mset[, which(!colnames(Mset) %in% mm)]
    GMset <- GMset[, which(!colnames(GMset) %in% mm)]
    pheno <- pheno[which(!pheno$SampleID %in% mm), ]
  }
  return(list(pheno = pheno, RGset = RGset, Mset = Mset, GMset = GMset, predicted_sex = predicted_sex))
}

#' Now call the function to check sex information
#' pheno_sex_col = the name of the sex column in the phenotype (Sex, Gender etc)
all_info <- check_sex_info(inputdata = output, pheno_file = pheno, pheno_sex_col = "Sex")

#' We will also check if all data have same samples
#' All shoudl be true
lapply(all_info[c(2:4)], function(x) all(colnames(x) %in% all_info$pheno$SampleID))

# Write predicted sex onto the phenotype file and flag the mismatches
pheno <- pheno[order(pheno$SampleID),]
all_info$predicted_sex <- all_info$predicted_sex[order(rownames(all_info$predicted_sex)),]
all(rownames(all_info$predicted_sex)==pheno$SampleID) # This should be TRUE
pheno$predictedSex <- all_info$predicted_sex$predictedSex
if(all(0|1%in%pheno$Sex) == T) {
  pheno$predictedSex[pheno$predictedSex == "M"] <- 0 
  pheno$predictedSex[pheno$predictedSex == "F"] <- 1
}
pheno$SexMismatch <- pheno$Sex != pheno$predictedSex 

# Extract preprocessed pval, beta, methylated and unmethylated signal files
pval <- detectionP(all_info$RGset)
beta <- getBeta(all_info$Mset)
signalA <- getUnmeth(all_info$Mset)
signalB <- getMeth(all_info$Mset)

# Save preprocessed pval, beta, methylated and unmethylated signal files
save(pval, beta, signalA, signalB, file = "IntermediaryFiles.RData")

# Save the new phenotype file with flagged sex mismatches
write.csv(pheno, file = "Pheno_QC.csv", row.names = F)

#' remove all objects from your workspace
rm(list=ls())
gc()

########################################################################################################
######################## 2) QC with CpGassoc ###########################################################
########################################################################################################

# Load the preprocessed pval, beta, methylated and unmethylated signal files

load("IntermediaryFiles.RData")

beta.new <- cpg.qc(beta, signalA, signalB, pval, p.cutoff=.01,cpg.miss=.1, sample.miss=.1)

#[1] "Removed 0 samples with low signal"
#[1] "Removed 1531 CpG sites with missing data for > 0.1 of samples"
#[1] "Removed 0 samples with missing data for > 0.1 of CpG sites"

## Remove cross reactive probes: Got the cross reactive probes from a paper (http://www.sciencedirect.com/science/article/pii/S221359601630071X).
#Processing QCd data to remove only Cross Reactive Probes
cross <- read.csv("cross.csv", stringsAsFactors = FALSE, header = TRUE) # This file is available with the pipeline package we've shared
cross <- cross[, 1]
beta.new <- beta.new[!(row.names(beta.new) %in% cross), ]

#' Run test if all cross cpgs are removed
stopifnot(all(!cross %in% rownames(beta.new)))

# Save Noob Normalized Cross-Reactive Probes removed beta values
write.csv(beta.new, file = "noob_qcd_crossReactiveProbesRemoved.csv", quote = FALSE, row.names = TRUE) # keep the noob normalized and QCed data

########################################################################################################
############# 3) ComBAT normalization to adjust for batch effects (chip and array position) ############
########################################################################################################
#' Clean
rm(list=ls())
gc()

# Load Noob Normalized Cross-Reactive Probes removed beta values
beta <- fread("noob_qcd_crossReactiveProbesRemoved.csv", data.table = F) # cpgs x samples
row.names(beta)<-beta$V1
beta<-beta[,-1]

# Load the new phenotype file with QC info
phen <- read.csv("Pheno_QC.csv", stringsAsFactors = FALSE, header = TRUE)
row.names(phen)<-pheno$SampleID

# Remove the failed samples and sex mismatches
phen <- subset(phen, failed == "FALSE" & SexMismatch == "FALSE")

# Methylation IDs (Row names of phenotype file and Column names of beta matrix) should match
beta <- beta[, which(names(beta) %in% row.names(phen))]
phen <- phen[row.names(phen) %in% names(beta), ]
beta <- beta[, order(names(beta))]
phen <- phen[order(row.names(phen)), ]
stopifnot(all(colnames(beta) == rownames(phen)))

## Define Variables for Combat step
## Define model.matrix, which includes the variables that you want to protect when adjusting for chip and position
## Generally variables that we use as covariates in the EWAS (sex, age, main phenotype -PTSD-, smoking) are included in the model.matrix
sex <- "Sex"         # Name of Sex Variable: Males coded as 0, females coded as 1
age <- "Age"         # Name of Age Variable
PTSD <- "PTSD"       # Name of PTSD Variable: Cases coded as 1, controls coded as 0
chip <- "SentrixID"
position <- "SentrixPosition"

## You should not have NAs in model matrix, so we remove subjects with no phenotype info
print(paste0("Samples with no PTSD information = ", sum(is.na(phen[,PTSD]))))
print(paste0("Samples with no Sex information = ", sum(is.na(phen[,sex]))))
print(paste0("Samples with no Age information = ", sum(is.na(phen[,age]))))

naindex <- (!is.na(phen[,PTSD]) & !is.na(phen[, age]) & !is.na(phen[, sex]))
phen <- phen[naindex, ]
beta <- beta[, naindex]
stopifnot(all(colnames(beta) == rownames(phen)))

chip <- as.factor(phen[,chip])
position <- as.factor(phen[,position])
ptsd<-as.factor(phen[,PTSD])
age<-as.numeric(phen[,age])
sex<-as.factor(phen[,sex])

moddata <- model.matrix(~ptsd+age+sex)

# Remaining samples
print(paste0("Remaining Samples = ", nrow(phen))) # N = 163

## ComBAT does not handle NAs in the methylation file,
## so we have to impute NAs in the methylation beta matrix
## Log transform to normalize the data (we'll reverse that later)
beta <- log((beta/(1-beta)))
beta.imputed <- impute.knn(as.matrix(beta))
beta.imputed <- beta.imputed$data

rm(beta)
gc()

# Run ComBat
SerialParam()
combat_beta <- ComBat(dat = beta.imputed, mod = moddata, batch = chip, BPPARAM = SerialParam())
combat_beta <- ComBat(dat = combat_beta, mod = moddata, batch = position, BPPARAM = SerialParam())

# Reverse Beta Values
reversbeta <- 1/(1+(1/exp(combat_beta)))

## We need to put NAs back (from the original matrix) to ComBAT adjusted beta matrix
## We don't want to use imputed beta values for missing data
norm <- fread("noob_qcd_crossReactiveProbesRemoved.csv", data.table = F)
rownames(norm) <- norm$V1
norm<-norm[,-1]

beta <- as.data.frame(reversbeta)

norm <- norm[, names(norm) %in% names(beta)]
norm <- norm[, order(names(norm))]
beta <- beta[, order(names(beta))]
table(names(norm) == names(beta))

beta[is.na(norm)] <- NA
beta[beta >= 0.9999998] <- 0.9999999

write.csv(beta,file="noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_sex.csv",quote=FALSE,row.names=TRUE)
# This is the final beta values that you'll use in EWAS
