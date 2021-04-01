#########################################################################
### PGC Pipeline - Script 1 - Quality Checks
#########################################################################
# Path to the source file (install_needed_packages.R) required to install the packages
# This install_needed_packages.R file is part of the pipeline package
source("/Volumes/Methyl/install_needed_packages.R") # PATH TO SOURCE FILE, this is an example path, change accordingly  

#' We need these packages
need_packages <- c("rlang", "stringi", "data.table", "purrr", "svd", "devtools","ewastools")

#' install by running the function
install_cran_pkgs(pkgs = need_packages)

#' Check if all you need is installed
check_installed(pkgs = need_packages)

#load all packages, if needed
load_pkgs <-  c("ewastools", "stringi", "data.table", "purrr", "svd")
lapply(load_pkgs, require, character.only = TRUE )

# Step 1) Load idats -----------------------------------------------------------

#' function to get idat files list
#' Input: path to the directory having idats (subfolders of each chip)
#' it will read names of all files, and create the paths to each file
#' Output: list of paths
file_list <- function(path){
  files <- list.files(path, recursive = TRUE, pattern = ".idat")
  files <- gsub(paste(c("_Grn.idat", "_Red.idat"), collapse = "|"), "", files)
  file_paths <- paste0(path, files)
  files <- gsub(".*/","",files)
  return(list(names = unique(files), paths = unique(file_paths)))
}

# Load sample sheet or phenotype file
# This file should include SentrixID and SentrixPosition as column names
# Then create SampleID by combining Sentrix barcode (Sentrix_ID) and Sentrix position (Sentrix_Position)
pheno <- read.csv(" ") # PHENOTYPE FILE - HERE
pheno$SampleID <- paste(pheno$SentrixID, pheno$SentrixPosition,sep = "_") # create SampleID by combining Sentrix barcode and Sentrix position

# Read idats
#' Path where you idat files are
#' It can have many subfolders containing idats
main_dir <- "G:/methylRawData2020/"  # PATH TO IDATS FOLDER HERE, this is an example path, change accordingly  
file_paths <- file_list(path = main_dir)
paths <- file_paths$paths[which(file_paths$names %in% pheno$SampleID)]
stopifnot(all(grepl(paste(c(pheno$SampleID), collapse = "|"), paths))) # stop if all samples are not in your path
meth <- read_idats(paths, quiet = FALSE)

#' Now we will sort the sample sheet according the read idat files
#' And we will check if all are sorted
pheno <- pheno[order(match(pheno$SampleID, meth$meta$sample_id)), ]
stopifnot(all(pheno$SampleID == meth$meta$sample_id))

# Step 2) Control metrics ------------------------------------------------------

# evaluates 17 control metrics which are describe in the BeadArray Controls Reporter Software Guide from Illumina.
# compares all 17 metrics against the thresholds recommended by Illumina.
ctrls <- control_metrics(meth)
pheno$failed = sample_failure(ctrls)
table(pheno$failed) # Samples with failed == TRUE should be identified and removed for further analysis
fails <- subset(pheno, failed == "TRUE")
message("Failed Samples: ", nrow(fails)) # No samples are failed

# Step 3) Genotype calling and outliers ----------------------------------------
# Need to do a quick pre-processing just for that step, but we won't be using these beta values for the further analysis
beta <- meth %>% detectionP %>% mask(0.01) %>% correct_dye_bias %>% dont_normalize
snps <- meth$manifest[probe_type == "rs", index]
snps <- beta[snps, ]

#' Check if you got only snps, if everything is fine, it should not give any error
stopifnot(all(grepl("rs", rownames(snps))))

# estimates the parameters of a mixture model consisting of three Beta distributions representing the heterozygous and the two homozygous genotypes,
# and a fourth component, a uniform distribution, representing outliers.
genotypes <- call_genotypes(snps, learn = FALSE)
# average log odds of belonging to the outlier component across all SNP probes.
# flagging samples with a score greater than -4 for exclusion is recommended, because of possible contamination.
pheno$outlier <- snp_outliers(genotypes)

#' This will raise an error if 'pheno' is not a data table
#' So we need to convert it to data table
if(!is.data.table(pheno)){
  pheno <- data.table(pheno)
}
pheno[outlier > -4, .(SampleID,outlier)]

# Flag the outlier samples: Flagged outlier samples are denoted as "Y", while good samples as "N"
pheno$outlierYN <- pheno$outlier > -4

# Check for duplicated samples
pheno$donor_id <- enumerate_sample_donors(genotypes)

# List duplicates
pheno[, n:=.N, by = donor_id]
pheno[n > 1, .(SampleID,donor_id)] # Check to see if there are any duplicates, and remove idats of the one from the folder

# Save the phenotype file with the QC details
write.csv(pheno, file = "Pheno_QC.csv", row.names = F)


