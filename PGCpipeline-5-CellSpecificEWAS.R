########################################################################
### PGC Pipeline - Script 5 - Cell Specific EWAS
#########################################################################

# Load the code to install packages
source("/Volumes/Seyma/Methyl/CrossSectionalEPIC/install_needed_packages.R") # PATH TO SOURCE FILE, this is an example path, change accordingly

#' Install packages if not already installed
bioc_packages <- c("TOAST", "RefFreeEWAS", "EpiDISH",
                    "limma")
cran_packages <- c("feather", "tibble", "tools",
                   "nnls", "data.table")

install_bioconductor_pkgs(pkgs = bioc_packages)
install_cran_pkgs(pkgs = cran_packages)

#' Run test
check_installed(pkgs = c(bioc_packages, cran_packages))

#' Load all packages, if needed
lapply(c(bioc_packages, cran_packages), require, character.only = TRUE)

#' Load beta values
# we need CpGs as rownames, so add the
# column that contains CpGs as rownames in the data
# In our data "V1" is the column containing CpGs
beta <- fread("noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_sex.csv", data.table = F)
beta <- column_to_rownames(beta, var = "V1")

## Load phenotype file with MethylationID (combination of Sentrix ID and Sentrix Position) as the 1st column
pheno <- read.csv("Pheno_QC.csv", row.names = 1)

# Now all the samples in beta and phenotype should match
# if samples are not matching, we will get the samples that
# are in phenotype file only
if(!all(colnames(beta) %in% rownames(pheno)) | !all(rownames(pheno) %in% colnames(beta))){
  beta <- beta[, which(colnames(beta) %in% rownames(pheno))]
  dim(beta)
  pheno <- pheno[rownames(pheno)%in%colnames(beta),]
}else message("Samples in both are matching")

# Now check if the columns in beta and pheno are in same order
# if not we can order them using the following code
table(colnames(beta) == rownames(pheno))
if(!all(colnames(beta) == rownames(pheno))){
  beta <- beta[, order(colnames(beta))]
  pheno <- pheno[order(rownames(pheno)), ]
  stopifnot(all(colnames(beta) == rownames(pheno))) # check the order again
}else message("Data is already ordered")

## Define Variable Names
study <- " " # name of the study, e.g. "GTP", "DNHS" etc.
ptsdVar <- " " # name of the ptsd variable, coded as: cases = 1 and controls = 0
sex <- " " # name of the sex variable
age <- " " # name of the age variable

# We need to make the columns that are categories as factors
# As an example, we have gender/PTSD as categories,
# so we will make them as factors
fact_cols <- c(ptsdVar,sex)
pheno[fact_cols] <- lapply(pheno[fact_cols], factor)

# Now we need pheno and cell proportions separate
# to use them in the model
pheno_final <- pheno[, c(sex,age,ptsdVar)]
head(pheno_final)

cell_prop <- pheno[, c("CD8T", "CD4T", "NK",
                       "Bcell", "Mono", "Neu")] # cell estimations, change names if needed

all(rownames(pheno_final) == rownames(cell_prop)) # checking order, all should be TRUE
all(colnames(beta) == rownames(pheno_final))

# Function to get the significant CpGs
# Only those CpGs with p < 0.05
get_significant <- function(results, cutoff){
  filtered <-  lapply(results, function(x){
    x <- rownames_to_column(x, var = "CpGs")
    x[which(x['fdr'] < cutoff), ]
  })
  return(filtered)
}

# Run TOAST csDM function
design <- makeDesign(pheno_final, Prop = cell_prop)
design$formula

# fit model
fitted_model <- fitModel(design, as.matrix(beta))
res_ptsd <- csTest(fitted_model, coef = ptsdVar)
res_ptsd_sig <-  get_significant(results = res_ptsd, cutoff = 0.05)

# Save the output in the directory
save(res_ptsd, file = paste(study,"csDM_TOAST.Rdata", sep = "_")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP
save(res_ptsd_sig, file = paste(study,"csDM_TOAST_significant(fdr=0.05).Rdata", sep = "_")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP

# ---------------------- END ------------------------------

