#########################################################################
### PGC Pipeline - Script 4 - EWAS
#########################################################################
#' Clean
rm(list=ls())
gc()

#' Load packages
load_pkgs <- c("CpGassoc", "data.table", "tibble", "feather")
lapply(load_pkgs, require, character.only = TRUE)

## Load Methylation Data
beta <- fread("G:/DNHS_2nd_Batch(using PGC pipeline)/noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_sex.csv", data.table = F)
rownames(beta)<-beta$V1
beta<-beta[,-1]

## Load phenotype file with MethylationID (combination of Sentrix ID and Sentrix Position) as the 1st column
pheno <- read.csv("G:/DNHS_2nd_Batch(using PGC pipeline)/Pheno_PCs_QC.csv", row.names = 1)

# A function here to get and order the required data
clean_order <- function(beta, pheno){
  cpg <- beta[, colnames(beta) %in% row.names(pheno)]
  cpg <- cpg[, order(colnames(cpg))]
  pheno <- pheno[rownames(pheno) %in% colnames(cpg), ]
  pheno <- pheno[order(rownames(pheno)), ]
  print(table(rownames(pheno) == colnames(cpg))) # should be TRUE
  stopifnot(all(rownames(pheno) == colnames(cpg)))
  return(list(pheno = pheno, cpg = cpg))
}

cleaned_df <- clean_order(beta = beta, pheno = pheno)
pheno <- cleaned_df$pheno
write.csv(pheno, file = "Pheno_QC_EWAS.csv",row.names = T)

## Define variables
study <- "DNHS" # name of the study, e.g. "GTP", "DNHS" etc.
ptsdVar <- "PTSDpm" # name of the ptsd variable, coded as: cases = 1 and controls = 0
ptsd <- pheno[, ptsdVar, FALSE]

## Define covariates to be adjusted for EWAS
## Covariates to be included:
##  - cell types from step 3 ("CD8T","CD4T","NK","Bcell","Mono")
##  - GWAS PCs PC1 and PC2 (if available), if not mPC2 (Comp.2) and mPC3 (Comp.3) from step 3.1
##  - age
##  - sex (if applicable)
covar <- data.frame(pheno[,c("Comp.2","Comp.3","CD8T","CD4T","NK","Bcell","Mono","female","Age")])

# Function to run EWAS with CpGAssoc
cpg_assoc_test <- function(cpg, pheno, covar){
  message("Running test, patience ...")
  test <- cpg.assoc(cpg, pheno, covar, logit.transform = T, large.data=TRUE)
  assoc <- test$results
  eff <- test$coefficients
  results <- cbind(assoc, eff)
  return(list(results = results, test = test))
}

# Run EWAS with CpGAssoc
results <- cpg_assoc_test(cpg = cleaned_df$cpg,
                          pheno = cleaned_df$pheno[, ptsdVar],
                          covar = covar)

save(results, file = paste0("G:/DNHS_2nd_Batch(using PGC pipeline)/", study,"_EWAS.Rdata")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP


#########################################################################
### PGC Pipeline - Script 4.1 - Smoking Sensitivity Analysis
#########################################################################

## Define covariates to be adjusted for EWAS
## Covariates to be included:
##  - cell types from step 3 ("CD8T","CD4T","NK","Bcell","Mono")
##  - GWAS PCs PC1 and PC2 (if available), if not mPC2 (Comp.2) and mPC3 (Comp.3) from step 3.1
##  - age
##  - sex (if applicable)
##  - smoking
covar2 <- data.frame(pheno[,c("Comp.2","Comp.3","CD8T","CD4T","NK","Bcell","Mono","female","Age","SmoS")])

results2 <- cpg_assoc_test(cpg = cleaned_df$cpg,
                          pheno = cleaned_df$pheno[, ptsdVar],
                          covar = covar2)

save(results2, file = paste0("G:/DNHS_2nd_Batch(using PGC pipeline)/",
                             study,"_EWAS_smoking.Rdata")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP


#########################################################################
### PGC Pipeline - Script 4.2 - Sex stratified analysis
#########################################################################
#' Clean
rm(results, results2)
gc()

# ## Define variables
# study <- "DNHS" # name of the study, e.g. "GTP", "DNHS" etc.
# ptsdVar <- "PTSDpm" # name of the ptsd variable, coded as: cases = 1 and controls = 0
sex <- "female" # name of the sex variable, coded as: females = 1 and males = 0

# Males -------------------------------------------------------------------
pheno.male <- pheno[pheno[sex] == 0,]

# Order Data
cleaned_df_male <- clean_order(beta = beta, pheno = pheno.male)
ptsd.male <- pheno.male[, ptsdVar, FALSE]

## Define covariates to be adjusted for EWAS
## Covariates to be included:
##  - cell types from step 3 ("CD8T","CD4T","NK","Bcell","Mono")
##  - GWAS PCs PC1 and PC2 (if available), if not mPC2 (Comp.2) and mPC3 (Comp.3) from step 3.1
##  - age
covar.male<-data.frame(pheno.male[,c("Comp.2","Comp.3","CD8T","CD4T","NK","Bcell","Mono","Age")])

# Run EWAS
results.male <- cpg_assoc_test(cpg = cleaned_df_male$cpg,
                          pheno = cleaned_df_male$pheno[, ptsdVar],
                          covar = covar.male)

save(results.male, file = paste0("G:/DNHS_2nd_Batch(using PGC pipeline)/", study,"_EWAS_males.Rdata")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP


# Females -------------------------------------------------------------------
pheno.female <- pheno[pheno[sex] == 1,]

# Order Data
cleaned_df_female <- clean_order(beta = beta, pheno = pheno.female)
ptsd.female <- pheno.female[, ptsdVar, FALSE]

## Define covariates to be adjusted for EWAS
## Covariates to be included:
##  - cell types from step 3 ("CD8T","CD4T","NK","Bcell","Mono")
##  - GWAS PCs PC1 and PC2 (if available), if not mPC2 (Comp.2) and mPC3 (Comp.3) from step 3.1
##  - age
covar.female<-data.frame(pheno.female[,c("Comp.2","Comp.3","CD8T","CD4T","NK","Bcell","Mono","Age")])

## Run EWAS with CpGAssoc
# test.female <- cpg.assoc(cpg.female, pheno.female[, ptsdVar], covar.female, logit.transform = T, large.data=TRUE)
# assoc.female <- test.female$results
# eff.female <- test.female$coefficients
# results.female <- cbind(assoc.female,eff.female)
results.female <- cpg_assoc_test(cpg = cleaned_df_female$cpg,
                               pheno = cleaned_df_female$pheno[, ptsdVar],
                               covar = covar.female)

save(results.female, file = paste0("G:/DNHS_2nd_Batch(using PGC pipeline)/", study,"_EWAS_females.Rdata")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP

#########################################################################
### PGC Pipeline - Script 4.3 - Race stratified analysis
#########################################################################
#' Clean
rm(results.male, results.female)
gc()

## Define variables
#study <- " " # name of the study, e.g. "GTP", "DNHS" etc.
#ptsdVar <- " " # name of the ptsd variable, coded as: cases = 1 and controls = 0
race <- " " # name of the race variable

# Europeans (White) -------------------------------------------------------------------
pheno.EA <- pheno[pheno[race] == "European",] # value for Europeans

# Order Data
cleaned_df_EA <- clean_order(beta = beta, pheno = pheno.EA)
ptsd.EA <- pheno.EA[, ptsdVar, FALSE]

## Define covariates to be adjusted for EWAS
## Covariates to be included:
##  - cell types from step 3 ("CD8T","CD4T","NK","Bcell","Mono")
##  - GWAS PCs PC1 and PC2 (if available), if not mPC2 (Comp.2) and mPC3 (Comp.3) from step 3.1
##  - age
##  - sex (if applicable)
covar.EA<-data.frame(pheno.EA[,c("C1","C2","CD8T","CD4T","NK","Bcell","Mono","Age","Sex")])

## Run EWAS with CpGAssoc
results.EA <- cpg_assoc_test(cpg = cleaned_df_EA$cpg,
                                 pheno = cleaned_df_EA$pheno[, ptsdVar],
                                 covar = covar.EA)

save(results.EA, file = paste(study,"EWAS_EA.Rdata", sep = "_")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP

# African Americans (Black) -------------------------------------------------------------------
pheno.AA <- pheno[pheno[race] == 5,] # value for African Americans

# Order Data
cleaned_df_AA <- clean_order(beta = beta, pheno = pheno.AA)
ptsd.AA <- pheno.AA[, ptsdVar, FALSE]

## Define covariates to be adjusted for EWAS
## Covariates to be included:
##  - cell types from step 3 ("CD8T","CD4T","NK","Bcell","Mono")
##  - GWAS PCs PC1 and PC2 (if available), if not mPC2 (Comp.2) and mPC3 (Comp.3) from step 3.1
##  - age
##  - sex (if applicable)
covar.AA<-data.frame(pheno.AA[,c("C1","C2","CD8T","CD4T","NK","Bcell","Mono","Age","Sex")])

## Run EWAS with CpGAssoc
results.AA <- cpg_assoc_test(cpg = cleaned_df_AA$cpg,
                                 pheno = cleaned_df_AA$pheno[, ptsdVar],
                                 covar = covar.AA)

save(results.AA, file = paste(study,"EWAS_AA.Rdata", sep = "_")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP

# Hispanic (Latino) -------------------------------------------------------------------
pheno.Lat <- pheno[pheno[race] == 1,] # value for Hispanic

# Order Data
cleaned_df_Lat <- clean_order(beta = beta, pheno = pheno.Lat)
ptsd.Lat <- pheno.Lat[, ptsdVar, FALSE]

## Define covariates to be adjusted for EWAS
## Covariates to be included:
##  - cell types from step 3 ("CD8T","CD4T","NK","Bcell","Mono")
##  - GWAS PCs PC1 and PC2 (if available), if not mPC2 (Comp.2) and mPC3 (Comp.3) from step 3.1
##  - age
##  - sex (if applicable)
covar.Lat<-data.frame(pheno.Lat[,c("C1","C2","CD8T","CD4T","NK","Bcell","Mono","Age","Sex")])

## Run EWAS with CpGAssoc
results.Lat <- cpg_assoc_test(cpg = cleaned_df_Lat$cpg,
                                 pheno = cleaned_df_Lat$pheno[, ptsdVar],
                                 covar = covar.Lat)

save(results.Lat, file = paste(study,"EWAS_Lat.Rdata", sep = "_")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP

