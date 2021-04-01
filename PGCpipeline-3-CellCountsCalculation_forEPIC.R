#########################################################################
### PGC Pipeline - Script 3 - Calculate cell type proportions
#########################################################################
#' remove all objects from your workspace
rm(list=ls())
gc()
# path to source file required to install the packages
source("/Volumes/Seyma/Methyl/CrossSectionalEPIC/install_needed_packages.R") # PATH TO SOURCE FILE, this is an example path, change accordingly

#' Install packages if not already installed
need_pkgs <- c("FlowSorted.Blood.EPIC", "EpiDISH")
install_bioconductor_pkgs(pkgs = need_pkgs)

#' Run test
check_installed(pkgs = need_pkgs)

#' Load all packages, if needed
load_pkgs <- c("FlowSorted.Blood.EPIC", "EpiDISH", "data.table", "tibble", "feather")
lapply(load_pkgs, require, character.only = TRUE)

#' This function will calculate the cell proportions and combine it with the phenotype information
#' input: beta values before Combat, phenotype file, reference cell proportions,
#' xcol- column name on which the data should be merged
#' ycol- column name
#' Output: Phenotye information with cell types
estimate_cell_proportion <- function(betavals, phenotyps, cell_prop_ref, xcol, ycol ){
  
  ## Calculate cell types using RPC method
  RPC <- epidish(betavals, as.matrix(cell_prop_ref), method = "RPC")
  
  cellTypes <- as.data.frame(RPC$estF) #RPC count estimates
  
  phen <- merge(phenotyps, cellTypes, by.x = xcol, by.y = ycol, all.x = T)
  
  return(phen)
}

## Load beta values (not ComBAT adjusted beta values, just normalized beta values)
beta <- fread("noob_qcd_crossReactiveProbesRemoved.csv", data.table = F)
rownames(beta) <- beta$V1
beta <- beta[,-1]

## Load ref dataset shared in the pipeline package
load("IDOLOptimizedCpGs_REF.Rdata")

## Load the phenotype file
pheno <- read.csv("Pheno_QC.csv")

## Run the function
output <- estimate_cell_proportion(betavals = beta, phenotyps = pheno, cell_prop_ref = ref, xcol = "SampleID", ycol = 'row.names' )
write.csv(output, file = "Pheno_QC.csv", row.names = FALSE) # Save the phenotype file with cell types


