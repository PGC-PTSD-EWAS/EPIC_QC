# EPIC_QC

The scripts are used to run cross-sectional EWAS for EPIC array using PGC pipeline. Following are the steps.

### 1.Cleaning data and performing quality control checks. 
**Script 1: PGCpipeline-1-QualityChecks.R**. 
Numerous quality control checks such as Illumina type 17 metrics (described in [BeadArray Controls Reporter Software Guide from Illumina](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium_hd_methylation/beadarray-controls-reporter-user-guide-1000000004009-00.pdf)), 
genotype calling for duplication and outlier check based on three categories **(AA, AB, BB)** were performed. 
Probes which didn't fit in these three categories were flaged as outliers.


### 2. Preprocessing, Normalization and ComBat adjustment. 
**Script 2: PGCpipeline-2-Preprocessing-Normalization.R**  
In this step, data preprocessing and normalization were performed. Quality control for sex discordance, removing probes and samples with > 10% missing values 
and detection p-values cutoff was set to 0.01. Also, [cross hybridized probes](http://www.sciencedirect.com/science/article/pii/S221359601630071X) were removed. 
Then beta values were imputed followed by ComBAT adjustment for chip id and position preserving variation for PTSD, age, and sex (if relevant). After comBAT adjustment, missing values were inserted back at the appropriate locations.


### 3.Cell proportion estimation and calculating principal components and smoking scores.  
**Script 3: PGCpipeline-3-CellCountsCalculation_forEPIC.R**  
We used [IDOL algorithm](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0943-7) 
to estimate the cell proportions.  
**Script 3.1: PGCpipeline-3.1-mPCs.R**.   
Principal components were calculated as described in [Barfield et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4090102/). If you have GWAS data, it is recommended to use that for calculating PCs. In that case you don't need to use *Barfield* approach.  
**Script 3.2: PGCpipeline-3.2-smoking.R**    
This script uses 39 probes to estimate a smoking score.


### 4.CpG association check (EWAS).  
**Script 4: PGCpipeline-4-EWAS.R**  
CpG association analysis was performed on the normalized and combat adjusted data. **PTSD** was
the response variable, methylation of individual CpG probes were the predictor of interest, and **PC2, PC3, CD8T, CD4T, NK, Bcell, Mono, Gender, Age ** were used as covariates.  
*4.1. Sensitivity analysis for smoking.* EWAS was performed by adding **Smoking Score** as a covariate.  
*4.2. Sex stratified analysis.* EWAS was run separetly for males and females.  
*4.2. Race stratified analysis.* EWAS was run separetly for Europens, African Americans and Hispanic/Latinos.  

### 5. Cell-specific differential methylation.
**Script 5: PGCpipeline-5-CellSpecificEWAS.R**  
To identify cell-specific differences, we used the [TOAST R package](https://www.bioconductor.org/packages/release/bioc/html/TOAST.html), using **PTSD** as the main effect and **Age/Sex** as covariates. 

### 6. Rmarkdown Summary.
**Script 6: PGC_Analysis_Summary.Rmd**  
R Markdown document for the summary of cross-sectional analysis using PGC pipeline.
