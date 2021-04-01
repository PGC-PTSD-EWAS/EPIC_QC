## R function location.based.pc() for computing principal components based on CpG sites near genetic variants
### Uses 1000 Genomes Project annotation data and "location-based pruning" approach described in: 
### Barfield et al. (2014) Accounting for population stratification in DNA methylation studies.  Gen Epid, in press. 

### This function does not require user modification. 
### Except in special circumstances, it is recommended to modify the code in "sample_code.R" and leave this function as is.

### beta.obj should be a data.frame or matrix of beta values or M-values with: 
###### one row per CpG site, one column per sample
###### row names = CpG names, column names = sample names  
### file.pathname is a character string providing the pathname for the annotation file
### See sample_code.R for examples.

### location.based.pc function
### This function will select only CpG sites in the location-based annotation file, 
### set missing values to the CpG-site-average, and compute principal components

location.based.pc<-function(beta.obj,file.pathname){
    load(file.pathname)

    beta.obj<-as.matrix(beta.obj)
    gc()
    if(sum(is.na(beta.obj))>0){
      cpgmeans<-t(t(rowMeans(beta.obj,na.rm=T))) %*% matrix(1,ncol=ncol(beta.obj))
      gc()
      missing.sites<-which(is.na(beta.obj))
      gc()
      beta.obj[missing.sites]<-cpgmeans[missing.sites]
      rm(cpgmeans)
      gc()
    }
    beta.obj<-beta.obj[which(rownames(beta.obj)%in%a2[,1]),]
    theresult<-princomp(x=beta.obj,cor=T)
    rm(a2)
    gc()
    return(theresult)
}

### The returned result will be a princomp object, with the principal components stored as object_name$loadings
### See sample_code.R for an example of how to use this function.
