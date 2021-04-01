
#' I think it is good idea to make the installation and loading process easy
#' Was 'rlang' package installed on your computer?. I got an error for ewastools without rlang,
#' so I included rlang package to the list
#'

#' For Cran packages
install_cran_pkgs <- function(pkgs){
  not_installed <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(not_installed)){
    install.packages(not_installed, dependencies = TRUE) # you may get a warning for ewas here, just ignore that
    if("ewastools" %in% pkgs & "ewastools" %in% not_installed){
      devtools::install_github("hhhh5/ewastools@master")
    }
  }else{
    message("Cran packages you need are installed")
  }
}

#' For bioconductor packages
install_bioconductor_pkgs <- function(pkgs){
  bioc_not_installed <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(bioc_not_installed)){
    BiocManager::install(c(bioc_not_installed))
  }else{
    message("Bioconductor packages you need are installed")
  }
}




#' Run test to check if all packages are installed
#' If not, figure out the problem and make sure all packages are istalled
check_installed <- function(pkgs){
  not_installed <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(not_installed)){
    message("Package(s) not installed: ", not_installed)
  }else{
    message("All good ...")
  }
}



