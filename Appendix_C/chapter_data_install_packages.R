if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","book_init.R"))

 # Install from CRAN
 fromCRAN <- c("knitr", "fslr", "oro.nifti", "oro.dicom",
   "dti", "aws", "adimpro", "jsonlite", "rgl",
   "igraph", "mritc", "fmri", "neuRosim", "glasso",
   "XML", "qMRI", "stringr", "KernSmooth", "remotes",
   "pixmap", "misc3d", "fastICA")
 update.packages(ask=FALSE)
 installed <- installed.packages()[,"Package"]
 notinstalled <- !fromCRAN%in%installed
 if(any(notinstalled)) install.packages(
               fromCRAN[!fromCRAN%in%installed])
 # Install from NeuroConductor
 source("https://neuroconductor.org/neurocLite.R")
 neuro_install(c("kirby21.base", "kirby21.t1",
               "kirby21.t2", "kirby21.flair", "mni"))
 # Install IRKR and  ANTsR
 mydeps <- c( "Rcpp", "tools", "methods", "magrittr" )
 install.packages( pkgs = mydeps, dependencies = TRUE )
 library(devtools)
 install_github("stnava/cmaker")
 install_github("stnava/ITKR")
 install_github("stnava/ANTsR")

