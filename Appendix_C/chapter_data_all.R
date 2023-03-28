

  # ----"set codeDir", echo=FALSE, eval=TRUE--------------------------------------------------------------
if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","book_init.R"))


  # ----"read code", echo=TRUE, eval=FALSE----------------------------------------------------------------
 setwd(baseDir)
 srcurl <- file.path("https://archive.wias-berlin.de",
                     "servlets","MCRFileNodeServlet",
                     "wias_derivate_00003583",
                     "MRIwithRcode.tgz")
 download.file(srcurl, file.path(baseDir,
                    "MRIwithRcode.tgz"), method="curl")
 untar(file.path(baseDir, "MRIwithRcode.tgz"))
 file.remove(file.path(baseDir, "MRIwithRcode.tgz"))


  # ----"Rbook example code", echo=TRUE, eval=FALSE-------------------------------------------------------
 cmd <- paste("git clone",
              "https://github.com/WIAS-BERLIN/MRIwithR")
 system(cmd)


  # ----"Install packages from CRAN",eval=FALSE,echo=TRUE-------------------------------------------------
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
 # Install ANTsR
 mydeps <- c( "Rcpp", "tools", "methods", "magrittr" )
 install.packages( pkgs = mydeps, dependencies = TRUE )
 library(devtools)
 install_github("stnava/cmaker")
 install_github("stnava/ITKR")
 install_github("stnava/ANTsR")


  # ----"set dataDir", echo=TRUE, eval=FALSE--------------------------------------------------------------
 codeDir <- file.path(baseDir, "MRIwithR")
 if(!dir.exists(codeDir))
   stop(paste("incorrect baseDir:",
           baseDir, "directory MRIwithR not found"))
 dataDir <- file.path(baseDir,"MRIwithRdata")
 if(!dir.exists(dataDir)) dir.create(dataDir)
 resDir <- file.path(baseDir,"MRIwithRresults")
 if(!dir.exists(resDir)) dir.create(resDir)


  # ----"Kirby data", echo=TRUE, eval=FALSE---------------------------------------------------------------
 library(kirby21.base)
 download_kirby21_data(modality = "T1",
   outdir = file.path(dataDir, "Kirby21", "T1"))
 download_kirby21_data(modality = "T2",
   outdir = file.path(dataDir, "Kirby21", "T2"))
 download_kirby21_data(modality = "FLAIR",
   outdir = file.path(dataDir, "Kirby21", "FLAIR"))


 if(!Sys.which("aws") == ""){
  # ----"Download fmri data 1", echo=TRUE, eval=FALSE-----------------------------------------------------
 dn01a <- file.path("ds000105", "sub-1", "func")
 dn1a <- file.path(dataDir, dn01a)
 cmd <- paste0("aws s3 sync --no-sign-request ",
               "s3://openneuro.org/", dn01a, " ", dn1a)
 system(cmd)
 dn01b <- file.path("ds000105", "sub-1", "anat")
 dn1b <- file.path(dataDir, dn01b)
 cmd <- paste0("aws s3 sync --no-sign-request ",
               "s3://openneuro.org/", dn01b, " ", dn1b)
 system(cmd)


  # ----"Download fmri data 2", echo=TRUE, eval=FALSE-----------------------------------------------------
 dn02 <- file.path("ds000117", "sub-01",
                   "ses-mri", "func")
 dn2 <- file.path(dataDir, dn02)
 cmd <- paste0("aws s3 sync --no-sign-request ",
               "s3://openneuro.org/", dn02, " ", dn2)
 system(cmd)


  # ----"Download MyConnectome data", echo=TRUE, eval=FALSE-----------------------------------------------
 dn04 <- file.path("sub-01", "ses-WashU")
 dn4 <- file.path(dataDir,"MyConnectome", dn04)
 cmd <- paste0("aws s3 sync --no-sign-request ",
         "s3://openneuro.org/ds000031/", dn04, " ", dn4)
 system(cmd)
 dn04 <- file.path("sub-01", "ses-2015")
 dn4 <- file.path(dataDir,"MyConnectome", dn04)
 cmd <- paste0("aws s3 sync --no-sign-request ",
         "s3://openneuro.org/ds000031/", dn04, " ", dn4)
 system(cmd)
} else {
   warning("Please install the Amazon Web Services CLI tools.")
}


  # ----"DICOM example data", echo=TRUE, eval=FALSE-------------------------------------------------------
 dn5 <-file.path(dataDir,"example-dicom")
 cmd <- paste("git clone",
   "https://github.com/datalad/example-dicom-structural",
   dn5)
 system(cmd)


  # ----"read MPM example data", echo=TRUE, eval=FALSE----------------------------------------------------
 srcurl <- paste0("https://owncloud.gwdg.de/index.php/",
                  "s/iv2TOQwGy4FGDDZ/download?path=",
                  "%2F&files=hmri_sample_dataset.zip")
 download.file(srcurl, file.path(dataDir, "MPM.zip"),
               method="curl")
 unzip(file.path(dataDir, "MPM.zip"), exdir=dataDir)
 file.remove(file.path(dataDir, "MPM.zip"))


  # ----"MNI atlas overview", echo=TRUE, eval=FALSE-------------------------------------------------------
 library(mni)
 mni_datasets()


  # ----"read MNI_ICBN_lin atlas data 2", echo=TRUE, eval=FALSE-------------------------------------------
 srcurl <- mni_datasets("nifti")[[5]][1]
 download.file(srcurl,
               file.path(dataDir,
               "mni_icbm152_lin_nifti.zip"))
 unzip(file.path(dataDir, "mni_icbm152_lin_nifti.zip"),
       exdir=file.path(dataDir,"mni"))
 file.remove(file.path(dataDir,
                       "mni_icbm152_lin_nifti.zip"))


  # ----"read shen268 atlas data", echo=TRUE, eval=FALSE--------------------------------------------------
 atlasd <- file.path(dataDir,"atlas")
 if(!dir.exists(atlasd)) dir.create(atlasd)
 fnshen <- file.path(atlasd,
                  "shen_1mm_268_parcellation.nii.gz")
 download.file(paste0("https://www.nitrc.org/frs/",
     "download.php/7976/shen_1mm_268_parcellation",
     ".nii.gz"), fnshen)


  # ----"read intermediate results", echo=TRUE, eval=FALSE------------------------------------------------
 setwd(baseDir)
 srcurl <- file.path("https://archive.wias-berlin.de",
                     "servlets","MCRFileNodeServlet",
                     "wias_derivate_00003583",
                     "MRIwithRresults.tgz")
 download.file(srcurl, file.path(baseDir,
               "MRIwithRresults.tgz"), method="curl")
 untar(file.path(baseDir, "MRIwithRresults.tgz"))
 file.remove(file.path(baseDir, "MRIwithRresults.tgz"))

