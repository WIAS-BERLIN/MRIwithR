if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","book_init.R"))

## ----"Kirby data", echo=TRUE, eval=FALSE---------------------------------------------------------------
 library(kirby21.base)
 download_kirby21_data(modality = "T1",
   outdir = file.path(dataDir, "Kirby21", "T1"))
 download_kirby21_data(modality = "T2",
   outdir = file.path(dataDir, "Kirby21", "T2"))
 download_kirby21_data(modality = "FLAIR",
   outdir = file.path(dataDir, "Kirby21", "FLAIR"))


## ----"Download fmri data 1", echo=TRUE, eval=FALSE-----------------------------------------------------
 if(!Sys.which("aws") == ""){
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


## ----"Download fmri data 2", echo=TRUE, eval=FALSE-----------------------------------------------------
 dn02 <- file.path("ds000117", "sub-01",
                   "ses-mri", "func")
 dn2 <- file.path(dataDir, dn02)
 cmd <- paste0("aws s3 sync --no-sign-request ",
               "s3://openneuro.org/", dn02, " ", dn2)
 system(cmd)


## ----"Download MyConnectome data", echo=TRUE, eval=FALSE-----------------------------------------------
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

## ----"DICOM example data", echo=TRUE, eval=FALSE-------------------------------------------------------
 dn5 <-file.path(dataDir,"example-dicom")
 cmd <- paste("git clone",
   "https://github.com/datalad/example-dicom-structural",
   dn5)
 system(cmd)


## ----"read MPM example data", echo=TRUE, eval=FALSE----------------------------------------------------
 srcurl <- paste0("https://owncloud.gwdg.de/index.php/",
                  "s/iv2TOQwGy4FGDDZ/download?path=",
                  "%2F&files=hmri_sample_dataset.zip")
 download.file(srcurl, file.path(dataDir, "MPM.zip"),
               method="curl")
 unzip(file.path(dataDir, "MPM.zip"), exdir=dataDir)
 file.remove(file.path(dataDir, "MPM.zip"))


## ----"MNI atlas overview", echo=TRUE, eval=FALSE-------------------------------------------------------
 library(mni)


## ----"read MNI_ICBN_lin atlas data 2", echo=TRUE, eval=FALSE-------------------------------------------
 srcurl <- mni_datasets("nifti")[[5]][1]
 download.file(srcurl,
               file.path(dataDir,
               "mni_icbm152_lin_nifti.zip"))
 unzip(file.path(dataDir, "mni_icbm152_lin_nifti.zip"),
       exdir=file.path(dataDir,"mni"))
 file.remove(file.path(dataDir,
                       "mni_icbm152_lin_nifti.zip"))


## ----"read shen268 atlas data", echo=TRUE, eval=FALSE--------------------------------------------------
 atlasd <- file.path(dataDir,"atlas")
 if(!dir.exists(atlasd)) dir.create(atlasd)
 fnshen <- file.path(atlasd,
                  "shen_1mm_268_parcellation.nii.gz")
 download.file(paste0("https://www.nitrc.org/frs/",
     "download.php/7976/shen_1mm_268_parcellation",
     ".nii.gz"), fnshen)

