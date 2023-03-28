if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","book_init.R"))
library(oro.nifti)
library(aws)
library(adimpro)
rimage.options(ylab = "z",
               zquantile = c(0.001, 0.999))
setCores(8) # number of compute cores to be used

MPMdataDir <- file.path(dataDir,"MPM")
npresDir <- file.path(resDir,"nparam")
if(!dir.exists(npresDir)) dir.create(npresDir)

t1Name <- file.path(MPMdataDir,
    "t1w_mfc_3dflash_v1i_R4_0015",
    "anon_s2018-02-28_18-26-190921-00001-00224-1.nii")
T1 <- as.array(readNIfTI(t1Name, reorient = FALSE))

setCores(16)