library(oro.nifti)
library(aws)
library(adimpro)
rimage.options(ylab = "z",
               zquantile = c(0.001, 0.999))
setCores(8) # number of compute cores to be used

if(!exists("baseDir")) baseDir <- "../.."
if(!dir.exists(baseDir)) stop("please define baseDir")
MPMdataDir <- file.path(baseDir,"data","MPM")
npresDir <- file.path(baseDir,"results","nparam")
if(!dir.exists(npresDir)) dir.create(npresDir)

t1Name <- file.path(MPMdataDir,
    "t1w_mfc_3dflash_v1i_R4_0015",
    "anon_s2018-02-28_18-26-190921-00001-00224-1.nii")
T1 <- as.array(readNIfTI(t1Name, reorient = FALSE))

