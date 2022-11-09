## load packages needed in qMRI chapter
MPMdataDir <- "../data/MPM"
MPMresDir <- "../results/MPM"
if(!dir.exists(MPMresDir)) dir.create(MPMresDir)
IRresDir <- "../results/IR"
if(!dir.exists(IRresDir)) dir.create(IRresDir)

installed <- installed.packages(fields="")
haveFSLr <- require(fslr)
library(qMRI)
haveANTsR <- require(ANTsR)
library(aws)
library(adimpro)
library(KernSmooth)
setCores(4)

t1Dir <- "t1w_mfc_3dflash_v1i_R4_0015"
pdDir <- "pdw_mfc_3dflash_v1i_R4_0009"
mtDir <- "mtw_mfc_3dflash_v1i_R4_0012"
t1Prefix <- "anon_s2018-02-28_18-26-190921-00001-"
pdPrefix <- "anon_s2018-02-28_18-26-185345-00001-"
mtPrefix <- "anon_s2018-02-28_18-26-190132-00001-"
t1ID <- c("00224-1.nii", "00448-2.nii", "00672-3.nii",
          "00896-4.nii", "01120-5.nii", "01344-6.nii",
          "01568-7.nii", "01792-8.nii")
t1Names <- paste0(t1Prefix, t1ID)
t1Files <- file.path(MPMdataDir, t1Dir, t1Names)

pdID <- c("00224-1.nii", "00448-2.nii", "00672-3.nii",
          "00896-4.nii", "01120-5.nii", "01344-6.nii",
          "01568-7.nii", "01792-8.nii")
pdNames <- paste0(pdPrefix, pdID)
pdFiles <- file.path(MPMdataDir, pdDir, pdNames)

mtID <- c("00224-1.nii", "00448-2.nii", "00672-3.nii",
          "00896-4.nii", "01120-5.nii", "01344-6.nii")
mtNames <- paste0(mtPrefix, mtID)
mtFiles <- file.path(MPMdataDir, mtDir, mtNames)

##
##  compute masks using fslbet (package fslr)
##
maskName <- paste0(t1Prefix,
                   "00224-1_brain")
maskFile <- file.path(MPMresDir, t1Dir, maskName)


rpdNames <- paste0(pdPrefix,"r", pdID)
rpdFiles <- file.path(MPMresDir, pdDir, rpdNames)

rmtNames <- paste0(mtPrefix,"r",  mtID)
rmtFiles <- file.path(MPMresDir, mtDir, rmtNames)


