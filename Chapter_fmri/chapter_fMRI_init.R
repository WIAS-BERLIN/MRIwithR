if(!exists("baseDir")) baseDir <- "../.."
if(!dir.exists(baseDir)) stop("please define baseDir")
dataDir <- file.path(baseDir,"data")
resDir <- file.path(baseDir,"results")
fMRIresDir <- file.path(resDir,"fMRI")
if(!dir.exists(fMRIresDir)) dir.create(fMRIresDir)
#  load packages and define filenames 

library(knitr)
library(adimpro)
haveANTsR <- require(ANTsR)
library(aws)
library(fmri)
haveFSLr <- require(fslr)
library(glasso)
library(igraph)
library(mritc)
library(neuRosim)
library(oro.nifti)
library(rgl)
library(XML)


f117df <- file.path("ds000117", "sub-01", "ses-mri", "func")
f117s1r1 <- file.path(dataDir, f117df,
 "sub-01_ses-mri_task-facerecognition_run-01_bold.nii")
if(!dir.exists(file.path(resDir,f117df))){
  dir.create(file.path(resDir, "ds000117"))
  dir.create(file.path(resDir, "ds000117", "sub-01"))
  dir.create(file.path(resDir, "ds000117", "sub-01",
                       "ses-mri"))
  dir.create(file.path(resDir, f117df))
}
f117s1r1stc <- file.path(resDir, f117df, paste0("stc",
    "sub-01_ses-mri_task-facerecognition_run-01_bold"))

f <- "sub-1_task-objectviewing_run-01_bold.nii.gz"
f105df <- file.path("ds000105", "sub-1", "func")
a105df <- file.path("ds000105", "sub-1", "anat")
f105s1r1 <- file.path(dataDir, f105df, f)
if(!dir.exists(file.path(resDir, "ds000105"))){
  dir.create(file.path(resDir, "ds000105"))
  dir.create(file.path(resDir, "ds000105", "sub-1"))
  dir.create(file.path(resDir, "ds000105", "sub-1",
                       "func"))
  dir.create(file.path(resDir, "ds000105", "sub-1",
                       "anat"))
}
f105dfwr <- file.path(resDir, f105df)
f105s1r1mc <- file.path(f105dfwr, paste0("mc", f))
f105s1r1amc <- file.path(f105dfwr, paste0("amc", f))
if(!dir.exists(file.path(resDir, "fMRI"))) dir.create(file.path(resDir, "fMRI"))
f105s1r1motion <- file.path(resDir, "fMRI", "mocoparams.rsc")
f105s1r1mc <- file.path(resDir, f105df, paste0("mc",f))
f105s1r1amc <- file.path(resDir, f105df, paste0("amc",f))
f105s1r1p <- file.path(resDir, f105df, paste0("p",f))
f105s1T1 <- file.path(resDir, a105df, "sub-1_T1w.nii.gz")
f105s1bold2T1 <- file.path(resDir, a105df, "sub-1_bold2T1w.nii.gz")
fileMNI <- file.path(dataDir, "mni",
                     "icbm_avg_152_t1_tal_lin.nii")
f105s1r1nmc <- file.path(resDir, f105df, paste0("nmc", f))
fileshen268 <- file.path(dataDir, "atlas",
                    "shen_1mm_268_parcellation.nii.gz")
if(!dir.exists(file.path(resDir, "atlas"))) dir.create(file.path(resDir, "atlas"))
fileshen268mni <- file.path(resDir, "atlas",
                    "shen_mni_268_parcellation.nii.gz")
f105s1r1atlas268 <- file.path(resDir, f105df, paste0("l268", f))
f105s1boldmask <- file.path(resDir, a105df, "sub-1_boldmask.nii.gz")
f105s1brain <- file.path(resDir, a105df, "sub-1_brain.nii.gz")
f105s1mask <- file.path(resDir, a105df, "sub-1_brain_mask.nii.gz")
f105s1brainBET <- file.path(resDir, a105df, "fslbet_brain")
f105s1brainFAST <- file.path(resDir, a105df, "fslfast_brain")
f1ev <- "sub-1_task-objectviewing_run-01_events.tsv"
f105s1r1events <- file.path(dataDir, f105df, f1ev)

if(!dir.exists(file.path(resDir, "MyConnectome"))) dir.create(file.path(dataDir, "MyConnectome"))
if(!dir.exists(file.path(resDir, "MyConnectome","sub-01")))
 dir.create(file.path(resDir, "MyConnectome","sub-01"))
if(!dir.exists(file.path(resDir, "MyConnectome","sub-01", "ses-105")))
  dir.create(file.path(resDir, "MyConnectome","sub-01", "ses-105"))
if(!dir.exists(file.path(resDir, "MyConnectome","sub-01", "ses-105", "func")))
  dir.create(file.path(resDir, "MyConnectome","sub-01", "ses-105", "func"))
MCf <- file.path("MyConnectome", "sub-01", "ses-105", "func")
g <- "sub-01_ses-105_task-rest_run-010_bold.nii.gz"
frest <- file.path(dataDir, MCf, g)
fprest <- file.path(resDir, MCf, paste0("prmc", g))
fprestmask <- file.path(resDir, MCf,
        "sub-01_ses-105_task-rest_run_010_mask.nii.gz")
fileMNI2 <- file.path(fsl_std_dir(),
                     "MNI152_T1_2mm_brain_mask.nii.gz")
fr0 <- "sub-01_ses-105_task-rest_run-001_bold"
mfprest0 <- file.path(resDir, MCf, paste0("meanpmc",fr0))
## writeNIfTI requires file names without extension
fr0 <- paste0(fr0,".nii.gz")
mnirestfile <- file.path(resDir, MCf, paste0("nprmc", fr0))
mfprest <- paste0(mfprest0,".nii.gz")
fprest1 <- file.path(resDir, MCf, paste0("prmc", fr0))
mnirestfile <- file.path(resDir, MCf, paste0("nprmc", fr0))
fnTalairach <- file.path(fsl_atlas_dir(),
       "Talairach","Talairach-labels-2mm.nii.gz")
