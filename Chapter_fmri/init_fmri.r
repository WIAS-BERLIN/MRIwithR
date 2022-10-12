library(knitr)
library(adimpro)
library(ANTsR)
library(aws)
library(fmri)
library(fslr)
library(glasso)
library(igraph)
library(mritc)
library(neuRosim)
library(oro.nifti)
library(rgl)
library(XML)

dataDir <- "data"
f1 <-
  "sub-01_ses-mri_task-facerecognition_run-01_bold.nii"
f117df <- file.path(dataDir,
               "ds000117", "sub-01", "ses-mri", "func")
f117s1r1 <- file.path(f117df, f1)
f117s1r1stc <- file.path(f117df, paste0("stc", f1))
f2 <- "sub-1_task-objectviewing_run-01_bold.nii.gz"
f105s1r1 <- file.path(dataDir, "ds000105",
                      "sub-1", "func", f2)
f105df <- file.path(dataDir, "ds000105",
                      "sub-1", "func")
f105s1r1mc <- file.path(f105df, paste0("mc", f1))
f105s1r1amc <- file.path(f105df, paste0("amc", f1))
f105s1r1motion <- file.path(f105df, "mocoparams.rsc")
f105s1r1 <- file.path(dataDir, "ds000105",
                      "sub-1", "func", f2)
f105s1r1mc <- file.path(dataDir, "ds000105",
            "sub-1", "func", paste0("mc",f2))
f105s1r1amc <- file.path(dataDir, "ds000105",
            "sub-1", "func", paste0("amc",f2))
f105s1r1motion <- file.path(dataDir, "ds000105",
            "sub-1", "func", "mocoparams.rsc")
f105s1r1p <- file.path(dataDir, "ds000105",
            "sub-1", "func", paste0("p",f2))
f105s1T1 <- file.path(dataDir, "ds000105",
            "sub-1", "anat", "sub-1_T1w.nii.gz")
f105s1bold2T1 <- file.path(dataDir, "ds000105",
            "sub-1", "anat", "sub-1_bold2T1w.nii.gz")
fileMNI <- file.path(dataDir, "mni",
                     "icbm_avg_152_t1_tal_lin.nii")
f105s1r1nmc <- file.path(dataDir, "ds000105",
                     "sub-1","func", paste0("nmc", f2))
fileshen268 <- file.path(dataDir, "atlas",
                       "shen_1mm_268_parcellation.nii")
fileshen268mni <- file.path(dataDir, "atlas",
                       "shen_mni_268_parcellation.nii")
f105s1r1atlas268 <- file.path(dataDir, "ds000105",
                   "sub-1", "func", paste0("l268", f2))
f105s1boldmask <- file.path(dataDir, "ds000105",
              "sub-1", "anat", "sub-1_boldmask.nii.gz")
f105s1brain <- file.path(dataDir, "ds000105",
                        "sub-1", "anat", "sub-1_brain")
f105s1mask <- file.path(dataDir, "ds000105",
                        "sub-1", "anat", "sub-1_mask")
f105s1boldmask <- file.path(dataDir, "ds000105",
              "sub-1", "anat", "sub-1_boldmask.nii.gz")
f105s1brainBET <- file.path(dataDir, "ds000105",
                       "sub-1", "anat", "fslbet_brain")
f105s1brainFAST <- file.path(dataDir, "ds000105",
                      "sub-1", "anat", "fslfast_brain")
f1events <- "sub-1_task-objectviewing_run-01_events.tsv"
f105s1r1events <- file.path(dataDir, "ds000105",
                             "sub-1", "func", f1events)

f3 <- "sub-01_ses-105_task-rest_run-010_bold.nii.gz"
frest <- file.path(dataDir, "MyConnectome",
                       "sub-01", "ses-105", "func", f3)
fprest <- file.path(dataDir, "MyConnectome", "sub-01",
                 "ses-105", "func", paste0("prmc", f3))
fprestmask <- file.path(dataDir, "MyConnectome",
        "sub-01", "ses-105", "func",
        "sub-01_ses-105_task-rest_run_010_mask.nii.gz")
mnirestfile <- file.path(dataDir, "MyConnectome",
      "sub-01", "ses-105", "func", paste0("nprmc", f3))
fr0 <- "meanpmcsub-01_ses-105_task-rest_run-001_bold"
mfrest <- file.path(dataDir, "MyConnectome", "sub-01",
                   "ses-105", "func", fr0)
fnTalairach <- file.path(fsl_atlas_dir(),
       "Talairach","Talairach-labels-2mm.nii.gz")
