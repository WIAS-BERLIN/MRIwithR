##
##   This code currently fails
##
source("../CodeSecondEdition/chapter_DWI_init.R")

resDir <- "../results"
fnmask <- file.path(rdwipd, "sub-01_ses-106_brain.nii.gz")
samples <- file.path(ldwipd, "MyConnectome", "sub-01",
                     "ses-106","logdir")
probtrackx(samples = samples,
           mask = fnmask, seed = fnmask)
