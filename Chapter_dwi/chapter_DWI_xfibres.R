source("../CodeSecondEdition/chapter_DWI_init.R")
gradientFile <- file.path(rdwipd,
                          "sub-01_ses-106_dwi_proc.eddy_rotated_bvecs")
bvalueFile <- file.path(rdwipd,
                        "sub-01_ses-106_dwi_proc.bval")
dataFile <- file.path(rdwipd,
                      "sub-01_ses-106_dwi_proc")
fnmask <- file.path(rdwipd,"sub-01_ses-106_brain_mask")
xfibres(dataFile, gradientFile, bvalueFile, fnmask, 2,
        verbose=FALSE, opts="--model=2")
