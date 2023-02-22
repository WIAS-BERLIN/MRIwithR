source(file.path("..","CodeSecondEdition","chapter_qMRI_initMPM.R"))

PD <- antsImageRead(pdFiles[1])


## ----label="registerBiasFieldMaps2", echo=TRUE, eval=TRUE---------------------------------------------------
b1prefix <- "anon_s2018-02-28_18-26-184837-00001-"
B1refName <- paste0(b1prefix, "00001-1_B1ref.nii")
B1refFile <- file.path(MPMdataDir, pdDir,
                       "Results", "Supplementary",
                       B1refName)
B1ref <- antsImageRead(B1refFile)


## ----label="registerBiasFieldMaps3", echo=TRUE, eval=TRUE---------------------------------------------------
B1mapName <- paste0(b1prefix, "00001-1_B1map.nii")
B1mapFile <- file.path(MPMdataDir, pdDir,
                       "Results", "Supplementary",
                       B1mapName)
B1map <- antsImageRead(B1mapFile)


## ----label="registerBiasFieldMaps4", echo=TRUE, eval=TRUE---------------------------------------------------
B1refReg <- antsRegistration(fixed = PD,
                             moving = B1ref,
                             typeofTransform = "SyN")
B1refRegName <- paste0(b1prefix,
                      "00001-1_B1refReg.nii")
B1refRegFile <- file.path(MPMresDir, pdDir,
                          "Results", "Supplementary",
                          B1refRegName)
if(!dir.exists(file.path(MPMresDir, pdDir, "Results"))) {
   dir.create(file.path(MPMresDir, pdDir, "Results"))
   dir.create(file.path(MPMresDir, pdDir, "Results", "Supplementary"))
}
antsImageWrite(B1refReg$warpedmovout, B1refRegFile)


## ----label="registerBiasFieldMaps5", echo=TRUE, eval=TRUE---------------------------------------------------
B1mapReg <- antsApplyTransforms(PD,
                                B1map,
                                B1refReg$fwdtransforms,
                                interpolator = "linear")
B1mapRegName <- paste0(b1prefix, "00001-1_B1mapReg.nii")
B1mapRegFile <- file.path(MPMresDir, pdDir,
                          "Results", "Supplementary",
                          B1mapRegName)
antsImageWrite(B1mapReg, B1mapRegFile)
