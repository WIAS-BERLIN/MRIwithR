source(file.path("..","CodeSecondEdition","chapter_qMRI_initMPM.R"))

if(file.exists(file.path(MPMresDir, "modelMPM.rsc"))) {
   load(file.path(MPMresDir, "modelMPM.rsc"))
} else source(file.path("..","CodeSecondEdition","chapter_qMRI_estESTATICS.R"))

modelMPMs <- smoothESTATICS(modelMPM,
                            kstar = 20,
                            alpha = 0.5,
                            verbose = TRUE)

modelMPMsp1 <- smoothESTATICS(modelMPM,
                              kstar = 20,
                              alpha = 0.5,
                              patchsize=1,
                              verbose = TRUE)

save(modelMPMs, modelMPMsp1, file=file.path(MPMresDir, "smoothModelMPM.rsc"))