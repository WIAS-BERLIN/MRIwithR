if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_MPM","chapter_MPM_init.R"))

if(file.exists(file.path(MPMresDir, "modelMPM.rsc"))) {
   load(file.path(MPMresDir, "modelMPM.rsc"))
} else source(file.path(baseDir,"MRIwithR","Chapter_MPM","chapter_MPM_estESTATICS.R"))

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
