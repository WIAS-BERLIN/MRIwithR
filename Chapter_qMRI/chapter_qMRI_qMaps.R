source(file.path("..","CodeSecondEdition","chapter_qMRI_initMPM.R"))

if(file.exists(file.path(MPMresDir, "modelMPM.rsc"))) {
   load(file.path(MPMresDir, "modelMPM.rsc"))
} else source(file.path("..","CodeSecondEdition","chapter_qMRI_estESTATICS.R"))


qMRImaps <- calculateQI(modelMPM, TR2 = 3.4)

b1prefix <- "anon_s2018-02-28_18-26-184837-00001-"
B1mapRegName <- paste0(b1prefix, "00001-1_B1mapReg.nii")
B1mapRegFile <- file.path(MPMresDir, pdDir, "Results",
                          "Supplementary", B1mapRegName)
if(!file.exists(B1mapRegFile)) source("chapter_qMRI_registerB1.R")

qMRImapsB1C <- calculateQI(modelMPM,
                        b1File = B1mapRegFile,
                        TR2 = 3.4)

qMRImapsQLB1C <- calculateQI(modelMPMQL,
                          b1File = B1mapRegFile,
                          TR2 = 3.4)

if(file.exists(file.path(MPMresDir, "smoothModelMPM.rsc"))) {
   load(file.path(MPMresDir, "smoothModelMPM.rsc"))
} else source(file.path("..","CodeSecondEdition","chapter_qMRI_smoothESTATICS.R"))

qMRISmoothedMaps <- calculateQI(modelMPMs, b1File =
                             B1mapRegFile, TR2 = 3.4)

qMRISmoothedp1Maps <- calculateQI(modelMPMsp1,
                                  b1File = B1mapRegFile,
                                  TR2 = 3.4)

save(qMRImaps, qMRImapsB1C, qMRImapsQLB1C, qMRISmoothedMaps, qMRISmoothedp1Maps, file=file.path(MPMresDir, "qMaps.rsc"))
