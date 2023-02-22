source(file.path("..","CodeSecondEdition","chapter_qMRI_initMPM.R"))

mpm <- readMPMData(t1Files, pdFiles, mtFiles,
                   maskFile, verbose = FALSE)


mpm <- readMPMData(t1Files, rpdFiles, rmtFiles,
                   maskFile, TR = mpm$TR, TE = mpm$TE,
                   FA = mpm$FA, verbose = FALSE)

save(mpm, file=file.path(MPMresDir, "MPMdata.rsc"))
