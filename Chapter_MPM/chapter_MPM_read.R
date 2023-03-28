if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_MPM","chapter_MPM_init.R"))

mpm <- readMPMData(t1Files, pdFiles, mtFiles,
                   maskFile, verbose = FALSE)


mpm <- readMPMData(t1Files, rpdFiles, rmtFiles,
                   maskFile, TR = mpm$TR, TE = mpm$TE,
                   FA = mpm$FA, verbose = FALSE)

save(mpm, file=file.path(MPMresDir, "MPMdata.rsc"))
