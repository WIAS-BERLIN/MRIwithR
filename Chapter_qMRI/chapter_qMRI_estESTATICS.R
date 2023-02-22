source(file.path("..","CodeSecondEdition","chapter_qMRI_initMPM.R"))

if(file.exists(file.path(MPMresDir, "MPMdata.rsc"))) {
   load(file.path(MPMresDir, "MPMdata.rsc"))
} else source(file.path("..","CodeSecondEdition","chapter_qMRI_readMPM.R"))

mpm <- readMPMData(t1Files, rpdFiles, rmtFiles,
                   maskFile, TR = mpm$TR, TE = mpm$TE,
                   FA = mpm$FA, verbose = FALSE)

modelMPM <- estimateESTATICS(mpm, verbose = TRUE)

img <- array(0,dim(mpm$mask))
img[mpm$mask] <- mpm$ddata[1, ]
sigma <- aws::awsLocalSigma(img, 16, mpm$mask,
                            ncoils=1, vext= c(1, 1, 1))
sigma <- mean(sigma$sigma[mpm$mask])

modelMPMQL <- estimateESTATICS(mpm, method = "QL",
                               sigma = sigma, L = 1)

save(sigma, modelMPM, modelMPMQL, file=file.path(MPMresDir, "modelMPM.rsc"))
