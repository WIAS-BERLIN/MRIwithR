source("../CodeSecondEdition/chapter_DWI_init.R")
## ----label="Reading bvecs", echo=TRUE, eval = TRUE, results=TRUE--------------------------------------------
run1File <- file.path(dwid,
                      "sub-01_ses-106_run-001_dwi.nii.gz")
run2File <- file.path(dwid,
                      "sub-01_ses-106_run-002_dwi.nii.gz")
fbval1 <- file.path(dwid,
                    "sub-01_ses-106_run-001_dwi.bval")
fbval2 <- file.path(dwid,
                    "sub-01_ses-106_run-002_dwi.bval")
bval1 <- scan(fbval1)
bval2 <- scan(fbval2)

fbvec1 <- file.path(dwid,
                    "sub-01_ses-106_run-001_dwi.bvec")
fbvec2 <- file.path(dwid,
                    "sub-01_ses-106_run-002_dwi.bvec")
bvec1 <- as.matrix(read.table(fbvec1))
bvec2 <- as.matrix(read.table(fbvec2))


## ----label="reading unprocessed data1", echo=TRUE, eval = TRUE, results='hide'------------------------------
dwobj1 <- readDWIdata(bvec1,
                      run1File,
                      bvalue = bval1,
                      format = "NIFTI")


## ----label="using afsigmc on unprocessed data", echo=TRUE, eval = TRUE, results='markup',cache=TRUE---------
v1 <- afsigmc(dwobj1@si[ , , , 1], level = 1.75,
              ncoils=1, method="modevn")
v2 <- afsigmc(dwobj1@si[ , , , 1], level = 1,
              ncoils=1, method="modem1chi")
v3 <- afsigmc(dwobj1@si[ , , , 1], level = 1,
              ncoils=1, method="bkm2chi")
v4 <- afsigmc(dwobj1@si[ , , , 1], level = 1,
              ncoils=1, method="bkm1chi")
cat("Estimated sigma: ", signif(c(v1, v2, v3, v4), 3))


## ----label="reading unprocessed data2", echo=TRUE, eval = TRUE, results='hide'------------------------------
dwobj2 <- readDWIdata(bvec2,
                      run2File,
                      bvalue = bval2,
                      format = "NIFTI")


## ----label="AF2013",echo=TRUE,eval = TRUE, cache=TRUE, results='hide'---------------------------------------
afsigma1b0 <- aflsigmc(dwobj1@si[, , , 1],
                       ncoils = 1,
                       level = 1)
afsigma2b0 <- aflsigmc(dwobj2@si[, , , 1],
                       ncoils = 1,
                       level = 1)


## ----label="Figure_5_10",echo=-1,eval = TRUE,cache=TRUE,fig.cap="Local noise estimation using \\code{aflsigmc}.",fig.width=8,fig.height=2.2,fig.align='center'----
par(mfrow = c(1, 4),
    mar = c(1, 1, 2, .1), mgp = c(2, 1, 0),
    xaxt = "n", yaxt = "n")
rimage(afsigma1b0$sigma[, , 41])
title("Effective sigma run 1")
rimage(afsigma1b0$Leff[, , 41])
title("Effective L run 1")
rimage(afsigma2b0$sigma[, , 41])
title("Effective sigma run 2")
rimage(afsigma2b0$Leff[, , 41])
title("Effective L run 2")


## ----"using awslsigmc", eval = TRUE, echo=TRUE, results='hide',cache=TRUE, warning=FALSE, message=FALSE-----
kstar <- 16; L <- 1; setCores(4)
sb1 <- awslsigmc(dwobj1@si[ , , , 1], kstar, ncoils = L)
sb2 <- awslsigmc(dwobj2@si[ , , , 1], kstar, ncoils = L)

save(sb1,sb2,v1,v2,v3,v4,afsigma1b0,afsigma2b0, file=file.path(rdwi,"sigmaest.rsc"))
## ----"using applytopup", eval = FALSE, echo=TRUE------------------------------------------------------------
acqParFileAll <- file.path(resDir, "MyConnectome",
                           "sourcedata", "dicom_headers",
                           "sub-01", "ses-106", "dwi",
                           "datain.txt")
fntopup <- file.path(rdwipd, "sub-01_ses-106_topupout")
fdn1 <- file.path(rdwipd, "tmpn1")
fdn2 <- file.path(rdwipd, "tmpn2")
fnsigma <- file.path(rdwipd, "sub-01_ses-106_sigma")
ds1 <- readNIfTI(run1File)
writeNIfTI(as.nifti(sb1$sigma, ds1), fdn1)
ds2 <- readNIfTI(run2File)
writeNIfTI(as.nifti(sb1$sigma, ds2), fdn2)
fdn1 <- paste0(fdn1,".nii.gz")
fdn2 <- paste0(fdn2,".nii.gz")
applytopup(c(fdn1, fdn2),
           acqParFileAll,
           c(1, 2),
           fntopup,
           fnsigma)
rm(dwobj1, dwobj2, ds1, ds2, sb1, sb2)

