#baseDir <- dirname(getwd())
if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_dwi","chapter_DWI_init.R"))

#
#   topup and eddy  run for hours-days 
#   if you don't have FSL installed or if you dont wand to send that much 
#   computing time you may want to use what is in the results dircetory
#
par(mfrow = c(1, 2),
    mar = c(1, 1, 2, 0.1), mgp = c(2, 1, 0))
run1FileName <- "sub-01_ses-2015_acq-g79_dir-pe0_dwi.nii.gz"
run2FileName <- "sub-01_ses-2015_acq-g79_dir-pe1_dwi.nii.gz"
run1File <- file.path(dwid, run1FileName)
run2File <- file.path(dwid, run2FileName)
ds1 <- readNIfTI(run1File)
ds2 <- readNIfTI(run2File)


## ----label="initialize dirs", echo=TRUE, eval = TRUE, results=FALSE, message=FALSE--------------------------
rdwipd <- file.path(resDir, "MyConnectome", "sub-01",
                   "ses-2015", "dwi-proc")
if (!dir.exists(rdwipd)) dir.create(rdwipd)


## ----label="Reading unprocessed data", echo=TRUE, eval=FALSE, results=FALSE---------------------------------
run3File <- file.path(dwid,
                      "sub-01_ses-2015_acq-g81_dir-pe0_dwi.nii.gz")
run4File <- file.path(dwid,
                      "sub-01_ses-2015_acq-g81_dir-pe1_dwi.nii.gz")
ds3 <- readNIfTI(run3File)
ds4 <- readNIfTI(run4File)


## ----label="Reading bvals", echo=TRUE, eval = TRUE, results=TRUE--------------------------------------------
fbval1 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g79_dir-pe0_dwi.bval")
fbval2 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g79_dir-pe1_dwi.bval")
fbval3 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g81_dir-pe0_dwi.bval")
fbval4 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g81_dir-pe1_dwi.bval")
bval1 <- scan(fbval1)
bval2 <- scan(fbval2)
bval3 <- scan(fbval3)
bval4 <- scan(fbval4)
bval1


## ----label="Reading bvecs", echo=TRUE, eval = TRUE, results=TRUE--------------------------------------------
fbvec1 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g79_dir-pe0_dwi.bvec")
fbvec2 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g79_dir-pe1_dwi.bvec")
fbvec3 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g81_dir-pe0_dwi.bvec")
fbvec4 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g81_dir-pe1_dwi.bvec")
bvec1 <- as.matrix(read.table(fbvec1))
bvec2 <- as.matrix(read.table(fbvec2))
bvec3 <- as.matrix(read.table(fbvec3))
bvec4 <- as.matrix(read.table(fbvec4))
bvec1[, 3]


## ----label="collecting-B0", echo=TRUE, eval=FALSE, results=FALSE--------------------------------------------
ddim <- c(130, 130, 81)
b0select <- c(2, 17, 34, 49, 66, 81)
ds1b0 <- array(0, c(ddim, 12))
ds1b0[, , , 1:6] <- ds1[, , , b0select]
ds1b0[, , , 7:12] <- ds2[, , , b0select]
ds1b0 <- as.nifti(ds1b0, ds1)
tmpfn <- file.path(tmpdir, "sub-01_ses-2015_topupin")
writeNIfTI(ds1b0, tmpfn)


## ----label="readingAcquisitionData", echo=TRUE, eval=TRUE, results=TRUE-------------------------------------
jsonFile1 <- file.path(dwid,
                       "sub-01_ses-2015_acq-g79_dir-pe1_dwi.json")
jsonFile2 <- file.path(dwid,
                       "sub-01_ses-2015_acq-g79_dir-pe1_dwi.json")
ttt1 <- read_json(jsonFile1)
ttt2 <- read_json(jsonFile2)
ttt1$InPlanePhaseEncodingDirection
ttt1$PhaseEncodingPolarityGE
ttt1$PhaseEncodingDirection 
ttt1$EchoTime
ttt2$InPlanePhaseEncodingDirection
ttt2$PhaseEncodingPolarityGE
ttt2$PhaseEncodingDirection 
ttt2$EchoTime


## ----label="creatingAcquisitionDataFile", echo=TRUE, eval=TRUE, results=TRUE--------------------------------
acqParFile <- file.path(tmpdir, "datainb0.txt")
dataIn <- matrix(0, 4, 12)
dataIn[, 1:6] <- c(0, 1, 0, ttt1$EchoTime)
dataIn[, 7:12] <- c(0, -1, 0, ttt2$EchoTime)
write(dataIn, acqParFile, ncolumns = 4)
dataIn


## ----label="creatingAcquisitionDataFile2", echo=TRUE, eval=TRUE, results=TRUE-------------------------------
acqParFileAll <- file.path(rdwipd, "datain.txt")
dataIn <- matrix(0, 4, 2)
dataIn[, 1] <- c(0, 1, 0, ttt1$EchoTime)
dataIn[, 2] <- c(0, -1, 0, ttt2$EchoTime)
write(dataIn, acqParFileAll, ncolumns = 4)

##  until now no use of fslr or ANTsR
## now use FSL/topup result in fntopup 
## ----label="Estimate_Distortions", echo=TRUE, eval=FALSE, message = FALSE, results=FALSE--------------------
fntopup <- file.path(rdwipd, "sub-01_ses-2015_topupout")
topio <- file.path(rdwipd, "sub-01_ses-2015_topupiout")
topup(infile = tmpfn,
      acqParFile,
      out = fntopup,
      subsamp = "1,1,1",
      fwhm = "8,4,0",
      miter = "4,4,8",
      iout = topio)


## ----label="generate brain mask", echo=TRUE, eval=FALSE, results=FALSE--------------------------------------
fnmask <- file.path(rdwipd, "sub-01_ses-2015_brain")
fslbet(topio, fnmask, opts = "-m")




## ----label="prepare for Eddy current correction", echo=TRUE, eval=FALSE, results=FALSE----------------------
dsall <- array(0, c(ddim, 328))
dsall[, , , 1:81] <- ds1
dsall[, , , 82:164] <- ds3
dsall[, , , 165:245] <- ds2
dsall[, , , 246:328] <- ds4
fndsall <- file.path(rdwipd, "sub-01_ses-2015_dwiall")
writeNIfTI(as.nifti(dsall, ds1), fndsall)
bval <- c(bval1, bval3, bval2, bval4)
bvec <- cbind(bvec1, bvec3, bvec2, bvec4)
fbval <- file.path(rdwipd, "sub-01_ses-2015_dwiall.bval")
fbvec <- file.path(rdwipd, "sub-01_ses-2015_dwiall.bvec")
write(bval, fbval, ncolumns = length(bval))
write.table(bvec, fbvec,
            row.names = FALSE, col.names = FALSE)
ind <- c(rep(1, 164), rep(2, 164))
fnindex <- file.path(rdwipd, "sub-01_ses-2015_dwi.ind")
write(ind, fnindex, ncolumns = length(ind))
rm(ds1,ds2,ds3,ds4,dsall)
gc()


## ----label="Eddy current correction", echo=TRUE, eval=FALSE, results=FALSE----------------------------------
fnout <- file.path(rdwipd, "sub-01_ses-2015_dwi_proc")
eddy(fndsall,
     fnmask,
     acqParFileAll,
     fnindex,
     fbvec,
     fbval,
     fntopup,
     fnout,
     eddy_cmd = "eddy_openmp")
# eddy_cmd = "eddy")


## ----label = "write bvals", echo = TRUE, eval = TRUE, message = FALSE, results = 'hide'---------------------
bval <- c(bval1, bval3, bval1, bval3)
fbval <- file.path(rdwipd,
                   "sub-01_ses-2015_dwi_proc.bval")
write(bval, fbval, ncolumns = length(bval))

