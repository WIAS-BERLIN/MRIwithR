## ----label="initQMRI", echo=FALSE, eval=TRUE-----------------------------------------------------------------------------------------------------

## ----label="packages and filenames", eval=TRUE, echo=-1------------------------------------------------------------------------------------------
if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_dwi","chapter_DWI_init.R"))

run1FileName <- "sub-01_ses-2015_acq-g79_dir-pe0_dwi.nii.gz"
run2FileName <- "sub-01_ses-2015_acq-g79_dir-pe1_dwi.nii.gz"
run1File <- file.path(dwid, run1FileName)
run2File <- file.path(dwid, run2FileName)


## ----label="Reading bvals", echo=TRUE, eval = TRUE, results=TRUE---------------------------------------------------------------------------------
fbval1 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g79_dir-pe0_dwi.bval")
fbval2 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g79_dir-pe1_dwi.bval")
bval1 <- scan(fbval1)
bval2 <- scan(fbval2)

## ----label="Reading bvecs", echo=TRUE, eval = TRUE, results=TRUE---------------------------------------------------------------------------------
fbvec1 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g79_dir-pe0_dwi.bvec")
fbvec2 <- file.path(dwid,
                    "sub-01_ses-2015_acq-g79_dir-pe1_dwi.bvec")
bvec1 <- as.matrix(read.table(fbvec1))
bvec2 <- as.matrix(read.table(fbvec2))

xind <- 24:108
yind <- 11:125
  

## ----label="reading unprocessed data1", echo=TRUE, eval = !havesigma, results='hide'-------------------------------------------------------------
 dwobj1 <- readDWIdata(bvec1,
                       run1File,
                       bvalue = bval1,
                       format = "NIFTI")

## ----label="using afsigmc on unprocessed data", echo=TRUE, eval = !havesigma, results='markup'---------------------------------------------------
 v1 <- afsigmc(dwobj1@si[ , , , 1], level = 1.75,
               ncoils=1, method="modevn")
 v2 <- afsigmc(dwobj1@si[ , , , 1], level = 1,
               ncoils=1, method="modem1chi")
 v3 <- afsigmc(dwobj1@si[ , , , 1], level = 1,
               ncoils=1, method="bkm2chi")
 v4 <- afsigmc(dwobj1@si[ , , , 1], level = 1,
               ncoils=1, method="bkm1chi")
 cat("Estimated sigma: ", signif(c(v1, v2, v3, v4), 3))

## ----label="reading unprocessed data2", echo=TRUE, eval = !havesigma, results='hide'-------------------------------------------------------------
 dwobj2 <- readDWIdata(bvec2,
                       run2File,
                       bvalue = bval2,
                       format = "NIFTI")

## ----label="AF2013",echo=TRUE,eval = !havesigma, results='hide'----------------------------------------------------------------------------------
 afsigma1b0 <- aflsigmc(dwobj1@si[, , , 1],
                        ncoils = 1, level = 1)
 afsigma2b0 <- aflsigmc(dwobj2@si[, , , 1],
                        ncoils = 1, level = 1)



## ----"using awslsigmc", eval = !havesigma, echo=TRUE, results='hide', warning=FALSE, message=FALSE-----------------------------------------------
 kstar <- 16
 sb1 <- awslsigmc(dwobj1@si[ , , , 1], kstar, ncoils = 1)
 sb2 <- awslsigmc(dwobj2@si[ , , , 1], kstar, ncoils = 1)

## ----"using applytopup", eval = !havesigma, echo=TRUE--------------------------------------------------------------------------------------------
 acqParFileAll <- file.path(rdwipd, "datain.txt")
 fntopup <- file.path(rdwipd, "sub-01_ses-2015_topupout")
 fdn1 <- file.path(rdwipd, "tmpn1")
 fdn2 <- file.path(rdwipd, "tmpn2")
 fnsigma <- file.path(rdwipd, "sub-01_ses-2015_sigma")
 ds1 <- readNIfTI(run1File)
 writeNIfTI(as.nifti(sb1$sigma, ds1), fdn1)
 ds2 <- readNIfTI(run2File)
 writeNIfTI(as.nifti(sb1$sigma, ds2), fdn2)
 applytopup(c(fdn1, fdn2),
            acqParFileAll,
            c(1, 2),
            fntopup,
            fnsigma)
 sigma <- as.array(readNIfTI(fnsigma))

## ----label="cleanup noise", echo=FALSE, eval = !havesigma, results='hide'------------------------------------------------------------------------

  rm(dwobj1, dwobj2, ds1, ds2)

save(sb1, sb2, sigma, afsigma1b0, afsigma2b0, v1, v2, v3, v4, 
     file=file.path(rdwi,"sigmaest.rsc"))

  
