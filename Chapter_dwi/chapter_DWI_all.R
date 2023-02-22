## ----label="initQMRI", echo=FALSE, eval=TRUE-----------------------------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(warning = FALSE,
                      message = FALSE)
options(width=50)


## ----label="packages and filenames", eval=TRUE, echo=-1------------------------------------------------------------------------------------------
baseDir <- "../.."
codeDirDWI <- file.path(baseDir,"MRIwithR","Chapter_dwi")
source(file.path(codeDirDWI,"chapter_DWI_init.R"))



## ----label="initMPM", echo=FALSE, eval=TRUE, message=FALSE, results=FALSE------------------------------------------------------------------------
rimage.options(zquantiles=c(0.001,0.999), 
               xlab="x", ylab="z", bty="n")
setCores(4) # number of cores for parallel computations


## ----"determine code to be executed, no need to reproduce results already computed in code snippets", eval=TRUE, echo=FALSE----------------------
preprocessed <- file.exists(file.path(rdwipd,"sub-01_ses-106_dwi_proc.nii.gz"))
# have results from chapter_DWI_preprocess.R
if(!haveFSLr&&!preprocessed) stop("You need to install fslr or directory results")
havedataobj <- file.exists(file.path(rdwi,"dataobj.rsc"))
# data will be saved after first run of the document
havesigma <- file.exists(file.path(rdwipd, "sub-01_ses-106_sigma.nii.gz"))
# have results from chapter_DWI_varest.R
havetensor <- file.exists(file.path(rdwi,"tensorObjects.rsc"))
# have results from chapter_DWI_tensorest.R
havedkitensor <- file.exists(file.path(rdwi,"dkiObjects.rsc"))
# have results from chapter_DWI_dkiest.R
haveQball <- file.exists(file.path(rdwi,"qBallObjects.rsc"))
# have results from chapter_DWI_qBallest.R
havemixtensor <- file.exists(file.path(rdwi,"dmtcomb.rsc"))
# have results from chapter_DWI_mixtens.R
havesmooth <- file.exists(file.path(rdwi,"dtiobjsmooth.rsc"))
# have results from chapter_DWI_smooth.R
havexfibres <- file.exists(file.path(ldwipd,"dyads2.nii.gz"))
# have results from chapter_DWI_xfibres.R
haveprobtrackx <- file.exists(file.path(ptdwipd,"fdt_paths.nii.gz"))
# have results from chapter_DWI_probtrackx.R
havefibers <- file.exists(file.path(rdwi,"fibertracks.rsc"))
# have results from chapter_DWI_fibertrack.R
haveconnect <- file.exists(file.path(rdwi,"Connectivity.rsc"))
if(!haveANTsR&&!haveconnect) stop("You need to install ANTsR or directory results")
# have results from chapter_DWI_ANTs_atlasreg.R


## ----label="Figure_5_1", echo=-1, eval = TRUE, results=FALSE, message=FALSE, fig.cap="Effect of phase encoding direction on a non-diffusion weighted image. Left: A-P encoding direction, right: P-A encoding direction. The distortions are rather large in, but not restricted to, the anterior brain regions (upper image part). The both images can be combined to a single distortions-corrected one.", fig.width=6, fig.height=3.3, out.width='66%', fig.align='center'----
par(mfrow = c(1, 2),
    mar = c(1, 1, 2, 0.1), mgp = c(2, 1, 0))
run1FileName <- "sub-01_ses-106_run-001_dwi.nii.gz"
run2FileName <- "sub-01_ses-106_run-002_dwi.nii.gz"
run1File <- file.path(dwid, run1FileName)
run2File <- file.path(dwid, run2FileName)
ds1 <- readNIfTI(run1File)
ds2 <- readNIfTI(run2File)
rimage(ds1[ , , 41, 2], main = "A-P phase encoding")
rimage(ds2[ , , 41, 2], main = "P-A phase encoding")


## ----label="initialize dirs", echo=TRUE, eval = TRUE, results=FALSE, message=FALSE---------------------------------------------------------------
rdwipd <- file.path(resDir, "MyConnectome", "sub-01",
                   "ses-106", "dwi-proc")
if (!dir.exists(rdwipd)) dir.create(rdwipd)

## ----label="Reading unprocessed data", echo=TRUE, eval=!preprocessed, results=FALSE--------------------------------------------------------------
if(!preprocessed){
 run3File <- file.path(dwid,
                   "sub-01_ses-106_run-003_dwi.nii.gz")
 run4File <- file.path(dwid,
                   "sub-01_ses-106_run-004_dwi.nii.gz")
 ds3 <- readNIfTI(run3File)
 ds4 <- readNIfTI(run4File)
}

## ----label="Reading bvals", echo=TRUE, eval = TRUE, results=TRUE---------------------------------------------------------------------------------
fbval1 <- file.path(dwid,
                    "sub-01_ses-106_run-001_dwi.bval")
fbval2 <- file.path(dwid,
                    "sub-01_ses-106_run-002_dwi.bval")
fbval3 <- file.path(dwid,
                    "sub-01_ses-106_run-003_dwi.bval")
fbval4 <- file.path(dwid,
                    "sub-01_ses-106_run-004_dwi.bval")
bval1 <- scan(fbval1)
bval2 <- scan(fbval2)
bval3 <- scan(fbval3)
bval4 <- scan(fbval4)
bval1

## ----label="Reading bvecs", echo=TRUE, eval = TRUE, results=TRUE---------------------------------------------------------------------------------
fbvec1 <- file.path(dwid,
                    "sub-01_ses-106_run-001_dwi.bvec")
fbvec2 <- file.path(dwid,
                    "sub-01_ses-106_run-002_dwi.bvec")
fbvec3 <- file.path(dwid,
                    "sub-01_ses-106_run-003_dwi.bvec")
fbvec4 <- file.path(dwid,
                    "sub-01_ses-106_run-004_dwi.bvec")
bvec1 <- as.matrix(read.table(fbvec1))
bvec2 <- as.matrix(read.table(fbvec2))
bvec3 <- as.matrix(read.table(fbvec3))
bvec4 <- as.matrix(read.table(fbvec4))
bvec1[, 3]


## ----label="collecting-B0", echo=TRUE, eval=!preprocessed, results=FALSE-------------------------------------------------------------------------
if(!preprocessed){
 ddim <- dim(ds1)[1:3]
 b0select <- c(2, 17, 34, 49, 66, 81)
 ds1b0 <- array(0, c(ddim, 12))
 ds1b0[, , , 1:6] <- ds1[, , , b0select]
 ds1b0[, , , 7:12] <- ds2[, , , b0select]
 ds1b0 <- as.nifti(ds1b0, ds1)
 tmpfn <- file.path(tmpdir, "sub-01_ses-106_topupin")
 writeNIfTI(ds1b0, tmpfn)
}

## ----label="readingAcquisitionData", echo=TRUE, eval=TRUE, results=TRUE--------------------------------------------------------------------------
jsonFile1 <- file.path(dwijson,
                      "sub-01_ses-106_run-001_dwi.json")
jsonFile2 <- file.path(dwijson,
                      "sub-01_ses-106_run-002_dwi.json")
ttt1 <- read_json(jsonFile1)
ttt2 <- read_json(jsonFile2)
ttt1$InPlanePhaseEncodingDirection
ttt1$GradientEncodingDirection
ttt1$EchoTime
ttt2$InPlanePhaseEncodingDirection
ttt2$GradientEncodingDirection
ttt2$EchoTime

## ----label="creatingAcquisitionDataFile", echo=TRUE, eval=!preprocessed, results=TRUE------------------------------------------------------------
if(!preprocessed){
 acqParFile <- file.path(dwijson, "datainb0.txt")
 dataIn <- matrix(0, 4, 12)
 dataIn[, 1:6] <- c(0, 1, 0, ttt1$EchoTime)
 dataIn[, 7:12] <- c(0, -1, 0, ttt2$EchoTime)
 write(dataIn, acqParFile, ncolumns = 4)
 dataIn

## ----label="creatingAcquisitionDataFile2", echo=TRUE, eval=!preprocessed, results=TRUE-----------------------------------------------------------
 acqParFileAll <- file.path(dwijson, "datain.txt")
 dataIn <- matrix(0, 4, 2)
 dataIn[, 1] <- c(0, 1, 0, ttt1$EchoTime)
 dataIn[, 2] <- c(0, -1, 0, ttt2$EchoTime)
 write(dataIn, acqParFileAll, ncolumns = 4)

## ----label="Estimate_Distortions", echo=TRUE, eval=!preprocessed, message = FALSE, results=FALSE-------------------------------------------------
 fntopup <- file.path(rdwipd, "sub-01_ses-106_topupout")
 topio <- file.path(rdwipd, "sub-01_ses-106_topupiout")
 topup(infile = tmpfn,
       acqParFile,
       out = fntopup,
       subsamp = "1,1,1",
       fwhm = "8,4,0",
       miter = "4,4,8",
       iout = topio)

## ----label="generate brain mask", echo=TRUE, eval=!preprocessed, results=FALSE-------------------------------------------------------------------
 fnmask <- file.path(rdwipd, "sub-01_ses-106_brain")
 rwid <- file.path(resDir, "MyConnectome", "sub-01",
                   "ses-106", "dwi")
 fslbet(topio, fnmask, opts = "-m")
} else {
  
## ----label="generate brain mask 2", echo=FALSE, eval=preprocessed, results=FALSE-----------------------------------------------------------------
fnmask <- file.path(rdwipd, "sub-01_ses-106_brain")
}

## ----label="Figure_5_2", echo=-1, eval = TRUE, results=FALSE, message=FALSE, fig.cap="Result of the correction for the susceptibility-induced artifacts and the derived brain map on a non-diffusion weighted image.", fig.width=8, fig.height=4.3, out.width='66%', fig.align='center'----
par(mfrow = c(1, 2),
    mar = c(1, 1, 2, 0.1), mgp = c(2, 1, 0))
topio <- file.path(rdwipd, "sub-01_ses-106_topupiout")
fnmask <- file.path(rdwipd, "sub-01_ses-106_brain")
dscorrected <- readNIfTI(topio)
mask <- readNIfTI(paste0(fnmask, "_mask"))
rimage(dscorrected[ , , 41, 1], main = "Corrected image")
rimage(mask[ , , 41], main = "Estimated brain mask")

## ----label="prepare for Eddy current correction", echo=-(18:19), eval=!preprocessed, results=FALSE-----------------------------------------------
if(!preprocessed){
 dsall <- array(0, c(ddim, 328))
 dsall[, , , 1:81] <- ds1
 dsall[, , , 82:164] <- ds3
 dsall[, , , 165:245] <- ds2
 dsall[, , , 246:328] <- ds4
 fndsall <- file.path(rdwipd, "sub-01_ses-106_dwiall")
 writeNIfTI(as.nifti(dsall, ds1), fndsall)
 bval <- c(bval1, bval3, bval2, bval4)
 bvec <- cbind(bvec1, bvec3, bvec2, bvec4)
 fbval <- file.path(rdwipd, "sub-01_ses-106_dwiall.bval")
 fbvec <- file.path(rdwipd, "sub-01_ses-106_dwiall.bvec")
 write(bval, fbval, ncolumns = length(bval))
 write.table(bvec, fbvec,
             row.names = FALSE, col.names = FALSE)
 ind <- c(rep(1, 164), rep(2, 164))
 fnindex <- file.path(rdwipd, "sub-01_ses-106_dwi.ind")
 write(ind, fnindex, ncolumns = length(ind))
 rm(ds1,ds2,ds3,ds4,dsall)
 gc()

## ----label="Eddy current correction", echo=TRUE, eval=!preprocessed, results=FALSE---------------------------------------------------------------
 fnout <- file.path(rdwipd, "sub-01_ses-106_dwi_proc")
 eddy(fndsall,
             fnmask,
             acqParFileAll,
             fnindex,
             fbvec,
             fbval,
             fnout,
             eddy_cmd = "eddy_openmp")

## ----label = "write bvals", echo = TRUE, eval = !preprocessed, message = FALSE, results = 'hide'-------------------------------------------------
 bval <- c(bval1, bval3, bval1, bval3)
 fbval <- file.path(rdwipd,
                    "sub-01_ses-106_dwi_proc.bval")
 write(bval, fbval, ncolumns = length(bval))
}

## ----label="dataobj", echo=FALSE, eval=havedataobj-----------------------------------------------------------------------------------------------
if(havedataobj){
 load(file.path(rdwi,"dataobj.rsc"))
}

## ----"set xind,yind", echo=FALSE, eval = TRUE---------------------------------------------------------------------------------------------
xind <- 24:108
yind <- 11:125
  
## ----label = "readDWIdata", echo = TRUE, eval = !havedataobj, message = FALSE, results = 'hide'--------------------------------------------------
if(!havedataobj){
 gradientFile <- file.path(rdwipd,
       "sub-01_ses-106_dwi_proc.eddy_rotated_bvecs")
 bvalueFile <- file.path(rdwipd,
       "sub-01_ses-106_dwi_proc.bval")
 dataFile <- file.path(rdwipd,
       "sub-01_ses-106_dwi_proc.nii.gz")
 gradient <- as.matrix(read.table(gradientFile))
 bvalue <- as.numeric(scan(bvalueFile))
 dwobj <- readDWIdata(gradient,
                      dataFile,
                      "NIFTI",
                      bvalue = bvalue,
                      xind = xind,
                      yind = yind)
 fnmask <- file.path(rdwipd,"sub-01_ses-106_brain_mask.nii.gz")
 dwobj <- setmask(dwobj, fnmask)

## ----label = "saveDWIdata", echo = TRUE, eval = !havedataobj, message = FALSE, results = 'hide'--------------------------------------------------
 save(dwobj, file=file.path(rdwi,"dataobj.rsc"))
}

## ----label = "showSummaryDWOBJ", echo = TRUE, eval = TRUE, message = FALSE-----------------------------------------------------------------------
summary(dwobj)

## ----label = "showDimensionOfObject", echo = TRUE------------------------------------------------------------------------------------------------
dim(dwobj@si)

## ----label = "Figure_5_3", echo = -2, eval = TRUE, results = 'markup', fig.cap = "An axial slice of the the first 6 volumes of the example diffusion weighted data. The first and second volume are non-diffusion weighted image. Images are scaled for maximum contrast. This demonstrates that the diffusion weighted images have a smaller signal-to-noise ratio.", fig.width = 6, fig.height = 5.7, out.width='75%',fig.align='center'----
rimage.options(swapy=TRUE)
par(mfrow = c(2, 3),
    mar = c(.5, .5, 2, 0.1), mgp = c(1.5, .5, 0))
for (vol in 1:6) {
  rimage(dwobj@si[, , 41, vol],
         main = paste("Volume", vol, ",",
                      "b-value", dwobj@bvalue[vol]))
}

## ----label = "Figure_5_4", echo = -(1:2), eval = TRUE, results = 'markup', fig.cap = "The same axial slice of three volumes acquired with b-values 0, 1500 and 3000 $s/mm^2$. The same intensity range is used for all images.", fig.width = 6, fig.height = 2.92, fig.align='center',out.width='75%'----
dti.options(swapy=TRUE)
par(mfrow = c(1, 3),
    mar=c(.5, .5, 2, 0.1), mgp = c(1.5, .5, 0),
    xaxt = "n", yaxt = "n")
for (i in 16:18) {
  img <- plot(dwobj, slice = 41, gradient = i)
  title(paste("b-value", dwobj@bvalue[i]))
}

## ----label = "writeSlice", echo = TRUE, eval = FALSE---------------------------------------------------------------------------------------------
# write.image(img, file="myslice.png")

## ----label = "Figure_5_5", echo = TRUE, eval = TRUE, results = 'hide', rgl = TRUE, fig.cap = "3D visualization of diffusion weighted data in the object \\code{dwobj} by the function \\code{show3d}. The data values are coded as distance from the center in the direction of the diffusion gradients. Color codes gradient directions.", fig.width = 6, fig.height = 3, fig.align='center'----
show3d(dwobj,
       xind = 41:44, yind = 55:56, zind = 41,
       what = "data",
       zoom = .5, scale= c(0.8, 0.8, 0.8),
       windowRect = c(0, 0, 1000, 400),
       userMatrix = rotationMatrix(0, 1, 0, 1))

## ----label = "Figure_5_6", echo = TRUE, eval = TRUE, results = 'markup', rgl = TRUE, fig.cap = "3D visualization of the regular polyeders \\code{icosa0}, \\code{icosa1}, \\code{icosa2}, \\code{icosa3} and \\code{icosa4} from the dataset \\code{polyeders} in package \\pkg{dti}.", fig.width = 10, fig.height = 1.88, fig.align='center'----
data(polyeders)
icosa <- list(icosa0, icosa1, icosa2, icosa3, icosa4)
open3d()
par3d(zoom = 0.15, windowRect = c(0, 0, 1300, 250))
for(i in 1:5){
  ni <- icosa[[i]]$ni * 3
  ind <- icosa[[i]]$indices
  triangles3d(icosa[[i]]$vertices[1, ind] + 3*(i-1),
                icosa[[i]]$vertices[2, ind],
                icosa[[i]]$vertices[3, ind],
                color = rep("red", ni))
}

## ----label = "Figure_5_7", echo = TRUE, eval = TRUE, results = 'hide',  message = FALSE, fig.cap = "Threshold definition using the function \\code{sdpar}.", fig.width = 8.7, fig.height = 4, fig.align='center'----
dwobj <- sdpar(dwobj, level = 1.75)

## ----label = "define mask", echo=TRUE, eval = TRUE, results='markup'-----------------------------------------------------------------------------
lmask <- getmask(dwobj)

## ----label = "Figure_5_8", echo = -1, eval = TRUE,  fig.cap = "Mask definitions, central axial slice. Left: Mean non-diffusion weighted image, center: mask defined using threshold and \\code{sdpar}, right: mask obtained from \\code{topup}.", fig.width = 6., fig.height = 2.9, fig.align='center', out.width='75%'----
par(mfrow=c(1,3),mar=c(1,1,2,.1),mgp=c(2,1,0))
rimage(lmask$s0[,,41],col=grey(0:255/255))
title("Mean image")
rimage(lmask$mask[,,41],col=grey(0:255/255))
title("Mask from thresholding ")
rimage(dwobj@mask[,,41],col=grey(0:255/255))
title("Mask from fslr::topup")

## ----"Figure_5_9",eval = TRUE,echo=FALSE,results='hide',fig.cap="Effect of low SNR for Rician distribution. Densities for $\\sigma=1$. Vertical dashed and dotted lines are located at theta and the corresponding expected value, the length of the horizontal lines indicate the gap between both values.", fig.width=8, fig.height=4----
Lhalf <- function(x,L=1){
   gamma(L+1/2)/gamma(1.5)/gamma(L)*dti:::hg1f1(-0.5,L, -x)
}
drician <- function(y,theta=0,sigma=1){
  y/sigma^2*exp(-(y^2+theta^2)/2/sigma^2)*besselI(y*theta/sigma^2,0)
}
Erician <- function(theta,sigma=1){
   sigma*sqrt(pi/2)*Lhalf(-theta^2/2/sigma^2)
}
par(mar=c(3,3,3,1),mgp=c(2,1,0))
x <- seq(0,5.5,.01)
plot(x,drician(x,0),type="l",ylab="Density value")
lines(x,drician(x,1),col=2)
lines(x,drician(x,2),col=3)
lines(x,drician(x,3),col=4)
lines(x,drician(x,4),col=5)
lines(c(0,0),c(0,.3),lty=2,col=1,lwd=2)
lines(rep(Erician(0),2),c(0,.3),lty=3,col=1)
lines(c(0,Erician(0)),c(.05,.05),col=1,lend=2,lwd=2)
lines(c(1,1),c(0,.3),lty=2,col=2,lwd=2)
lines(rep(Erician(1),2),c(0,.3),lty=3,col=2)
lines(c(1,Erician(1)),c(.1,.1),col=2,lend=2,lwd=2)
lines(c(2,2),c(0,.3),lty=2,col=3,lwd=2)
lines(rep(Erician(2),2),c(0,.3),lty=3,col=3)
lines(c(2,Erician(2)),c(.15,.15),col=3,lend=2,lwd=2)
lines(c(3,3),c(0,.3),lty=2,col=4,lwd=2)
lines(rep(Erician(3),2),c(0,.3),lty=3,col=4)
lines(c(3,Erician(3)),c(.2,.2),col=4,lend=2,lwd=2)
lines(c(4,4),c(0,.3),lty=2,col=5)
lines(rep(Erician(4),2),c(0,.3),lty=3,col=5)
lines(c(4,Erician(4)),c(.25,.25),col=5,lend=2,lwd=2)
legend(4.5,.6,paste0("theta=",0:4),col=1:5,lty=rep(1,5))
title("Rician distributions")

## ----label="have variances", echo=FALSE, eval=havesigma------------------------------------------------------------------------------------------
if(havesigma){
 fnsigma <- file.path(rdwipd, "sub-01_ses-106_sigma")
 load(file.path(rdwi,"sigmaest.rsc"))
 sigma <- as.array(readNIfTI(fnsigma))
} else {

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
}

## ----label="Figure_5_10",echo=-c(1,10),eval = TRUE,fig.cap="Local noise estimation using \\code{aflsigmc}.",fig.width=8,fig.height=2.2,fig.align='center'----
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
rm(afsigma1b0,afsigma2b0)

## ----"using awslsigmc", eval = !havesigma, echo=TRUE, results='hide', warning=FALSE, message=FALSE-----------------------------------------------
if(!havesigma){
 kstar <- 16
 sb1 <- awslsigmc(dwobj1@si[ , , , 1], kstar, ncoils = 1)
 sb2 <- awslsigmc(dwobj2@si[ , , , 1], kstar, ncoils = 1)

## ----"using applytopup", eval = !havesigma, echo=TRUE--------------------------------------------------------------------------------------------
 fntopup <- file.path(dwipd, "sub-01_ses-106_topupout")
 fdn1 <- file.path(dwipd, "tmpn1")
 fdn2 <- file.path(dwipd, "tmpn2")
 fnsigma <- file.path(dwipd, "sub-01_ses-106_sigma")
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
 gc()
}

## ----"Figure_5_11", eval = TRUE, echo= -c(1,15),fig.cap="Local estimates of sigma and theta for non-diffusion  weighted (A-P and P-A) images. Right: Combined estimate of sigma (top) and densities of the estimated sigma values (bottom).",fig.width=8,fig.height=5.9,out.width='75%',fig.align='center'----
par(mfrow = c(2, 3), mar = c(1.5, 1.5, 2, .1),
    mgp = c(1, .5, 0))
rimage(sb1$sigma[, , 41], main = "sigma bv=0 A-P")
rimage(sb2$sigma[, , 41], main = "sigma bv=0 P-A")
rimage(sigma[, , 41], main = "sigma bv=0 combined")
rimage(sb1$theta[, , 41], main = "theta A-P")
rimage(sb2$theta[, , 41], main = "theta P-A")
d0 <- density(sb1$sigma, to = .6)
d1 <- density(sb2$sigma)
d2 <- density(sigma)
ylim <- range(d0$y, d1$y, d2$y)
plot(d0, ylim = ylim, ylab = "", xlab = "",
     main = "density(sigma)")
lines(d1, col = 2)
lines(d2, col = 3)
legend(0.3, ylim[2], c("A-P", "P-A", "combined"),
       col = 1:3, lty = rep(1,3))
rm(sb1, sb2)

## ----label="Figure_5_12",eval = TRUE, echo=TRUE,fig.cap="ADC for selected voxel.",fig.width=8,fig.height=3,rgl=TRUE,results='hide'---------------
show3d(dwobj,
       what = "adc",
       xind = 31:33, yind = 53, zind = 41,
       scale = 1, zoom = .4)

## ----"have tensor",  echo=FALSE, eval=havetensor-------------------------------------------------------------------------------------------------
if(havetensor){
 load(file.path(rdwi,"tensorObjects.rsc"))
 load(file.path(rdwi,"tensorIndices.rsc"))
} else {
  
## ----"DTI model",eval = !havetensor,echo=TRUE,results='hide'-------------------------------------------------------------------------------------
 dtiobj <- dtiTensor(dwobj, method = "linear")
}

## ----"DTI model 1",eval = TRUE-------------------------------------------------------------------------------------------------------------------
dim(dtiobj@D)

## ----"DTI model 2",eval = TRUE-------------------------------------------------------------------------------------------------------------------
signif(dtiobj@D[, 56, 56, 41],2)

## ----"DTI model 3",eval = TRUE-------------------------------------------------------------------------------------------------------------------
ind <- c(1, 2, 3, 2, 4, 5, 3, 5, 6)
DTensor <- matrix(dtiobj@D[ind, 56, 56, 41], 3, 3)
adcDT <- diag(t(dtiobj@gradient) %*%
                DTensor %*% dtiobj@gradient)
adcDT[3:5]

## ----"nlDTI model",eval = !havetensor,echo=TRUE,results='hide'-----------------------------------------------------------------------------------
if(!havetensor){
 dtiobjnl <- dtiTensor(dwobj, method = "nonlinear")
}

## ----"nl DTI model 2",eval = TRUE----------------------------------------------------------------------------------------------------------------
signif(matrix(dtiobjnl@D[, 56, 56, 41][ind],
              c(3, 3)), 3)

## ----"qlDTI model",eval=!havetensor,echo=TRUE,results='hide'-------------------------------------------------------------------------------------
if(!havetensor){
 sigma <- sigma[xind, yind, ]
 ## sigma contains negative values introduced in topup
 thresh <- quantile(sigma[mask], .1)
 sigma[sigma < thresh] <- thresh
 dtiobjql <- dtiTensor(dwobj,
                       method = "quasi-likelihood",
                       sigma = sigma)
}

## ----"qlDTI model res",eval = TRUE, results='markup'---------------------------------------------------------------------------------------------
signif(matrix(dtiobjql@D[, 56, 56, 41][ind],
              c(3, 3)), 3)

## ----"qlDTI model summary",eval = TRUE, results='markup'-----------------------------------------------------------------------------------------
summary(dtiobjql)

## ----"meanDiffusivity",echo=TRUE, eval=FALSE-----------------------------------------------------------------------------------------------------
## md <- dwiMD(dwobj)

## ----"dtiIndices", echo=TRUE, eval = !havetensor, results='hide'---------------------------------------------------------------------------------
if(!havetensor){
 dtiind   <- dtiIndices(dtiobj)
 dtiindnl <- dtiIndices(dtiobjnl)
 dtiindql <- dtiIndices(dtiobjql)
}

## ----label="Figure_5_13",echo=-1,eval = TRUE, fig.cap="Left: FA map of an axial slice of the example dMRI data. Right:  Corresponding map of the geodesic anisotropy. The latter has been up-scaled in the value range for more image contrast. Thus, large GA values are saturated.", fig.width=5.3, fig.height=3.85,out.width='66%',fig.align='center'----
 par(mfrow = c(1, 2),
     mar = c(1, 1, 2, .1), mgp = c(2, 1, 0))
 rimage(dtiindnl@fa[, , 41], main = "FA map")
 rimage(dtiindnl@ga[, , 41], zlim = c(0, 1.5),
        main = "GA map")

## ----"fig:dti3dprep",echo=-c(2,4,6), eval=TRUE,results='hide'------------------------------------------------------------------------------------
show3d(dtiobjql,
       xind = 30:39, yind = 28:34, zind = 41,
       scale = .7, zoom = .6, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0, 1, 0, 0))
snapshot3d("figure/Figure_5_14a.png")
show3d(dtiobjql,
       xind = 37, yind = 32, zind = 41,
       subdivide = 4, zoom = .7, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700))
snapshot3d("figure/Figure_5_14b.png")
show3d(dtiindql,
       xind = 30:39, yind = 28:34, zind = 41,
       zoom = .6, lwd = 4, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0, 1, 0, 0))
snapshot3d("figure/Figure_5_14c.png")

## ----"Figure_5_15",echo=-c(1,14),eval = TRUE, fig.cap="Upper row: FA maps of an axial slice of the example dMRI data obtained using the three available methods for tensor estimation. Lower row:  Corresponding map of the geodesic anisotropy. The latter has been up-scaled in the value range for more image contrast. Thus, large GA values are saturated.",fig.width=6.6,fig.height=6.4, out.width='75%',fig.align='center'----
par(mfrow = c(2, 3),
    mar = c(.1, .1, 2, .1), mgp = c(2, 1, 0),
    xaxt = "n", yaxt = "n")
plot(dtiind, slice = 41)
title("FA linear model")
plot(dtiindnl, slice = 41)
title("FA non-linear model")
plot(dtiindql, slice = 41)
title("FA quasi-likelihood")
plot(dtiind, what = "ga", slice = 41)
title("GA linear model")
plot(dtiindnl, what = "ga", slice = 41)
title("GA non-linear model")
plot(dtiindql, what = "ga", slice = 41)
title("GA quasi-likelihood")
rm(dtiobj,dtiobjnl,dtiobjql,dtiind,dtiindnl)

## ----"generate muenster color scale",eval = TRUE,echo=FALSE,results='hide'-----------------------------------------------------------------------
showFAColorScale("figure/Figure_5_16.png")

## ----"Figure_5_17",eval = TRUE,echo=-1,fig.cap="Color coded FA using various available color schemes.",fig.width=8,fig.height=2.92---------------
par(mfrow = c(1, 4), mar = c(.5, .5, 2, .1),
    mgp = c(2, 1, 0))
for (i in c(1, 3, 4, 6)){
  plot(dtiindql, slice = 41, method = i,
       xaxt = "n", yaxt = "n", mar = par("mar"))
  title(paste0("method=", i))
}

## ----"havedkitensor",eval=havedkitensor,echo=FALSE-----------------------------------------------------------------------------------------------
if(havedkitensor){
 load(file.path(rdwi,"dkiObjects.rsc"))
} else {
  
## ----"Compute dkitensor",eval = !havedkitensor,echo=TRUE, results='hide'-------------------------------------------------------------------------
 dkiobj <- dkiTensor(dwobj, method = "CLLS-QP")
}

## ----"DTensor from dki",echo=TRUE,eval = TRUE----------------------------------------------------------------------------------------------------
dim(dkiobj@D)
dim(dkiobj@W)

## ----"show3d on dkiobj",echo=-c(2,4),eval=TRUE,results='hide'------------------------------------------------------------------------------------
show3d(dkiobj,
       xind = 30:39, yind = 28:34, zind = 41,
       what = "DT", zoom = .6, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0, 1, 0, 0))
snapshot3d("figure/Figure_5_18a.png")
show3d(dkiobj,
       xind = 30:39, yind = 28:34, zind = 41,
       what = "KT", zoom = .6, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0, 1, 0, 0))
snapshot3d("figure/Figure_5_18b.png")

## ----"Compute dkiIndices",eval = !havedkitensor, echo=TRUE, results='hide',warning=FALSE,message=FALSE-------------------------------------------
if(!havedkitensor){
 dkiind <- dkiIndices(dkiobj)
}

## ----"Figure_5_19",eval = TRUE,echo=-c(1,6),fig.cap="Indices defined for the DKI model. From left to right: mean diffusivity, fractional (tensor) anisotropy, apparent kurtosis, and fractional kurtosis.", fig.width=8,fig.height=2.9----
par(mfrow = c(1, 4),
    mar = c(.5, .5, 2, .1), mgp = c(2, 1, 0),
    xaxt = "n", yaxt = "n")
plot(dkiind, slice = 41, what = "md", main="MD")
plot(dkiind, slice = 41, what = "fa", main="FA")
plot(dkiind, slice = 41, what = "mk2",
     zlim = c(0, .025), main="Kapp")
plot(dkiind, slice = 41, what = "fak", main="FAK")
rm(dkiind,dkiobj)

## ----"arttensor",echo= -c(13,16),eval = TRUE,results='hide'--------------------------------------------------------------------------------------
data(polyeders)
D <- diag(c(1, .3, .3))[c(1, 2, 3, 5, 6, 9)]
grad <- cbind(c(0, 0, 0), icosa2$vertices)
btb <- dti:::create.designmatrix.dti(grad)
bv <- c(0, rep(1, 162))
si <- exp(-bv*D%*%btb)
si <- as.nifti(array(si, c(1, 1, 1, 163)))
fn <- file.path(tmpdir,"tensor")
writeNIfTI(si, fn)
dstens <- readDWIdata(grad,
                      file.path(tmpdir,"tensor.nii.gz"),
                      "NIFTI", bvalue = bv)
odf0 <- dwiQball(dstens,
                 what = "ODF",
                 mask = array(TRUE, c(1, 1, 1)),
                 order = 6)
show3d(odf0, bgcolor = "white")
snapshot3d("figure/Figure_5_20a.png")
odfw <- dwiQball(dstens,
                 what = "wODF",
                 mask = array(TRUE, c(1, 1, 1)),
                 order = 6)
show3d(odfw, bgcolor = "white")
snapshot3d("figure/Figure_5_20b.png")

## ----"haveQball",eval=haveQball,echo=FALSE-------------------------------------------------------------------------------------------------------
if(haveQball){
 load(file.path(rdwi,"qBallObjects.rsc"))
} else {

## ----"dwiQball1",eval = !haveQball, echo=TRUE, results='hide'------------------------------------------------------------------------------------
 qballw4 <- dwiQball(dwobj, what = "wODF",
                     order = 4, lambda = 1e-2)
}

## ----"Figure_5_21",eval = TRUE, echo=TRUE, results='hide',fig.cap="Qball-reconstruction order 4 in a region from a central slice.", fig.width=10, fig.height=7, rgl=TRUE----
show3d(qballw4,
       xind = 30:39, yind = 28:34, zind = 41,
       scale = .7, zoom = .6, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0, 1, 0, 0))


## ----"dwiQball2",eval = !haveQball, echo=FALSE, results='hide'-----------------------------------------------------------------------------------
if(haveQball){
 qballw4a <- dwiQball(dwobj, what="wODF",
                      order=4, lambda=1e-3)
 qballw4c <- dwiQball(dwobj, what="wODF",
                      order=4, lambda=1e-1)
 qballw6a <- dwiQball(dwobj, what="wODF",
                      order=6, lambda=1e-3)
 qballw6b <- dwiQball(dwobj, what="wODF",
                      order=6, lambda=1e-2)
 qballw6c <- dwiQball(dwobj, what="wODF",
                      order=6, lambda=1e-1)
 qballw8a <- dwiQball(dwobj, what="wODF",
                      order=8, lambda=1e-3)
 qballw8b <- dwiQball(dwobj, what="wODF",
                      order=8, lambda=1e-2)
 qballw8c <- dwiQball(dwobj, what="wODF",
                      order=8, lambda=1e-1)
}

## ----"Qball figures", eval=TRUE, echo=FALSE, results='hide'--------------------------------------------------------------------------------------
show3d(qballw4,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
snapshot3d("figure/Figure_5_22w4e2.png")
show3d(qballw4a,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
snapshot3d("figure/Figure_5_22w4e3.png")
show3d(qballw4c,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
snapshot3d("figure/Figure_5_22w4e1.png")
show3d(qballw6a,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
snapshot3d("figure/Figure_5_22w6e3.png")
show3d(qballw6b,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
snapshot3d("figure/Figure_5_22w6e2.png")
show3d(qballw6c,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
snapshot3d("figure/Figure_5_22w6e1.png")
show3d(qballw8a,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
snapshot3d("figure/Figure_5_22w8e3.png")
show3d(qballw8b,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
snapshot3d("figure/Figure_5_22w8e2.png")
show3d(qballw8c,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
snapshot3d("figure/Figure_5_22w8e1.png")

## ----"Overlay", eval=TRUE, echo=TRUE-------------------------------------------------------------------------------------------------------------
overlay <- function(fn1, fn2, thresh = 100){
  img1 <- read.image(fn1)
  img2 <- read.image(fn2)
  img1 <- shrink.image(img1,
                       xt = img2$dim[1],
                       yt = img2$dim[2],
                       ratio = FALSE)
  over <- function(x1, x2, thresh = 100){
    if(length(dim(x2)) == 2){
       x2[x2 < thresh] <- x1[x2 < thresh]
    } else {
       ind <- apply(x2 < thresh, c(1, 2), "all")
       dim(x2) <- c(prod(dim(x2)[1:2]), dim(x2)[3])
       for(i in 1:dim(x1)[3])
         x2[ind, i] <- x1[, , i][ind]
       dim(x2) <- dim(x1)
    }
    x2
  }
  combine(img1, img2, over, thresh = thresh)
}
img <- plot(dtiindql[30:39, 28:34, 41], show = FALSE)
write.image(img, file="figure/coloredFA.png")

## ----"Figure_5_23",eval = TRUE,echo=-2,fig.cap="ODF's overlayed on a tensor FA map.", fig.width=6, fig.height=4.2,results='hide'-----------------
show3d(qballw4,
       xind = 30:39, yind = 28:34, zind = 41,
       scale = 0.7, zoom = 0.59,
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0*pi, 1, 0, 0))
snapshot3d("figure/wODFoverlay.png")
par(mar = c(0, 0, 0, 0))
show.image(overlay("figure/coloredFA.png",
                   "figure/wODFoverlay.png", 100),
           xaxt = "n", yaxt = "n")

## ----"rm Qball",eval=TRUE,echo=FALSE,results=FALSE-----------------------------------------------------------------------------------------------
rm(qballw4,qballw4a,qballw4c,qballw6a,qballw6b,qballw6c,qballw8a,qballw8b,qballw8c)
gc()

## call to xfibres needs > 1 week
## ----"Using fsl/xfibres", eval=FALSE, echo=TRUE--------------------------------------------------------------------------------------------------
## gradientFile <- file.path(rdwipd,
##           "sub-01_ses-106_dwi_proc.eddy_rotated_bvecs")
## bvalueFile <- file.path(rdwipd,
##                         "sub-01_ses-106_dwi_proc.bval")
## dataFile <- file.path(rdwipd,
##                       "sub-01_ses-106_dwi_proc")
## fnmask <- file.path(rdwipd,"sub-01_ses-106_brain_mask")
## xfibres(dataFile, gradientFile, bvalueFile, fnmask, 1,
##         verbose=FALSE, opts="--model=2")

## ----label="initadimpro3", echo=FALSE, eval=TRUE, message=FALSE----------------------------------------------------------------------------------
rimage.options(zquantiles=c(0.001,0.98),xlab="x",ylab="z",bty="n",xaxt="n",yaxt="n")

## ----"Figure_5_24", eval = havexfibres, echo=-2, fig.cap="Results for the ball and (2)-sticks model: unweighted image, diffusivity, compartment size, six components of estimated directions.",fig.width=8,fig.height=4.7----
rimage.options(swapy = FALSE)
par(mfrow = c(2, 5), mar = c(.5, .5, 2, .1), 
    mgp = c(2, 1, 0), xaxt = "n", yaxt = "n")
imgdyad1 <- readNIfTI(file.path(ldwipd,"dyads1.nii.gz"))
imgdyad2 <- readNIfTI(file.path(ldwipd,"dyads2.nii.gz"))
imgf1 <- readNIfTI(file.path(ldwipd,
                             "mean_f1samples.nii.gz"))
imgf2 <- readNIfTI(file.path(ldwipd,
                             "mean_f2samples.nii.gz"))
imgS0 <- readNIfTI(file.path(ldwipd,
                             "mean_S0samples.nii.gz"))
imgd <- readNIfTI(file.path(ldwipd,
                            "mean_dsamples.nii.gz"))
rimage(as.array(imgS0)[xind, yind, 41], main = "S0")
rimage(as.array(imgf1)[xind, yind, 41], main = "f1")
rimage(as.array(imgdyad1)[xind, yind, 41, 1], 
       main = "x_11")
rimage(as.array(imgdyad1)[xind, yind, 41, 2],
       main = "x_12")
rimage(as.array(imgdyad1)[xind, yind, 41, 3],
       main = "x_13")
rimage(as.array(imgd)[xind, yind, 41], 
       main = "Diffusivity")
rimage(as.array(imgf2)[xind, yind, 41],
       main = "f2")
rimage(as.array(imgdyad2)[xind, yind, 41, 1],
       main = "x_21")
rimage(as.array(imgdyad2)[xind, yind, 41, 2],
       main = "x_22")
rimage(as.array(imgdyad2)[xind, yind, 41, 3],
       main = "x_23")

## ----"havemixtensor",eval=havemixtensor,echo=FALSE-----------------------------------------------------------------------------------------------
if(havemixtensor){
 load(file.path(rdwi,"dmtcomb.rsc"))
} else {
  
## ----"Tensor mixtures", eval=!havemixtensor, echo=-c(15,16), results='hide'----------------------------------------------------------------------
 dmtobj5 <- dwiMixtensor(dwobj, maxcomp = 5,
                         model = "MTiso")
 dmtobj4 <- dwiMixtensor(dwobj, maxcomp = 4,
                         model = "MTiso")
 dmtobj3 <- dwiMixtensor(dwobj, maxcomp = 3,
                         model = "MTiso")
 dmtobj2 <- dwiMixtensor(dwobj, maxcomp = 2,
                         model = "MTiso")
 dmtobj1 <- dwiMixtensor(dwobj, maxcomp = 1,
                         model = "MTiso")
 dmtcomb <- dwiMtCombine(dmtobj5, dmtobj4)
 dmtcomb <- dwiMtCombine(dmtcomb, dmtobj3)
 dmtcomb <- dwiMtCombine(dmtcomb, dmtobj2)
 dmtcomb <- dwiMtCombine(dmtcomb, dmtobj1)
 rm(dmtobj1,dmtobj2,dmtobj3,dmtobj4,dmtobj5)
 gc()
}

## ----"Figure_5_25", eval = TRUE, echo=TRUE, results='hide',fig.cap="Weighted ODF for the tensor mixture model described in the text for the voxel selection of Figure~\\ref{fig:Figure_5_23} estimated using  the dwiMixtensor function of package \\pkg{dti}.", fig.width=8,fig.height=5,rgl=TRUE----
show3d(dmtcomb,
       xind = 30:39, yind = 28:34, zind = 41,
       subdivide = 3, scale = .4, zoom = .6,
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0*pi, 1, 0, 0))

## ----"extract from mtobj",echo=TRUE,eval = !havemixtensor----------------------------------------------------------------------------------------
if(!havemixtensor){
 mtindices <- extract(dmtcomb, what =
                      c("w0", "fa", "eorder", "order"))
}

## ----"Figure_5_26", eval = TRUE, echo=-1, results='hide',fig.cap="Tensor mixture model of maximal order five: Isotropic compartment size, fractional anisotropy (FA), order of mixture and effective order (EO) of mixture (from left to right) for the central axial slice.",fig.width=10.,fig.height=3.6----
par(mar = c(.5, .5, 2, .1), mgp = c(2, 1, 0))
plot(dmtcomb, what = c("w0", "fa", "order", "eorder"),
     slice = 41, xaxt = "n", yaxt = "n")

## ----"havesmooth", eval=havesmooth, echo=FALSE---------------------------------------------------------------------------------------------------
if(havesmooth){
 load(file.path(rdwi,"dtiobjsmooth.rsc"))
} else {

## ----"Gaussian smoothing of dwi data",echo=TRUE,eval = !havesmooth,results='hide'----------------------------------------------------------------
 dwobj.gauss <- dwobj
 for(i in 1:dwobj@ngrad){
   dwobj.gauss@si[, , , i] <-
         kernsm(dwobj@si[, , , i], h = 1)@yhat
 }

## ----"qlDTI model gauss",eval = !havesmooth, echo=-5, results='hide'-----------------------------------------------------------------------------
 dtiobjql.gauss <- dtiTensor(dwobj.gauss,
                          method = "quasi-likelihood",
                          sigma = sigma)
 dtiindql.gauss <- dtiIndices(dtiobjql.gauss)
 rm(dwobj.gauss, dtiobjql.gauss)

## ----"smooth dtiTensor object", eval=!havesmooth, echo=-3----------------------------------------------------------------------------------------
 dtiobj.sm <- dti.smooth(dwobj, hmax=3)
 dtiind.sm <- dtiIndices(dtiobj.sm)
 rm(dtiobj.sm,dtiind.sm)

## ----"POAS multishell",eval=!havesmooth,echo=-5--------------------------------------------------------------------------------------------------
 sigmap <- awslsigmc(dwobj@si[ , , , 1],
                     steps = 16)$sigma
 dwobj.poas <- dwi.smooth.ms(dwobj, kstar = 12,
                             sigma = sigmap)
 rm(sigmap)

## ----"estimate diffusion tensor model QL",echo=-c(5:6),eval=!havesmooth--------------------------------------------------------------------------
 dtiobjql.poas <- dtiTensor(dwobj.poas,
                            method = "quasi-likelihood",
                            sigma = sigma)
 dtiindql.poas <- dtiIndices(dtiobjql.poas)
 rm(dwobj.poas,dtiobjql.poas)
 gc()
}

## ----"Figure_5_27",eval = TRUE,echo=-c(1,5),fig.cap="Color coded FA maps for tensor estimates using the original data (left), data smoothed by POAS (center) and Gaussian smoothing (right).",fig.width=8.,fig.height=3.83, out.width='80%',fig.align='center'----
par(mfrow = c(1, 3),
    mar = c(.5, .5, 2, .1), mgp = c(2, 1, 0))
plot(dtiindql, slice = 41, xaxt = "n", yaxt = "n")
plot(dtiindql.poas, slice = 41, xaxt = "n", yaxt = "n")
plot(dtiindql.gauss, slice = 41, xaxt = "n", yaxt = "n")
rm(sigma,dtiindql.gauss)

## ----"havefibers", echo=FALSE, eval=havefibers-------------------------------------------------------------------
if(havefibers){
  load(file.path(rdwi,"fibertracks.rsc"))
} else {
  
## ----"Fiber tracking from tensors",eval=!havefibers,echo=TRUE---------------------------------------------------
  trxql <- tracking(dtiindql)
  trxql <- selectFibers(trxql, minlength = 20)
  trxql <- reduceFibers(trxql)
  
## ----"Fiber tracking from tensors smooth",eval=!havefibers,echo=TRUE------------------------------------------
  trxql.poas <- tracking(dtiindql.poas)
  trxql.poas <- selectFibers(trxql.poas, minlength = 20)
  trxql.poas <- reduceFibers(trxql.poas)
  
## ----"Fiber tracking from tensor rmixtures",eval=!havefibers,echo=TRUE----------------------------------------
  trxcomb5 <- tracking(dmtcomb, mincompartsize = .05)
  trxcomb5 <- trxcomb50 <- reduceFibers(trxcomb5)
  trxcomb5 <- selectFibers(trxcomb5, minlength = 20)
}
  
## ----"rm dmtcomb", eval=TRUE, echo=FALSE---------------------------------------------------------------------------------------------------------
rm(dmtcomb)

## ----"Tractography illustrations",eval = TRUE,echo=c(1:2,4,6)------------------------------------------------------------------------------------
wRec <- c(0, 0, 1000, 1000)
show3d(trxql,
       zoom = .5, windowRect = wRec, bgcolor = "white")
snapshot3d("figure/Figure_5_29a.png")
show3d(trxql.poas,
       zoom = .5, windowRect = wRec, bgcolor = "white")
snapshot3d("figure/Figure_5_29b.png")
show3d(trxcomb5,
       zoom = .5, windowRect = wRec, bgcolor = "white")
snapshot3d("figure/Figure_5_29c.png")

## probabilistic fibertracking needs > 1 week
## ----"seedmask", eval=FALSE, echo=TRUE-----------------------------------------------------------------------------------------------------------
## fast(file.path(rdwipd,"sub-01_ses-106_brain.nii.gz"),
##      outfile=rdwipd)
## img2 <- readNIfTI(file.path(rdwipd,"_pve_2.nii.gz"))
## img1 <- readNIfTI(file.path(rdwipd,"_pve_1.nii.gz"))
## writeNIfTI(img1*img2>0,file.path(rdwipd,"seeds"))

## ----"Probablilistic fiber tracking",eval=FALSE,echo=TRUE----------------------------------------------------------------------------------------
## sdwipd <- file.path(ldwipd,"merged")
## fnmask <- file.path(rdwipd,
##                     "sub-01_ses-106_brain.nii.gz")
## seeds <- file.path(rdwipd,"seeds.nii.gz")
## opt <- paste("--forcedir -V 2 -l --onewaycondition -c
##    0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01
##    --distthresh=25.0 --sampvox=0.0 --opd --omatrix1
##    --dir=",ptdwipd)
## probrackx(samples = sdwipd, mask = fnmask,
##           seed = fnmask, opt=opt)

## ----"probtrackx call",echo=TRUE,eval=FALSE------------------------------------------------------------------------------------------------------
## cmd <- paste("probtrackx2 -s",sdwipd,"-m",fnmask,"-x",seeds,opt)
## system(cmd)

## ----"Figure_5_30", eval=haveprobtrackx, echo=TRUE,out.width='75%', fig.align="center",fig.cap="Orthographic view of fiber density estimated by the function \\code{probtrackx}.",fig.width=10,fig.height=8--------------------------------------
if(haveprobtrackx){ 
 fdtpaths <- readNIfTI(file.path(ptdwipd,
                                  "fdt_paths.nii.gz"))
 orthographic(fdtpaths, xyz = c(65.5, 65.5, 46),
             zlim = quantile(fdtpaths, c(0, .99)))
}
  
## ----"haveconnect", eval=haveconnect, echo=FALSE--------------------------------------
if(haveconnect){
  load(file.path(rdwi,"Connectivity.rsc"))
  fnHOatlas <- file.path(rdwipd,
                       "sub-01_ses-106_HOatlas.nii.gz")
} else {
  
## ----"Produce a mask for HarvardOxford atlas",eval=!haveconnect,echo=TRUE--------------------------------------
  pho <- readNIfTI(file.path(fsl_atlas_dir(),
                             "HarvardOxford",
                             "HarvardOxford-cort-prob-1mm.nii.gz"))
 img <- as.nifti(apply(pho, 1:3, sum) > 0, pho)
 fnhomask <- file.path(tmpdir, "HOmask")
 writeNIfTI(img, fnhomask)
  
## ----"Register HarvardOxford atlas to DWI subject space compute transform",eval=!haveconnect,echo=TRUE--------------------------------------
  fnmask <- file.path(rdwipd,
                      "sub-01_ses-106_brain_mask.nii.gz")
 mask <- antsImageRead(fnmask)
 homask <- antsImageRead(file.path(tmpdir,
                                  "HOmask.nii.gz"))
 ho2dwi <- antsRegistration(fixed = mask,
                           moving = homask,
                           typeofTransform = "Affine")

## ----"Register HarvardOxford atlas to DWI subject space apply transform",eval=!haveconnect,echo=TRUE--------------------------------------
  fnatlas <- file.path(fsl_atlas_dir(), "HarvardOxford",
                       "HarvardOxford-cort-maxprob-thr25-1mm.nii.gz")
 atlas <- antsImageRead(fnatlas)
 hoindwi <- antsApplyTransforms(mask, atlas,
                               ho2dwi$fwdtransforms, interpolator = "genericLabel")
 fnHOatlas <- file.path(rdwipd,
                       "sub-01_ses-106_HOatlas.nii.gz")
 antsImageWrite(hoindwi, fnHOatlas)
  
## ----"Calculate adjacency matrices",eval=!haveconnect,echo=TRUE,warning=FALSE,message=FALSE--------------------------------------
  HOatlas <- readNIfTI(fnHOatlas)
 zMT <- AdjacencyMatrix(trxcomb50, HOatlas)
}

## ----"rm tracks", eval=TRUE, echo=FALSE--------------------------------------
  rm(trxql, trxql.poas, trxcomb5)

## ----"Figure_5_31", eval=TRUE, echo=TRUE,out.width='100%', fig.align="center",fig.cap="Anatomic connectivity for the HarvardOxford cortical atlas: Adjacency matrices obtained from fiber tracking results for the tensor mixture model.",fig.width=10,fig.height=10----
heatmap(zMT, symm = TRUE, col = grey(0:255/255))
title("Tensor Mixture Model")

## ----"construct connectivity graphs", eval=!haveconnect, echo=TRUE--------------------------------------
if(!haveconnect){
  set.seed(1)
  gMT <- graph_from_adjacency_matrix(zMT,
                                   "undirected",
                                   weighted = TRUE)
}
  
## ----"Figure_5_32", eval=TRUE, echo=-1,out.width='100%', fig.align="center",fig.cap="Anatomic connectivity for the HarvardOxford cortical atlas: Visualization of connectivity graphs obtained from fiber tracking results for the tensor mixture model.",fig.width=11,fig.height=6----
par(mfrow = c(1, 2),
    mar = c(3, 3, 3, 1), mgp = c(2, 1, 0))
plot(permute(gMT,canonical_permutation(gMT)$labeling),
     layout = layout_on_sphere,
     edge.color = grey(33:100/100))
plot(permute(gMT,canonical_permutation(gMT)$labeling),
     layout = layout_in_circle,
     edge.color = grey(33:100/100))
title("Connectivity graph from MT")

## ----"read atlas metadata",echo=TRUE,eval=!haveconnect,warning=FALSE,meessage=FALSE--------------------------------------
if(!haveconnect){
  HOmeta <- xmlTreeParse(file.path(fsl_atlas_dir(),
                                   "HarvardOxford-Cortical.xml"))
 HOmeta <- xmlToList(HOmeta)
 regionNames <- as.character(HOmeta$data[1, ])
 coords <- matrix(as.numeric(unlist(HOmeta$data[2, ])),
                 48, 4)[, 2:4]
}
  