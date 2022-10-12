## ----label="init", echo=FALSE, eval=TRUE--------------------------------------------------------------------
library(knitr)
library(rgl)
knitr::opts_chunk$set(cache = FALSE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE)
knit_hooks$set(rgl = hook_rgl)
options(width=50)


## ----label="packages", eval=TRUE, echo=-1-------------------------------------------------------------------
library(knitr)
library(adimpro)
library(ANTsR)
library(aws)
library(fmri)
library(fslr)
library(glasso)
library(igraph)
library(mritc)
library(neuRosim)
library(oro.nifti)
library(rgl)
library(XML)


## ----label="initadimpro3", echo=FALSE, eval=TRUE, message=FALSE---------------------------------------------
rimage.options(zquantiles=c(0.001,0.98),xlab="x",ylab="z",bty="n",xaxt="n",yaxt="n")


## ----label="filename declarations", eval=TRUE---------------------------------------------------------------
dataDir <- "../data"
resDir <- "../results"
f117s1r1 <- file.path(dataDir, "ds000117", "sub-01",
 "ses-mri", "func",
 "sub-01_ses-mri_task-facerecognition_run-01_bold.nii")
f117s1r1stc <- file.path(resDir, "ds000117", "sub-01",
    "ses-mri", "func", paste0("stc",
    "sub-01_ses-mri_task-facerecognition_run-01_bold"))
f <- "sub-1_task-objectviewing_run-01_bold.nii.gz"
f105df <- file.path(dataDir, "ds000105",
                      "sub-1", "func")
f105dfwr <- file.path(resDir, "ds000105",
                      "sub-1", "func")
f105s1r1 <- file.path(dataDir, f105df, f)
f105s1r1mc <- file.path(f105dfwr, paste0("mc", f))
f105s1r1amc <- file.path(f105dfwr, paste0("amc", f))
f105s1r1motion <- file.path(f105df, "mocoparams.rsc")
f105s1r1 <- file.path(dataDir, "ds000105",
                      "sub-1", "func", f)
f105s1r1mc <- file.path(resDir, "ds000105",
            "sub-1", "func", paste0("mc",f))
f105s1r1amc <- file.path(resDir, "ds000105",
            "sub-1", "func", paste0("amc",f))
f105s1r1motion <- file.path(resDir, "ds000105",
            "sub-1", "func", "mocoparams.rsc")
f105s1r1p <- file.path(dataDir, "ds000105",
            "sub-1", "func", paste0("p",f))
f105s1T1 <- file.path(resDir, "ds000105",
            "sub-1", "anat", "sub-1_T1w.nii.gz")
f105s1bold2T1 <- file.path(resDir, "ds000105",
            "sub-1", "anat", "sub-1_bold2T1w.nii.gz")
fileMNI <- file.path(dataDir, "mni",
                     "icbm_avg_152_t1_tal_lin.nii")
f105s1r1nmc <- file.path(resDir, "ds000105",
                     "sub-1","func", paste0("nmc", f))
fileshen268 <- file.path(dataDir, "atlas",
                    "shen_1mm_268_parcellation.nii.gz")
fileshen268mni <- file.path(resDir, "atlas",
                    "shen_mni_268_parcellation.nii.gz")
f105s1r1atlas268 <- file.path(resDir, "ds000105",
                   "sub-1", "func", paste0("l268", f))
f105s1boldmask <- file.path(resDir, "ds000105",
              "sub-1", "anat", "sub-1_boldmask.nii.gz")
f105s1brain <- file.path(resDir, "ds000105",
                "sub-1", "anat", "sub-1_brain.nii.gz")
f105s1mask <- file.path(resDir, "ds000105",
            "sub-1", "anat", "sub-1_brain_mask.nii.gz")
f105s1brainBET <- file.path(resDir, "ds000105",
                       "sub-1", "anat", "fslbet_brain")
f105s1brainFAST <- file.path(resDir, "ds000105",
                      "sub-1", "anat", "fslfast_brain")
f1ev <- "sub-1_task-objectviewing_run-01_events.tsv"
f105s1r1events <- file.path(dataDir, "ds000105",
                             "sub-1", "func", f1ev)
g <- "sub-01_ses-105_task-rest_run-010_bold.nii.gz"
frest <- file.path(dataDir, "MyConnectome",
                       "sub-01", "ses-105", "func", g)
fprest <- file.path(resDir, "MyConnectome", "sub-01",
                 "ses-105", "func", paste0("prmc", g))
fprestmask <- file.path(resDir, "MyConnectome",
        "sub-01", "ses-105", "func",
        "sub-01_ses-105_task-rest_run_010_mask.nii.gz")
fileMNI2 <- file.path(fsl_std_dir(),
                     "MNI152_T1_2mm_brain_mask.nii.gz")
fr0 <- "sub-01_ses-105_task-rest_run-001_bold"
mfprest0 <- file.path(resDir, "MyConnectome", "sub-01",
              "ses-105", "func", paste0("meanpmc",fr0))
## writeNIfTI requires file names without extension
fr0 <- paste0(fr0,".nii.gz")
mnirestfile <- file.path(resDir, "MyConnectome",
     "sub-01", "ses-105", "func", paste0("nprmc", fr0))
mfprest <- paste0(mfprest0,".nii.gz")
fprest1 <- file.path(resDir, "MyConnectome", "sub-01",
                "ses-105", "func", paste0("prmc", fr0))
mnirestfile <- file.path(resDir, "MyConnectome",
                         "sub-01", "ses-105",
                         "func", paste0("nprmc", fr0))
fnTalairach <- file.path(fsl_atlas_dir(),
       "Talairach","Talairach-labels-2mm.nii.gz")


## ----label="slice timing version 1", echo=TRUE, eval=FALSE, warning=FALSE-----------------------------------
niirun1.o <- readNIfTI(f117s1r1, reorient = FALSE)
niirun1 <- oro2fmri(niirun1.o, setmask = FALSE)
sltimes <- c(seq(1, 33, 2), seq(2, 32, 2))
stcniirun1 <- slicetiming(niirun1, sltimes)
stcdata <- extractData(stcniirun1)
stcdata.nii <- as.nifti(stcdata, niirun1.o)
writeNIfTI(stcdata.nii, file = f117s1r1stc)


## ----label="Use ANTsR for motion correction1", echo=TRUE, eval=FALSE, results=FALSE-------------------------
fmrirun1 <- antsImageRead(f105s1r1)
fmrirun1MC <- antsrMotionCalculation(fmrirun1,
                    typeofTransform = "BOLDRigid")


## ----label="Use ANTsR for motion correction2", echo=TRUE, eval=FALSE, results=FALSE-------------------------
badtimes <- which(fmrirun1MC$fd$MeanDisplacement > .2 )
## motion corrected time series
antsImageWrite(fmrirun1MC$moco_img, f105s1r1mc)
## average corrected image
if(packageVersion("ANTsR") == '0.6.1'){
  avimg <- as.antsImage(fmrirun1MC$moco_avg_img,
      spacing = antsGetSpacing(fmrirun1)[1:3],
      origin = antsGetOrigin(fmrirun1)[1:3],
      direction = antsGetDirection(fmrirun1)[1:3, 1:3])
} else{
  avimg <- fmrirun1MC$moco_avg_img
}
antsImageWrite(avimg, f105s1r1amc)
## save estimated parameters
mocoparams <- fmrirun1MC$moco_params
save(badtimes, mocoparams, file = f105s1r1motion)


## ----label="filenames with ANTsR 1", echo=FALSE, eval=TRUE--------------------------------------------------
avimg <- antsImageRead(f105s1r1amc)


## ----label="Figure_4_1", echo=FALSE, eval=TRUE, fig.cap="ANTsR rigid motion correction parameters. The left column shows the three parameters related to rotation, the right column shows the time course of the translation parameters.",fig.width=7,fig.height=5----
load(f105s1r1motion)
par(mfcol=c(3,2),mar=c(3,3,3,1),mgp=c(2,1,0))
scans <- dim(mocoparams)[1]
xl <- "frame"
yl <- "Rotation"
plot( 1:scans, mocoparams[,1], type="l",
      main="MOCOparam1 A(1,2)", xlab = xl, ylab= yl)
plot( 1:scans, mocoparams[,2], type="l",
      main="MOCOparam2 A(1,3)", xlab = xl, ylab= yl)
plot( 1:scans, mocoparams[,3], type="l",
      main="MOCOparam3 A(2,3)", xlab = xl, ylab= yl)
yl <- "Translation"
plot( 1:scans, mocoparams[,4], type="l",
      main="MOCOparam4 b(1)", xlab = xl, ylab= yl)
plot( 1:scans, mocoparams[,5], type="l",
      main="MOCOparam5 b(2)", xlab = xl, ylab= yl)
plot( 1:scans, mocoparams[,6], type="l",
      main="MOCOparam6 b(3)", xlab = xl, ylab= yl)


## ----label="ANTsR Registration", echo=TRUE, eval=FALSE, results=FALSE, message=FALSE------------------------
T1 <- antsImageRead(f105s1T1)
avimg <- antsImageRead(f105s1r1amc)
amcf2T1 <- antsRegistration(fixed = T1,
                            moving = avimg,
                  typeofTransform = "SyNBoldAff")
antsImageWrite(amcf2T1$warpedmovout, f105s1bold2T1)


## ----label="Figure_4_2", echo=-(1:3), eval=TRUE, results=FALSE, fig.cap="ANTsR-registration: Mean fMRI image (top), T1 image (center) and Mean fMRI image registered to T1 image (bottom).",fig.width=7,fig.height=8----
avimg <- antsImageRead(f105s1r1amc)
T1 <- antsImageRead(f105s1T1)
bold2T1 <- antsImageRead(f105s1bold2T1)
par(mfrow = c(3, 1),
    mar = c(1, 2, 4, 1), mgp = c(3, 2, 0))
invisible(plot(avimg, slices = c(24, 29, 34, 39),
               nslices = 4))
invisible(plot(T1, slices = c(98, 118, 138, 158),
               nslices = 4))
invisible(plot(bold2T1, slices = c(98, 118, 138, 158),
               nslices = 4))


## ----label="ANTsR Normalization", echo=TRUE, eval=FALSE, results=FALSE--------------------------------------
avimg0 <- antsImageRead(f105s1r1mc)
# T1w volume has incorrect direction and origin
antsSetDirection(T1, antsGetDirection(avimg)[1:3, 1:3])
antsSetOrigin(T1, antsGetOrigin(avimg)[1:3])
mni <- antsImageRead(fileMNI)
T1toMNI <- antsRegistration(fixed = mni,moving = T1,
                            typeofTransform = "SyN")
mcfmri2T1 <- antsApplyTransforms(T1, avimg0,
                  amcf2T1$fwdtransforms, imagetype = 3)
mcfmri2MNI <- antsApplyTransforms(mni, mcfmri2T1,
                  T1toMNI$fwdtransforms, imagetype = 3)
antsImageWrite(mcfmri2MNI, f105s1r1nmc)


## ----label="Brain_mask using fslr",echo=TRUE, eval=FALSE, cache=FALSE, results="hide"-----------------------
fslbet(f105s1T1, outfile = f105s1brain, retimg = TRUE,
         reorient = FALSE, betcmd = "bet2", opt = "-m")


## ----label="Brain mask in subject space", echo=TRUE, eval=FALSE, cache=FALSE--------------------------------
mask <- antsImageRead(f105s1mask)
fmrimask <- antsApplyTransforms(avimg, mask,
                                amcf2T1$invtransforms)
antsImageWrite(fmrimask, f105s1boldmask)


## ----label="T1segmentWithFSL", echo= TRUE, eval=FALSE-------------------------------------------------------
fast(f105s1brain, outfile = f105s1brainFAST)


## ----label="ANTsRCore segmentation", echo=-(1:2), eval=FALSE, cache=FALSE-----------------------------------
T1 <- antsImageRead(f105s1T1)
mask <- antsImageRead(f105s1mask)
T1 <- n3BiasFieldCorrection(T1, 4)
simg <- kmeansSegmentation(T1, 3, mask)


## ----label="T1segmentWithMRITC", echo=-c(1:2), eval=FALSE---------------------------------------------------
T1 <- readNIfTI(f105s1T1)
mask <- readNIfTI(f105s1mask)
tc.icm <- mritc(T1, mask, method = "ICM")
gm <- wm <- csf <- array(0, dim = dim(mask))
csf[mask == 1] <- tc.icm$prob[, 1]
gm[mask == 1] <- tc.icm$prob[, 2]
wm[mask == 1] <- tc.icm$prob[, 3]


## ----label="Atlas mapping to subject space", echo=TRUE, eval=FALSE, results=FALSE, cache=FALSE--------------
T1 <- antsImageRead(f105s1T1)
mni <- antsImageRead(fileMNI)
T1toMNI <- antsRegistration(fixed = mni,
              moving = T1, typeofTransform = "SyN")
shen268 <- antsImageRead(fileshen268)
shen268mni <- antsRegistration(fixed = mni,
              moving = shen268,
              typeofTransform = "Affine")$warpedmovout
antsImageWrite(shen268mni, fileshen268mni)
# Shen atlas has a different orientation than MNI
shen268sub1T1 <- antsApplyTransforms(T1,
                     shen268mni, T1toMNI$invtransforms,
                     interpolator = "nearestNeighbor")
avimg <- antsImageRead(f105s1r1amc)
shen268sub1 <- antsApplyTransforms(avimg,
                  shen268sub1T1, amcf2T1$invtransforms,
                  interpolator = "nearestNeighbor")
antsImageWrite(shen268sub1, f105s1r1atlas268)


## ----label="Figure_4_3", echo=-5, eval=TRUE, fig.cap="ANTsR-registration: Shen268 atlas parcellation in MNI (upper row) and fMRI subject space (lower row).",fig.width=7,fig.height=3.6----
shen268sub1 <- antsImageRead(f105s1r1atlas268)
shen268mni <- antsImageRead(fileshen268mni)
par(mfrow=c(2,1),mar=c(1,2,4,1),mgp=c(3,2,0))
invisible( plot(shen268mni,
                slices=c(61,81,101,121), nslices=4))
invisible( plot(shen268sub1,
                slices=c(24,29,34,39), nslices=4))


## ----label="Gaussian smoothing", echo=TRUE, eval=TRUE, results=FALSE, cache=TRUE----------------------------
nii <- snii <- readNIfTI(f105s1r1mc, reorient = FALSE)
bw <- 8/nii@pixdim[2:4]
for (i in 1:dim(nii)[4])
   snii@.Data[ , , , i] <- kernsm(nii[, , , i],
                              h = bw, unit="FWHM")@yhat


## ----label="Figure_4_4", echo=FALSE, eval=TRUE, results=FALSE,out.width='66%',fig.cap="Left: original slice. Right: after smoothing with a Gaussian filter with FWHM bandwidth of 8mm.", fig.width=6, fig.height=5.2, fig.align='center'----
par(mfrow=c(1,2),mar=c(1,1,3,.1),mgp=c(2,1,0),cex=1)
image(nii@.Data[ , , 32, 100],col=grey(0:255/255),xaxt="n",yaxt="n",xlab="",ylab="",asp=1.6)
title("Slice 32 original")
image(snii@.Data[ , , 32, 100],col=grey(0:255/255),xaxt="n",yaxt="n",xlab="",ylab="",asp=1.6)
title("Slice 32 smoothed")


## ----label="Figure_4_5", echo=FALSE, eval=TRUE, results=FALSE,fig.cap="Schematic view of the characteristic hemodynamic response function (HRF) $h(t)$ following a short stimulus. This principle behavior of the BOLD-signal has been first observed in the occipital lobe of the brain for visual stimuli, see \\citet{Huettel2009}.", fig.width=10,fig.height=4.3----
par(mfrow=c(1,1),mar=c(3,3,3,.1),mgp=c(2,1,0))
canonicalHRF <- function(t, par = c(6, 12, 0.9, 0.9, 0.35)) {
        ttpr <- par[1] * par[3]
        ttpu <- par[2] * par[4]
        ((t/ttpr)^par[1] * exp(-(t - ttpr)/par[3]) - par[5] *
            (t/ttpu)^par[2] * exp(-(t - ttpu)/par[4]))/2.885802
    }
gammaHRF <- function(t, par = 4) {
        th <- 0.242 * par[1]
        1/(th * factorial(3)) * (t/th)^3 * exp(-t/th)
    }
tt <- seq(0,20,.1)
plot(tt,canonicalHRF(tt),type="l",ylab="HRF",xlab="time (s)",lwd=2)
lines(tt,gammaHRF(tt),col=2,lwd=2)
legend(10,.3,c("Canonical HRF","Gamma HRF"),col=1:2,lwd=c(2,2),lty=c(1,1))


## ----label="Create expected BOLD", echo=TRUE, eval=FALSE, results=FALSE-------------------------------------
## bold <- convolve(stimulus, rev(hrf), type="open")


## ----label="Create expected BOLD with fmri", echo=TRUE, eval=TRUE, results=FALSE----------------------------
hrfc <- fmri.stimulus(105, c(16, 46, 76), 15, 2,
                      type = "canonical")
hrfg <- fmri.stimulus(105, c(16, 46, 76), 15, 2,
                      type = "gamma")
hrfb <- fmri.stimulus(105, c(16, 46, 76), 15, 2,
                      type = "boxcar")


## ----label="Figure_4_6", echo=FALSE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Expected BOLD response for gamma HRF and canonical HRF. For the gamma HRF the effect of different choices for $\\tau_h$ are shown. The green line is the stimulus indicator function.",fig.width=6.5,fig.height=2.8----
par(mfrow=c(1,3),mar=c(3,3,3,1),mgp=c(2,1,0))
stim <- .25*c(rep(0,15),rep(1,15),rep(0,15),rep(1,15),rep(0,15),rep(1,15),rep(0,15))
scans <- 1:105
plot(scans,hrfc,type="l",main="canonical HRF",ylab="stimulus/expected BOLD response")
lines(scans,stim,col="green")
plot(scans,hrfg,type="l",main="gamma HRF",ylab="stimulus/expected BOLD response")
lines(scans,stim,col="green")
plot(scans,hrfb,type="l",main="boxcar HRF",ylab="stimulus/expected BOLD response")
lines(scans,stim,col="green")


## ----label="create linear model1", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE-------------------------
nii <- readNIfTI(f105s1r1mc, reorient = FALSE)
TR <- 2.5
scans <- dim(nii)[4]


## ----label="create linear model2", echo=TRUE, eval=TRUE, results=FALSE,cache=FALSE--------------------------
ttt <- read.table(f105s1r1events, header = TRUE)
ntrials <- dim(ttt)[1]
indScissors <- (1:ntrials)[ttt$trial_type == "scissors"]
indFace <- (1:ntrials)[ttt$trial_type == "face"]
indCat <- (1:ntrials)[ttt$trial_type == "cat"]
indShoe <- (1:ntrials)[ttt$trial_type == "shoe"]
indHouse <- (1:ntrials)[ttt$trial_type == "house"]
indScrb <- (1:ntrials)[ttt$trial_type == "scrambledpix"]
indBottle <- (1:ntrials)[ttt$trial_type == "bottle"]
indChair <- (1:ntrials)[ttt$trial_type == "chair"]
onsets <- ttt$onset
duration <- ttt$duration


## ----label="create linear model3", echo=TRUE, eval=TRUE, results=FALSE--------------------------------------
hrfScissors <- fmri.stimulus(scans, onsets[indScissors],
                             duration[indScissors],
                             TR = TR, times = TRUE)
hrfFace <- fmri.stimulus(scans, onsets[indFace],
                         duration[indFace],
                         TR = TR, times = TRUE)
hrfCat <- fmri.stimulus(scans, onsets[indCat],
                        duration[indCat],
                        TR = TR, times=TRUE)
hrfShoe <- fmri.stimulus(scans, onsets[indShoe],
                         duration[indShoe],
                         TR = TR, times = TRUE)
hrfHouse <- fmri.stimulus(scans, onsets[indHouse],
                          duration[indHouse],
                          TR = TR, times = TRUE)
hrfScrambled <- fmri.stimulus(scans, onsets[indScrb],
                              duration[indScrb],
                              TR = TR, times = TRUE)
hrfBottle <- fmri.stimulus(scans, onsets[indBottle],
                           duration[indBottle],
                           TR = TR, times = TRUE)
hrfChair <- fmri.stimulus(scans, onsets[indChair],
                          duration[indChair],
                          TR = TR, times = TRUE)
hrf <- cbind(hrfScissors,hrfFace,hrfCat,hrfShoe,
             hrfHouse,hrfScrambled,hrfBottle,hrfChair)


## ----label="create linear model4", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE-------------------------
load(f105s1r1motion)
xdesign <- fmri.design(hrf, order = 2, cef = mocoparams)


## ----label="Figure_4_7", echo=FALSE, eval=TRUE, results=FALSE, fig.cap="Design matrix: HRF, nuisance parameters, trends.", fig.width=7.5,fig.height=3----
par(mfrow=c(1,3),mar=c(3,2,3,1),mgp=c(2,1,0))
plot(c(1, scans), range(xdesign[,1:8]), type = "n", xlab="Frame", ylab="",main="Expected BOLD response")
for (i in 1:8) lines(1:scans, xdesign[, i], col = i, lwd = 2)
plot(c(1, scans), range(xdesign[,9:14]), type = "n", xlab="Frame", ylab="",main="Nuisance")
for (i in 1:6) lines(1:scans, xdesign[, i+8], col = i, lwd = 2)
plot(c(1, scans), range(xdesign[,15:17]), type = "n", xlab="Frame", ylab="",main="Trend")
for (i in 1:3) lines(1:scans, xdesign[, i+14], col = i, lwd = 2)


## ----label="evaluate GLM", echo=TRUE, eval=TRUE, results=FALSE----------------------------------------------
nii <- readNIfTI(f105s1r1mc, reorient = FALSE)
mask <- readNIfTI(f105s1boldmask) > 0
ds <- oro2fmri(nii)
ds$mask <- mask
contrast <- rep(1, 8) # Image versus rest
spm <- fmri.lm(ds, xdesign, contrast = contrast, 
               mask=mask)


## ----label="design for simulated data", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE--------------------
nscan <- 105
TR <- 2
on1 <- on2 <- c(16, 46, 76) * TR
onsets <- list(list(on1), list(on2))
dur <- list(list(15*TR), list(15*TR))


## ----label="temporal noise", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE-------------------------------
design <- simprepTemporal(totaltime = nscan*TR, 
                          regions = 2, onsets = onsets,
                          durations = dur, TR = TR,
                          effectsize = list(40, 40),
                          hrf = "double-gamma")


## ----label="spatial noise", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE--------------------------------
rc1 <- c(5, 5, 5)
rc2 <- c(12, 12, 12)
regions <- simprepSpatial(regions = 2,
                    coord = list(rc1, rc2),
                    radius = c(2, 4), form = "sphere")


## ----label="create simulated data", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE------------------------
ds.sim <- simVOLfmri(design = design, image = regions,
                     base = 100, dim = c(20, 20, 20),
                     noise = "mixture", type = "rician",
                     spat = "gaussRF",
                     weights = c(0.1, 0.3, 0.2,
                                 0.1, 0.1, 0.2),
                     rho.temp = 0.3, SNR = 7)


## ----label="Figure_4_8", echo=FALSE, eval=TRUE, results=FALSE, fig.cap="Artificial data created by \\pkg{neuRosim} in a voxel outside the activated regions and in a voxel within.",fig.width=10,fig.height=4----
par(mfrow = c(1, 2), mar = c(3, 3, 3, 0.1), mgp=c(2,1,0))
plot(ds.sim[1, 1, 1,], type = "l",xlab="Frame",ylab="BOLD signal",main="Background")
plot(ds.sim[12, 12, 12,], type = "l",xlab="Frame",ylab="BOLD signal",main="Active region")


## ----label="voxelwise inference", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE--------------------------
alpha <- 0.05
thr <- qt(1-alpha, spm$df)
theta <- spm$cbeta/sqrt(spm$var)
pvalue <- pt(-theta, spm$df)
pvalue[pvalue >= alpha] <- NA
dim(pvalue) <- spm$dim[1:3]


## ----label="Figure_4_9", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Voxelwise signal detection.", fig.width=10,fig.height=3.3----
bimg <- readNIfTI(f105s1T1)
pvalue <- fmri.pvalue(spm, "voxelwise")
plot(pvalue, bimg,
     view = "axial", ncol = 4, slices = 26:29)


## ----label="bonferroni", echo=TRUE, eval=TRUE, results=FALSE------------------------------------------------
alphaC <- alpha/sum(spm$mask)
thrC <- qt(1-alphaC, spm$df)


## ----label="Figure_4_10", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Voxelwise signal detection with Bonferrroni correction.", fig.width=10,fig.height=3.3----
pvalue <- fmri.pvalue(spm, "Bonferroni")
plot(pvalue, bimg,
     view = "axial", ncol = 4, slices = 26:29)


## ----label="smooth and reevaluate GLM", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, warning=FALSE-----
nii <- readNIfTI(f105s1r1mc, reorient = FALSE)
nii.s <- nii
for (i in 1:dim(nii.s)[4])
  nii.s@.Data[ , , , i] <- kernsm(nii[ , , , i],
                              h = 8/nii@pixdim[2:4],
                              unit="FWHM")@yhat
# thats 8mm spatial FWHM bandwidth
ds.s <- oro2fmri(nii.s)
ds.s$mask <- mask
spm.s <- fmri.lm(ds.s, xdesign, contrast = contrast, mask=mask)


## ----label="Figure_4_11", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Signal detection by RFT.", fig.width=10,fig.height=3.3----
pvaluesm <- fmri.pvalue(spm.s)
plot(pvaluesm, bimg,
     view = "axial", ncol = 4, slices = 26:29)


## ----label="Figure_4_12", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Signal detection by RFT (orthographic view).", fig.width=8,fig.height=2.9----
plot(pvaluesm, bimg, mask, view = "orthographic")


## ----label="FDR", echo=TRUE, eval=TRUE, results=FALSE-------------------------------------------------------
 fdr <- 0.05
 ind <- spm$cbeta > 0 & spm$mask
 theta <- spm$cbeta[ind]/sqrt(spm$var[ind])
 pvalue <- pt(-theta, spm$df)
 n <- length(pvalue)
 ind <- 1:n
 oind <- order(pvalue)
 nind <- length(ind[pvalue[oind] <= ind/n*fdr])
 oind <- oind[1:nind]
 thresh <- min(theta[oind])


## ----label="Figure_4_13", echo=TRUE, eval=TRUE, results=FALSE, fig.cap="Voxelwise signal detection by FDR.", fig.width=10,fig.height=3.3----
pvaluefdr <- fmri.pvalue(spm, mode = "FDR", alpha = .01)
plot(pvaluefdr, bimg,
     view = "axial", ncol = 4, slices = 26:29)


## ----label="p-value with cluster thresholding", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE------------
pvaluecl <-fmri.cluster(spm, alpha =.01, 
                        ncmin =8, ncmax=15)


## ----label="Figure_4_14", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Cluster based signal detection.", fig.width=10,fig.height=3.3----
plot(pvaluecl, bimg,
     view = "axial", ncol = 4, slices = 26:29)


## ----label="adaptive smoothing with fmri", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE-----------------
spm.smooth <- fmri.smooth(spm, hmax = 3,
                          adaptation = "awsfull")


## ----label="p-values with aws fmri", echo=TRUE, eval=TRUE, results=FALSE------------------------------------
pvalueaws <- fmri.pvalue(spm.smooth, alpha=.01)


## ----label="Figure_4_15", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Signal detection using structural adaptive smoothing \\citep{fmri_tabelow06}.", fig.width=10,fig.height=3.3----
plot(pvalueaws, bimg, view="axial", 
     ncol=4, slices=26:29)


## ----label="segment with fmri", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE----------------------------
spm.segment <- fmri.smooth(spm, hmax = 3,
                  adaptation = "segment", alpha = 0.01)


## ----label="segment results", echo=TRUE, eval=TRUE, results=FALSE-------------------------------------------
 signal <- spm.segment$cbeta
 signal <- signal * spm.segment$segm
 signal[signal <= 0] <- NA
  ## image with a suitable colorscale


## ----label="Figure_4_16", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Signal detection using structural adaptive segmentation \\citep{polzehletal10a}.", fig.width=8,fig.height=3.6----
par(mfrow = c(1, 4),
    mar = c(1, 1, 3, .1), mgp = c(1.5, .75, 0))
x <- 1:40; y <- 1:64
for (i in 26:29) {
  rimage(x, y, nii[, , i, 1], main = paste("Slice", i))
  rimage(x, y, signal[, , i], add = TRUE,
         col = heat.colors(256))
}


## ----"Figure_4_17", eval=TRUE, echo=FALSE, rgl = TRUE, cache=TRUE, fig.cap = "3D visualization of searchlights with radii 1.5, 2, 2.5, 3 and 3.5.", fig.width = 6, fig.height = 1.2, fig.align='center',results='hide',messages=FALSE----
open3d()
par3d(zoom=.3)
obj <- cube3d()
coords <- 2*t(fmri:::searchlight(1.5))
for(i in 1:dim(coords)[1]) shade3d(translate3d(obj, coords[i,1]-20, coords[i,2], coords[i,3]), color="gold")

coords <- 2*t(fmri:::searchlight(2))
for(i in 1:dim(coords)[1]) shade3d(translate3d(obj, coords[i,1], coords[i,2], coords[i,3]), color="red")

coords <- 2*t(fmri:::searchlight(2.5))
for(i in 1:dim(coords)[1]) shade3d(translate3d(obj, coords[i,1]+20, coords[i,2], coords[i,3]), color="cyan")

coords <- 2*t(fmri:::searchlight(3))
for(i in 1:dim(coords)[1]) shade3d(translate3d(obj, coords[i,1]+40, coords[i,2], coords[i,3]), color="blue")

coords <- 2*t(fmri:::searchlight(3.5))
for(i in 1:dim(coords)[1]) shade3d(translate3d(obj, coords[i,1]+60, coords[i,2], coords[i,3]), color="green")


## ----label="MVPA with fmri", echo=TRUE, eval=TRUE, results=FALSE--------------------------------------------
nii <- readNIfTI(f105s1r1mc, reorient = FALSE)
ds <- oro2fmri(nii)
ds$mask <- mask
cfaces <- c(-1/7, 1, rep(-1/7,6))
chouses <- c(rep(-1/7, 4), 1, rep(-1/7, 3))
spm.faces <- fmri.lm(ds, xdesign, contrast = cfaces, mask=mask)
spm.houses <- fmri.lm(ds, xdesign, contrast = chouses, mask=mask)
pvaluesl.faces <- fmri.searchlight(spm.faces,
                                   alpha = .01,
                                   radius = 2,
                                   kind = "abs")
pvaluesl.houses <- fmri.searchlight(spm.houses,
                                    alpha = .01,
                                    radius = 2,
                                    kind = "abs")


## ----label="Figure_4_18", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Results of the information-based activation detection for the faces contrast.", fig.width=10,fig.height=3.3----
plot(pvaluesl.faces, bimg,
     view = "axial", ncol = 4, slices = 26:29)


## ----label="Figure_4_19", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Results of the information-based activation detection for the houses contrast.", fig.width=10,fig.height=3.3----
plot(pvaluesl.houses, bimg,
     view = "axial", ncol = 4, slices = 26:29)


## ----"get pattern in voxel v", echo=TRUE, eval=TRUE---------------------------------------------------------
pattern.faces <- getSearchlightPattern(spm.faces,
                      pvaluesl.faces$pvalue < 0.05,
                      radius = 2)
dim(pattern.faces)


## ----label="ANTsR pre-processing for group ICA1", echo=TRUE, eval=FALSE, results=FALSE----------------------
fmrirest <- antsImageRead(frest)
pfmrirest <- preprocessfMRI(fmrirest)
antsImageWrite(pfmrirest$cleanBoldImage,fprest)
antsImageWrite(pfmrirest$maskImage, fprestmask)
avimg <-getAverageOfTimeSeries(pfmrirest$cleanBoldImage)


## ----label="ANTsR pre-processing for group ICA2", echo=TRUE, eval=FALSE, results=FALSE----------------------
for( i in 1:9){
  fr <- paste0("sub-01_ses-105_task-rest_run-00",
               i, "_bold.nii.gz")
  file <- file.path(dataDir, "MyConnectome",
              "sub-01", "ses-105", "func", fr)
  fmrirest <- antsImageRead(file)
  pfmrirest <- preprocessfMRI(fmrirest,
            meanBoldFixedImageForMotionCorrection=avimg)
  pfile <- file.path(resDir, "MyConnectome",
        "sub-01", "ses-105", "func", paste0("prmc", fr))
  antsImageWrite(pfmrirest$cleanBoldImage, pfile)
  fileprestmask <- file.path(resDir, "MyConnectome",
                    "sub-01", "ses-105", "func",
              paste0("sub-01_ses-105_task-rest_run_00",
                     i, "_mask.nii.gz"))
  antsImageWrite(pfmrirest$maskImage, fileprestmask)
}


## ----label="sICA", echo=TRUE, eval=TRUE, results=FALSE, cache=TRUE------------------------------------------
mask <- as.array(readNIfTI(fprestmask, reorient=FALSE))
nii <- readNIfTI(frest, reorient = FALSE)
ds <- oro2fmri(nii)
ICAresult <- fmri.sICA(ds, mask, ncomp = 40, 
          degree = 3, bws = 8, bwt = 4, unit = "FWHM")


## ----label="ICAfingerprint", echo=TRUE, eval=TRUE, results=FALSE, cache=TRUE--------------------------------
ICAresult <- ICAfingerprint(ICAresult)


## ----label="Figure_4_20", echo=TRUE, eval=TRUE, results=FALSE, fig.cap="Diagnostic plot for IC evaluation.", fig.width=10,fig.height=5.,warning=FALSE----
plot(ICAresult, 21)


## ----label="sICA group study", echo=TRUE, eval=TRUE, results=FALSE, cache=TRUE------------------------------
ICAresults <- list(NULL)
for(i in 1:9){
  fr <- paste0("prmcsub-01_ses-105_task-rest_run-00",
               i, "_bold.nii.gz")
  file <- file.path(resDir, "MyConnectome",
                    "sub-01", "ses-105",
                    "func", fr)
  nii <- readNIfTI(file, reorient = FALSE)
  ds <- oro2fmri(nii)
  ICAresults[[i]] <- fmri.sICA(ds, mask,
                               ncomp = 40, degree = 3,
                               bws = 8, bwt = 4,
                               unit = "FWHM",
                               alg.typ = "parallel",
                               alpha = 1.)
}
ICAresults[[10]] <- ICAresult
ICallruns <- fmri.sgroupICA(ICAresults, thresh = .9,
                           minsize = 6)


## ----label="Figure_4_21", echo=TRUE, eval=TRUE, results=FALSE, fig.cap="First independent component from combined IC analysis of 10 runs.", fig.width=9,fig.height=3.8,cache=FALSE,warning=FALSE----
plot(ICallruns, 1)


## ----"get mean BOLD image for registration",eval=TRUE,echo=TRUE,results=FALSE-------------------------------
pmcrest <- readNIfTI(fprest1, reorient = FALSE)
ameanpmcrest <- as.nifti(apply(abs(pmcrest), 1:3, mean),
                        pmcrest)
writeNIfTI(ameanpmcrest, mfprest0)


## ----"Register preprocessed resting data to MNI",eval=FALSE, echo=TRUE--------------------------------------
fmrirest <- antsImageRead(fprest1)
mni <- antsImageRead(fileMNI2)
prest2MNI <- antsRegistration(mni, mfprest, "Affine")
mnirest <- antsApplyTransforms(mni, fmrirest,
                prest2MNI$fwdtransforms, imagetype = 3)
antsImageWrite(mnirest, mnirestfile)


## ----"Access the Talairach atlas",eval=TRUE,echo=TRUE-------------------------------------------------------
atlas <- readNIfTI(fnTalairach)
rsdata <- readNIfTI(mnirestfile)
mask <- apply(abs(rsdata),1:3,max) > 0


## ----"read atlas metadata",echo=TRUE,eval=TRUE,warnings=FALSE-----------------------------------------------
Tmeta <- xmlTreeParse(file.path(fsl_atlas_dir(),
                                "Talairach.xml"))
Tmeta <- xmlToList(Tmeta)
regionNames <- as.character(Tmeta$data[1, ])
ind <- rep(0, length(regionNames))
getregion <- function(regionNames){
  ind <- rep(0, length(regionNames))
  for (i in 1:length(regionNames)){
    z <- strsplit(regionNames[i], "Brodmann area")
    if(length(z[[1]]) == 2){
      ind[i] <- as.numeric(z[[1]][2])
      if(length(grep("Left", z[[1]][1])) > 0)
        ind[i] <- ind[i] + 50
    }
  }
  ind
}
ind <- getregion(regionNames)
regions <- sort(unique(ind[ind > 0]))


## ----"aggregate atlas",eval=TRUE,echo=TRUE------------------------------------------------------------------
batlas <- atlas@.Data
batlas[, , ] <- 0
for(i in 1:97){
  jind <- (1:length(ind))[ind == i]
  batlas[atlas %in% jind] <- i
}
batlas[!mask] <- 0


## ----"create region names", eval=TRUE, echo=TRUE------------------------------------------------------------
ni <- numeric(length(regions))
for(i in 1:length(regions))
  ni[i] <- sum(batlas == regions[i])
regions <- regions[ni > 1]
nregions <- length(regions)
nL <- sum(regions > 50)
regionNames <- paste0(regions %% 50,
                      c(rep("R", nregions-nL),
                        rep("L", nL)))


## ----"create mean time series",eval=TRUE, echo=TRUE---------------------------------------------------------
mts <- matrix(0, 240, nregions)
dim(rsdata) <- c(prod(rsdata@dim_[2:4]),
                 rsdata@dim_[5])
for(i in 1:nregions)
   mts[, i] <- apply(rsdata[batlas == regions[i], ],
                     2, mean)
colnames(mts) <- regionNames


## ----"get correlation matrix", eval=TRUE, echo=TRUE---------------------------------------------------------
rcor <- cor(mts)


## ----"Figure_4_22",echo=TRUE,eval=TRUE,out.width='100%', fig.align="center",fig.cap="Functional  connectivity  between Brodmann areas: Heatmap of correlation matrix from resting state fMRI.",fig.width=10,fig.height=10----
heatmap(rcor, symm = TRUE)


## ----"get sparse partial correlation matrix",eval=TRUE,echo=TRUE--------------------------------------------
smat <- cov(mts)
pmat <- glasso(smat, .002)$wi
dsmath <- diag(1/sqrt(diag(pmat)))
pcor <- -dsmath %*% pmat %*% dsmath
diag(pcor) <- 1
dimnames(pcor) <- dimnames(rcor)


## ----"Construction of graphs",eval=TRUE,echo=TRUE-----------------------------------------------------------
diag(rcor) <- diag(pcor) <- 0
gcor <- graph_from_adjacency_matrix(rcor,
                                    "undirected",
                                    weighted = TRUE)
gpc <- graph_from_adjacency_matrix(pcor,
                                   "undirected",
                                   weighted = TRUE)


## ----"Figure_4_23",eval=TRUE,echo=TRUE,out.width='100%', fig.align="center",fig.cap="Functional connectivity between Brodmann areas: Connectivity graphs obtained from resting state fMRI using correlation (left) and partial correlation matrices (right).",fig.width=14,fig.height=8----
par(mfrow = c(1, 2),
    mar = c(3, 1, 3, .5), mgp = c(2, 1, 0))
colind <- min(255*(rcor+1)/2):max(255*(rcor+1)/2)
plot(permute(gcor,
             canonical_permutation(gcor)$labeling),
     layout = layout_in_circle, vertex.size = 12,
     edge.color = cm.colors(256)[colind])
title("functional connectivity (correlation)")
colind <- min(255*(pcor+1)/2):max(255*(pcor+1)/2)
plot(permute(gpc, canonical_permutation(gpc)$labeling),
     layout = layout_in_circle, vertex.size = 12,
     edge.color = cm.colors(256)[colind])
title("functional connectivity (partial correlation)")

