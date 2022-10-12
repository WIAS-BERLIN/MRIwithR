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
