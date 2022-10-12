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

