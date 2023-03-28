if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_fmri","chapter_fMRI_init.R"))

if(!haveFSLr||!haveANTsR) stop("need packages fslr and ANTsR installed \n precomputed objects  available in results")
## ----label="filenames with ANTsR 1", echo=FALSE, eval=TRUE--------------------------------------------------
avimg <- antsImageRead(f105s1r1amc)

## ----label="ANTsR Registration", echo=TRUE, eval=FALSE, results=FALSE, message=FALSE------------------------
T1 <- antsImageRead(f105s1T1)
avimg <- antsImageRead(f105s1r1amc)
amcf2T1 <- antsRegistration(fixed = T1,
                            moving = avimg,
                  typeofTransform = "SyNBold")
antsImageWrite(amcf2T1$warpedmovout, f105s1bold2T1)

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


