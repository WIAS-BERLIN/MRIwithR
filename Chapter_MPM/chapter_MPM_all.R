## ----label="initQMRI", echo=FALSE, eval=TRUE----
library(knitr)
knitr::opts_chunk$set(warning = FALSE,
                      message = FALSE)
options(width=50, digits=3)

## ----label="packages and filenames", eval=TRUE, echo=-1----
if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_MPM","chapter_MPM_init.R"))

## ----label="initMPM", echo=FALSE, eval=TRUE, message=FALSE, results=FALSE----
rimage.options(zquantiles=c(0.001,0.999), 
               xlab="x", ylab="z", bty="n")

## ----label="MPMdata", echo=TRUE, eval=TRUE------
t1Dir <- "t1w_mfc_3dflash_v1i_R4_0015"
pdDir <- "pdw_mfc_3dflash_v1i_R4_0009"
mtDir <- "mtw_mfc_3dflash_v1i_R4_0012"
t1Prefix <- "anon_s2018-02-28_18-26-190921-00001-"
pdPrefix <- "anon_s2018-02-28_18-26-185345-00001-"
mtPrefix <- "anon_s2018-02-28_18-26-190132-00001-"

## ----label="MPMdataT1", echo= TRUE, eval=TRUE----
t1ID <- c("00224-1.nii", "00448-2.nii", "00672-3.nii",
          "00896-4.nii", "01120-5.nii", "01344-6.nii",
          "01568-7.nii", "01792-8.nii")
t1Names <- paste0(t1Prefix, t1ID)
t1Files <- file.path(MPMdataDir, t1Dir, t1Names)

## ----label="MPMdataPD", echo= TRUE, eval=TRUE----
pdID <- c("00224-1.nii", "00448-2.nii", "00672-3.nii",
          "00896-4.nii", "01120-5.nii", "01344-6.nii",
          "01568-7.nii", "01792-8.nii")
pdNames <- paste0(pdPrefix, pdID)
pdFiles <- file.path(MPMdataDir, pdDir, pdNames)

## ----label="MPMdataMT", echo= TRUE, eval=TRUE----
mtID <- c("00224-1.nii", "00448-2.nii", "00672-3.nii",
          "00896-4.nii", "01120-5.nii", "01344-6.nii")
mtNames <- paste0(mtPrefix, mtID)
mtFiles <- file.path(MPMdataDir, mtDir, mtNames)

## ----label="defineBrainMaskForMPM", echo= TRUE, eval=TRUE, results=FALSE----
maskName <- paste0(t1Prefix, "00224-1_brain")
maskFile <- file.path(MPMresDir, t1Dir, maskName)

## ----"Precomputed", echo=FALSE, eval=TRUE-------
havemask <- file.exists(paste0(maskFile,".nii.gz"))
if(!havemask&&!haveFSLr) stop("No fslr installed, need precomputed mask")
registered <- dir.exists(file.path(MPMresDir,"pdw_mfc_3dflash_v1i_R4_0009"))
if(!registered&&!haveANTsR) stop("No ANTsR installed, need precomputed registered images")
registeredB1 <- file.exists(file.path(MPMresDir,"pdw_mfc_3dflash_v1i_R4_0009",
         "Results","Supplementary","anon_s2018-02-28_18-26-184837-00001-00001-1_B1mapReg.nii"))
haveMPMdata <- file.exists(file.path(MPMresDir,"MPMdata.rsc"))
havemodelMPM <- file.exists(file.path(MPMresDir,"modelMPM.rsc"))
havesmoothedModelMPM <- file.exists(file.path(MPMresDir,"smoothModelMPM.rsc"))
haveMPMmaps <- file.exists(file.path(MPMresDir,"qMaps.rsc"))

## ----label="ComputeBrainMaskForMPM", echo= TRUE, eval=!havemask, results=FALSE----
if(!havemask){
fslbet(t1Files[1], outfile = maskFile, retimg = TRUE,
       reorient = FALSE, opts = "-f 0.3 -m",
       betcmd = "bet2", verbose = FALSE)
}

## ----label="readMPMdata", echo= TRUE, eval=TRUE, cache=TRUE, message=FALSE----
mpm <- readMPMData(t1Files, pdFiles, mtFiles,
                   maskFile, verbose = FALSE)

## ----label="getMPMDataParameters", echo=TRUE, eval=TRUE----
mpm$TE
mpm$FA
mpm$TR
mpm$sdim

## ----label="Figure_6_1", echo= 2:6, eval=TRUE, fig.cap="Unregistered first T1, MT and PD images." , fig.width=7, fig.height=2.4----
oldpar <- par(mfrow = c(1, 3), mar = c(3, 3, 3, .1), mgp = c(2, 1, 0))
ddata <- extract(mpm,"ddata")
iy <- 160
rimage(ddata[1, , iy, ], main="first T1")
rimage(ddata[9, , iy, ], main="first MT")
rimage(ddata[15, , iy, ], main="first PD")
par(oldpar)

## ----label="MPMdataPDreg", echo= TRUE, eval=TRUE----
rpdNames <- paste0(pdPrefix,"r", pdID)
rpdFiles <- file.path(MPMresDir, pdDir, rpdNames)

## ----label="MPMdataMTreg", echo= TRUE, eval=TRUE----
rmtNames <- paste0(mtPrefix,"r",  mtID)
rmtFiles <- file.path(MPMresDir, mtDir, rmtNames)

## ----label="RegisterMPMdata1", echo=TRUE, eval=!registered----
if(!registered){
 T1 <- antsImageRead(t1Files[1])
 MT <- antsImageRead(mtFiles[1])
 PD <- antsImageRead(pdFiles[1])
 mask <- antsImageRead(paste0(maskFile, ".nii.gz"))
 aregMT2T1 <- antsRegistration(fixed = mask,
         moving = MT, typeofTransform = "DenseRigid",
         verbose = TRUE)

## ----label="RegisterMPMdata2", echo=TRUE, eval=!registered----
 for (id in 1:6){
   img <- antsImageRead(mtFiles[id])
   warpedMT <- antsApplyTransforms(fixed = T1,
                moving = img,
                transformlist = aregMT2T1$fwdtransforms)
   antsCopyImageInfo(img, warpedMT)
   antsImageWrite(warpedMT, rmtFiles[id])
 }

## ----label="RegisterMPMdata3", echo=TRUE, eval=!registered----
 aregPD2T1 <- antsRegistration(fixed = mask,
            moving = PD, typeofTransform = "DenseRigid",
            verbose = TRUE)
 for (id in 1:8){
   img <- antsImageRead(pdFiles[id])
   warpedPD <- antsApplyTransforms(fixed = T1,
                moving = img,
                transformlist = aregPD2T1$fwdtransforms)
   antsCopyImageInfo(img, warpedPD)
   antsImageWrite(warpedPD, rpdFiles[id])
 }
} else {
  
## ----"loadMPMdata", eval=haveMPMdata, echo=FALSE----
load(file.path(MPMresDir,"MPMdata.rsc"))
}

## ----label="readMPMRegdata", echo= TRUE, eval=!haveMPMdata, message=FALSE----
if(!haveMPMdata){
 mpm <- readMPMData(t1Files, rpdFiles, rmtFiles,
                    maskFile, TR = mpm$TR, TE = mpm$TE,
                    FA = mpm$FA, verbose = FALSE)
 save(mpm, file=file.path(MPMresDir, "MPMdata.rsc"))
}

## ----label="Figure_6_2", echo= 2:5, eval=TRUE, fig.cap="Registered first $T_1$, MT and PD images, cf. Fig. \\ref{fig:Figure_6_1}.", fig.width=7, fig.height=2.4----
oldpar <- par(mfrow = c(1, 3), mar = c(3, 3, 3, .1), mgp = c(2, 1, 0))
ddata <- extract(mpm,"ddata")
rimage(ddata[1, , iy, ], main = "first T1")
rimage(ddata[9, , iy, ], main = "first MT")
rimage(ddata[15, , iy, ], main = "first PD")
par(oldpar)

## ----label="Figure_6_3", echo=-1, eval=TRUE, fig.cap="Mean signal (over all voxel within a brain mask) for the T1w, MTw, and PDw multi-echo sequence.", fig.width=7, fig.height=3----
oldpar <- par(mfrow = c(1, 3), mar = c(3,3,3,1), mgp=c(2,1,0))
dim(ddata) <- c(length(mpm$TE), prod(mpm$sdim))
mask <- as.vector(extract(mpm,"mask"))
meanData <- apply(ddata[, mask], 1, mean)
plot(mpm$TE[1:8], meanData[1:8],
     xlab="TE in ms", ylab="Average signal", main="T1w")
plot(mpm$TE[9:14], meanData[9:14],
     xlab="TE in ms", ylab="Average signal", main="MTw")
plot(mpm$TE[15:22], meanData[15:22],
     xlab="TE in ms", ylab="Average signal", main="PDw")

## ----label="Reset",echo=FALSE,eval=TRUE---------
par(oldpar)

## ----label="loadModelMPM", echo=FALSE, eval=havemodelMPM----
if(havemodelMPM){
 load(file.path(MPMresDir,"modelMPM.rsc"))
} else {

## ----label="estimateMPM", echo=TRUE, eval=!havemodelMPM----
 modelMPM <- estimateESTATICS(mpm, verbose = TRUE)
}

## ----label="inspectMPMobject1", echo=TRUE, eval=TRUE----
dim(extract(modelMPM,"modelCoeff"))
dim(extract(modelMPM,"invCov"))

## ----label="inspectMPMobject", echo=TRUE, eval=TRUE----
modelMPM$TEScale
modelMPM$dataScale

## ----label="loadMaps", echo=FALSE, eval=haveMPMmaps----
if(haveMPMmaps){
 load(file.path(MPMresDir,"qMaps.rsc"))
} else {

## ----label="calculateMaps", echo=TRUE, eval=!haveMPMmaps, results=FALSE----
 qMRImaps <- calculateQI(modelMPM, TR2 = 3.4)
}

## ----label="inspectQMRIobject", echo=TRUE, eval=TRUE----
names(qMRImaps)
dim(extract(qMRImaps,"R1"))

## ----label="Figure_6_4", echo=2:10, eval=TRUE, fig.cap="qMRI maps generated from the ESTATICS model.", fig.width=7, fig.height=6.5----
oldpar <- par(mfrow=c(2, 2), mar=c(3,3,3,1), mgp=c(2,1,0))
ix <- 28:245; iy <- 160; iz <- 23:197
R1 <- extract(qMRImaps, "R1")[ix, iy, iz]
rimage(ix, iz, R1, zlim = c(0, 1.75), main = "R1")
R2star <- extract(qMRImaps, "R2star")[ix, iy, iz]
rimage(ix, iz, R2star, zlim = c(0, 50), main = "R2*")
PD <- extract(qMRImaps, "PD")[ix, iy, iz]
rimage(ix, iz, PD, zlim = c(0, 15000), main = "A")
MT <- extract(qMRImaps, "MT")[ix, iy, iz]
rimage(ix, iz, MT, zlim = c(0, 2), main = "MT", 
       col = colMT)
par(oldpar)

## ----label="registerBiasFieldMaps1", echo=TRUE, eval=!registeredB1----
if(!registeredB1){
 PD <- antsImageRead(pdFiles[1])

## ----label="registerBiasFieldMaps2", echo=TRUE, eval=!registeredB1----
 b1prefix <- "anon_s2018-02-28_18-26-184837-00001-"
 B1refName <- paste0(b1prefix, "00001-1_B1ref.nii")
 B1refFile <- file.path(MPMresDir, pdDir,
                        "Results", "Supplementary",
                        B1refName)
 B1ref <- antsImageRead(B1refFile)

## ----label="registerBiasFieldMaps3", echo=TRUE, eval=!registeredB1----
 B1mapName <- paste0(b1prefix, "00001-1_B1map.nii")
 B1mapFile <- file.path(MPMresDir, pdDir,
                        "Results", "Supplementary",
                        B1mapName)
 B1map <- antsImageRead(B1mapFile)

## ----label="registerBiasFieldMaps4", echo=TRUE, eval=!registeredB1----
 B1refReg <- antsRegistration(fixed = PD,
                              moving = B1ref,
                              typeofTransform = "SyN")
 B1refRegName <- paste0(b1prefix,
                       "00001-1_B1refReg.nii")
 B1refRegFile <- file.path(MPMresDir, pdDir,
                           "Results", "Supplementary",
                           B1refRegName)
 antsImageWrite(B1refReg$warpedmovout, B1refRegFile)

## ----label="registerBiasFieldMaps5", echo=TRUE, eval=!registeredB1----
 B1mapReg <- antsApplyTransforms(PD,
                                 B1map,
                                 B1refReg$fwdtransforms,
                                 interpolator = "linear")
 B1mapRegName <- paste0(b1prefix, "00001-1_B1mapReg.nii")
 B1mapRegFile <- file.path(MPMresDir, pdDir,
                           "Results", "Supplementary",
                           B1mapRegName)
 antsImageWrite(B1mapReg, B1mapRegFile)
}

## ----label="calculateMapsFTcorr", echo=TRUE, eval=!haveMPMmaps, results=FALSE----
if(!haveMPMmaps){
 b1prefix <- "anon_s2018-02-28_18-26-184837-00001-"
 B1mapRegName <- paste0(b1prefix, "00001-1_B1mapReg.nii")
 B1mapRegFile <- file.path(MPMresDir, pdDir,
                           "Results", "Supplementary",
                           B1mapRegName)
 qMRImapsB1C <- calculateQI(modelMPM,
                           b1File = B1mapRegFile,
                           TR2 = 3.4)
}

## ----label="Figure_6_5", echo=2:9, eval=TRUE, fig.cap="Bias-corrected qMRI maps.", fig.width=7, fig.height=6.5----
oldpar <- par(mfrow=c(2, 2), mar=c(3,3,3,1), mgp=c(2,1,0))
R1 <- extract(qMRImapsB1C, "R1")[ix, iy, iz]
rimage(ix, iz, R1, zlim = c(0, 1.75), main = "R1")
R2star <- extract(qMRImapsB1C, "R2star")[ix, iy, iz]
rimage(ix, iz, R2star, zlim = c(0, 50), main = "R2*")
PD <- extract(qMRImapsB1C, "PD")[ix, iy, iz]
rimage(ix, iz, PD, zlim = c(0, 15000), main = "A")
MT <- extract(qMRImapsB1C, "MT")[ix, iy, iz]
rimage(ix, iz, MT, zlim = c(0, 2), main = "MT", 
       col = colMT)
par(oldpar)


## ----label="estimateNoiseInMPM", echo=TRUE, eval=!havemodelMPM, results=FALSE----
if(!havemodelMPM){
 sigma <- aws::awsLocalSigma(mpm$ddata[1, , , ], 16,
               mpm$mask, ncoils = 1, vext = c(1, 1, 1))
 sigma <- mean(sigma$sigma[mpm$mask])

## ----label="estimateESTATICSQL", echo=TRUE, eval=!havemodelMPM----
 modelMPMQL <- estimateESTATICS(mpm, method = "QL",
                                sigma = sigma, L = 1)

## ----label="calculateMapsQL", echo=TRUE, eval=!havemodelMPM, results=FALSE----
 qMRImapsQL <- calculateQI(modelMPMQL,
                           b1File = B1mapRegFile,
                           TR2 = 3.4)
 save(sigma, modelMPM, modelMPMQL, file=file.path(MPMresDir, "modelMPM.rsc"))
}

## ----label="Parameter extraction", echo=TRUE, eval=TRUE----
mask <- extract(qMRImapsB1C, "mask")[ , iy, ]
R2star <- extract(qMRImapsB1C, "R2star")[ , iy, ]
R2starQL <- extract(qMRImapsQLB1C, "R2star")[ , iy, ]
rR2star <- R2starQL / R2star
x <- R2star[mask]
y <- rR2star[mask]
ind <- is.finite(y) & is.finite(x)
xR2s <- x[ind]; yR2s <- y[ind]

R1 <- extract(qMRImapsB1C, "R1")[, iy, ]
R1QL <- extract(qMRImapsQLB1C, "R1")[, iy, ]
rR1 <- R1QL / R1
x <- R1[mask]
y <- rR1[mask]
ind <- is.finite(y) & is.finite(x)
xR1 <- x[ind]
yR1 <- y[ind]

PD <- extract(qMRImapsB1C, "PD")[, iy, ]
PDQL <- extract(qMRImapsQLB1C, "PD")[, iy, ]
rPD <- PDQL / PD
x <- PD[mask] / mean(PD[mask]) * 69
y <- rPD[mask]
ind <- is.finite(y) & is.finite(x)
xPD <- x[ind]
yPD <- y[ind]

MT <- extract(qMRImapsB1C, "MT")[, iy, ]
MTQL <- extract(qMRImapsQLB1C, "MT")[, iy, ]
rMT <- MTQL / MT
x <- MT[mask] / mean(MT[mask])
y <- rMT[mask]
ind <- is.finite(y) & is.finite(x) & x>0 & x<2 & y>0
xMT <- x[ind]
yMT <- y[ind]

## ----label="Figure_6_6", echo=-1, eval=TRUE, fig.cap="Relative effect of bias correction on the estimated parameters.", fig.width=8, fig.height=4----
oldpar <- par(mfrow = c(1, 4), mar = c(3, 3, 3, 1), mgp = c(2, 1, 0))
plot(c(0, 50), c(0.8, 1.2), xlab = "R2* (1/s)",
     ylab = "Ratio of QL to LS estimate", type = "n")
points(xR2s, yR2s)
lines(c(0, 50), c(1, 1), col = 2)
title("Correction effect on R2*")
meanR2star <- locpoly(xR2s, yR2s, bandwidth = 1.25)
lines(meanR2star, col = 3, lwd = 2)

plot(c(0, 2), c(0.9, 1.1), xlab = "R1 (1/s)",
     ylab = "Ratio of QL to LS estimate", type = "n")
points(xR1, yR1)
lines(c(0, 2), c(1, 1), col = 2)
title("Correction effect on R1")
meanR1 <- locpoly(xR1, yR1, bandwidth = 0.05)
lines(meanR1, col = 3, lwd = 2)

plot(c(0, 100), c(0.9, 1.1), xlab = "PD (p.u.)",
     ylab = "Ratio of QL to LS estimate", type = "n")
points(xPD, yPD)
lines(c(0, 100), c(1, 1), col = 2)
title("Correction effect on PD")
meanPD <- locpoly(xPD, yPD, bandwidth = 2.5)
lines(meanPD, col = 3, lwd = 2)

plot(c(0, 2), c(0.925, 1.4), xlab = "MT (none)",
     ylab = "Ratio of QL to LS estimate", type = "n")
points(xMT, yMT)
lines(c(0, 2.), c(1, 1), col = 2)
title("Correction effect on MT")
meanMT <- locpoly(xMT, yMT, bandwidth = .05)
lines(meanMT, col = 3, lwd = 2)

## ----label="reset2",echo=FALSE,eval=TRUE--------
par(oldpar)

## ----label="Figure_6_7", echo= TRUE, eval=TRUE, fig.cap="Spatial distribution of the relative effect of the bias correction.", fig.width=7, fig.height=7----
rimage.options(xlab = "", ylab = "", yaxt = "n")
layout(matrix(1:8, 4, 2),
       heights = c(0.425, 0.075, 0.425, 0.075))
par(mar = c(.1, 1, 3.1, 1), mgp = c(2, 1, 0) )
img <- rR2star[ix,iz]
rimage(ix,iz,img, zlim = c(0.98, 1.075), xaxt = "n",
       main = list("Correction factor on R2*"))
par(mar = c( 3.1, 1, .1, 1), mgp = c(2, 1, 0) )
rimage(seq(0.98, 1.075, length = 256), 1,
       matrix(1:256, 256, 1))
par(mar = c(0.1, 1, 3.1, 1), mgp = c(2, 1, 0) )
img <- rR1[ix,iz]
rimage(ix,iz,img, zlim = c(0.94, 1.005), xaxt = "n",
       main = list("Correction factor on R1"))
par(mar = c(3.1, 1, .1, 1), mgp = c(2, 1, 0) )
rimage(seq(0.94, 1.005, length = 256), 1,
       matrix(1:256, 256, 1))
par(mar = c(0.1, 1, 3.1, 1), mgp = c(2, 1, 0) )
img <- rPD[ix,iz]
rimage(ix,iz,img, zlim = c(0.95, 1.03), xaxt = "n",
       main = list("Correction factor on PD"))
par(mar = c(3.1, 1, .1, 1), mgp = c(2, 1, 0) )
rimage(seq(0.95, 1.03, length = 256), 1,
       matrix(1:256, 256, 1))
par(mar = c(0.1, 1, 3.1, 1), mgp = c(2, 1, 0) )
img <- rMT[ix,iz]
rimage(ix,iz,img, zlim = c(0.97, 1.075), xaxt = "n",
       main = list("Correction factor on MT"))
par(mar = c(3.1, 1, .1, 1), mgp = c(2, 1, 0) )
rimage(seq(0.97, 1.075, length = 256), 1,
       matrix(1:256, 256, 1))
rimage.options(xlab="x", ylab="z", yaxt="s")

## ----label="havesmoothModel", echo= FALSE, eval=havesmoothedModelMPM----
if(havesmoothedModelMPM){
load(file.path(MPMresDir,"smoothModelMPM.rsc"))
} else {
  
## ----label="smoothModel", echo= TRUE, eval=!havesmoothedModelMPM, results=FALSE----
 modelMPMs <- smoothESTATICS(modelMPM,
                             kstar = 20,
                             alpha = 0.5,
                             verbose = TRUE)
}

## ----label="calculateSmoothedMapsFTcorr", echo=TRUE, eval=!haveMPMmaps, results=FALSE----
if(!haveMPMmaps){
 B1mapRegName <- paste0(b1prefix,
                       "00001-1_B1mapReg.nii")
 B1mapRegFile <- file.path(MPMdataDir, pdDir, "Results",
                           "Supplementary", B1mapRegName)
 qMRISmoothedMaps <- calculateQI(modelMPMs, b1File =
                              B1mapRegFile, TR2 = 3.4)
}

## ----label="Figure_6_8", echo=2:9, eval=TRUE, fig.cap="Smoothed and bias-corrected qMRI maps. Here, the statistical penalty from Eq.~\\eqref{eq:QMRIStatPen} was used.", fig.width=7, fig.height=6.4----
oldpar <- par(mfrow=c(2, 2), mar=c(3,3,3,1), mgp=c(2,1,0))
R1 <- extract(qMRISmoothedMaps, "R1")[ix, iy, iz]
rimage(ix, iz, R1, zlim = c(0, 1.75), main = "R1")
R2star <- extract(qMRISmoothedMaps,"R2star")[ix, iy, iz]
rimage(ix, iz, R2star, zlim = c(0, 50), main = "R2*")
PD <- extract(qMRISmoothedMaps, "PD")[ix, iy, iz]
rimage(ix, iz, PD, zlim = c(0, 15000), main = "PD")
MT <- extract(qMRISmoothedMaps, "MT")[ix, iy, iz]
rimage(ix, iz, MT, zlim = c(0, 2), main = "MT", 
       col = colMT)
par(oldpar)


## ----label="smoothModelp1", echo=TRUE, eval=!havesmoothedModelMPM, results=FALSE----
if(!havesmoothedModelMPM){
 modelMPMsp1 <- smoothESTATICS(modelMPM,
                               kstar = 20,
                               alpha = 0.5,
                               patchsize=1,
                               verbose = TRUE)
 save(modelMPMs, modelMPMsp1, file=file.path(MPMresDir, "smoothModelMPM.rsc"))
}

## ----label="calculateSmoothedp1MapsFTcorr", echo=TRUE, eval=!haveMPMmaps, results=FALSE----
if(!haveMPMmaps){
 B1mapRegName <- paste0(b1prefix,
                        "00001-1_B1mapReg.nii")
 B1mapRegFile <- file.path(dataDir, "MPM", pdDir,
                           "Results", "Supplementary",
                           B1mapRegName)
 qMRISmoothedp1Maps <- calculateQI(modelMPMsp1,
                                   b1File = B1mapRegFile,
                                   TR2 = 3.4)
 save(qMRImaps, qMRImapsB1C, qMRImapsQLB1C, qMRISmoothedMaps, qMRISmoothedp1Maps, file=file.path(MPMresDir, "qMaps.rsc"))
}

## ----label="Figure_6_9", echo=2:9, eval=TRUE, fig.cap="Smoothed and bias-corrected qMRI maps. Here the statistical penalty from Eq.~\\eqref{eq:QMRIStatPenPAWS} of the patch wise refinement of the adaptive smoothing procedure was used.", fig.width=7, fig.height=6.4----
oldpar <- par(mfrow=c(2, 2), mar=c(3,3,3,1), mgp=c(2,1,0))
R1 <- extract(qMRISmoothedp1Maps, "R1")[ix, iy, iz]
rimage(ix, iz, R1, zlim = c(0, 1.75), main = "R1")
R2star <-
  extract(qMRISmoothedp1Maps, "R2star")[ix, iy, iz]
rimage(ix, iz, R2star,  zlim = c(0, 50), main = "R2*")
PD <- extract(qMRISmoothedp1Maps, "PD")[ix, iy, iz]
rimage(ix, iz, PD, zlim = c(0, 15000), main = "PD")
MT <- extract(qMRISmoothedp1Maps, "MT")[ix, iy, iz]
rimage(ix, iz, MT, zlim = c(0, 2), main = "MT", 
       col = colMT)
par(oldpar)


## ----label="Figure_6_10", echo=2:22, eval=TRUE, fig.cap="Comparison of un-smoothed and smoothed quantitative maps. The adaptive smoothing used the statistical penalty from Eq.~\\eqref{eq:QMRIStatPenPAWS}. The procedure is able to reduce the noise in the maps without blurring fine anatomical structures.", fig.width=7, fig.height=4----
oldpar <- par(mfrow=c(2, 4), mar=c(3,3,3,.1), mgp=c(2,1,0))
ix <- 151:200; iz <- 121:165
R1 <- extract(qMRImaps, "R1")[ix, iy, iz]
rimage(ix, iz, R1, zlim = c(0, 1.75), main="R1")
R2star <- extract(qMRImaps, "R2star")[ix, iy, iz]
rimage(ix, iz, R2star, zlim = c(0, 50), main="R2*")
PD <- extract(qMRImaps, "PD")[ix, iy, iz]
rimage(ix, iz, PD, zlim = c(3000, 11000), main="PD")
MT <- extract(qMRImaps, "MT")[ix, iy, iz]
rimage(ix, iz, MT, zlim = c(0, 2), main="MT", 
       col = colMT)
R1 <- extract(qMRISmoothedp1Maps, "R1")[ix, iy, iz]
rimage(ix, iz, R1, zlim = c(0, 1.75),
       main = "smoothed R1")
R2star <- extract(qMRISmoothedp1Maps, "R2star")[ix,iy,iz]
rimage(ix, iz, R2star, zlim = c(0, 50),
       main = "smoothed R2*")
PD <- extract(qMRISmoothedp1Maps, "PD")[ix, iy, iz]
rimage(ix, iz, PD, zlim = c(3000, 11000),
       main = "smoothed PD")
MT <- extract(qMRISmoothedp1Maps, "MT")[ix, iy, iz]
rimage(ix, iz, MT, zlim = c(0, 2), col = colMT,
       main = "smoothed MT")

