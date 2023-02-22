## ----label="initnpsm", echo=FALSE, eval=TRUE----------------------------------------------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(cache = FALSE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE)
options(width=50)


## ----label="initnpara", echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------------------------------
library(oro.nifti)
library(aws)
library(adimpro)
rimage.options(ylab = "z",
               zquantile = c(0.001, 0.999))
setCores(4)

baseDir <- "../.."
codeDirNP <- file.path(baseDir,"MRIwithR","Chapter_nonpara")
source(file.path(codeDirNP,"chapter_nonpara_init.R"))


## ----"precomputed", eval=TRUE, echo=FALSE-------------------------------------------------------------------------------------------------------------------------
havepropagation <- file.exists(file.path(npresDir, "propagation.rsc"))
havesmoothed <- file.exists(file.path(npresDir, "awsestimate.rsc"))


## ----label="MPMT1data", echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------------------------------
t1Name <- file.path(dataDir, "MPM",
    "t1w_mfc_3dflash_v1i_R4_0015",
    "anon_s2018-02-28_18-26-190921-00001-00224-1.nii")
T1 <- as.array(readNIfTI(t1Name, reorient = FALSE))


## ----label="smoothT1", echo=TRUE, eval=TRUE-----------------------------------------------------------------------------------------------------------------------
T1sm <- kernsm(T1, h = 2., kern = "Epanechnikov")@yhat


## ----label="Figure_A_1", echo= -1, eval=TRUE, fig.cap="Original and smoothed $T_1$ image.", fig.width=10, fig.height=4.55-----------------------------------------
par(mfrow=c(1,2), mar=c(3,3,3,1), mgp=c(2,1,0))
rimage(T1[,160,], main="Original",
       zlim= c(0,2600))
rimage(T1sm[,160,], main="Smoothed image",
       zlim= c(0,2600))


## ----label="Figure_A_2", echo= -1, eval=TRUE, fig.cap="Effect of bandwidth for a detail of the $T_1$ map.", fig.width=7.5, fig.height=5.33, cache=TRUE------------
par(mfrow=c(2,3), mar=c(3,3,3,.1), mgp=c(2,1,0))
ix <- 151:200; iy <- 160; iz <- 121:165
for (h in c(1, 1.15, 1.7, 2.7, 4.2, 6.6)){
  T1sm <- kernsm(T1, h, kern = "Epanechnikov")@yhat
  rimage(ix, iz, T1sm[ix, iy, iz],
         main = paste("bw=", h))
}


## ----"clear",eval=TRUE, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------
rm(T1sm)


## ----label="initadimpro", echo=FALSE, eval=TRUE, message=FALSE----------------------------------------------------------------------------------------------------
rimage.options(ylab = "z",
               zquantile = c(0.001, 0.999))


## ----"load propagation", echo=FALSE, eval=havepropagation---------------------------------------------------------------------------------------------------------
if(havepropagation){
 load(file.path(npresDir, "propagation.rsc"))
} else {

## ----"Propagationcondition",echo=TRUE,eval=!havepropagation,results=FALSE,fig.keep="none"-------------------------------------------------------------------------
 z1 <- awstestprop(c(200, 200, 200), hmax = 6)
 z2 <- awstestprop(c(200, 200, 200), hmax = 6,
                   ladjust = 1.5)
}

## ----"Figure_A_3",echo=-1,eval=TRUE,results=FALSE,fig.cap="Contour plots to evaluate the propagation condition.",fig.width=10,fig.height=5------------------------
par(mfrow=c(1,2),mar=c(3,3,3,.5),mgp=c(2,1,0))
step <- 1:length(z1$h)
contour(z1$z, step, z1$prob, levels = z1$levels,
        xlab = "z", ylab = "step", lwd = 2)
title("Exceedence probabilities ladjust=1")
contour(z1$z, step, z1$probna, levels = z1$levels,
        add = TRUE, col = 2, lty = 2)
contour(z2$z, step, z2$prob, levels = z2$levels,
        xlab = "z", ylab = "step", lwd = 2)
title("Exceedence probabilities ladjust=1.5")
contour(z2$z, step, z2$probna, levels = z2$levels,
        add = TRUE, col = 2, lty = 2)


## ----"rm propagation", eval=TRUE, echo=FALSE----------------------------------------------------------------------------------------------------------------------
rm(z1, z2)


## ----"load smoothed", echo=FALSE, eval=havesmoothed---------------------------------------------------------------------------------------------------------------
if(havesmoothed){
 load(file.path(npresDir, "awsestimate.rsc"))
} else {

## ----"estimating sigma",echo=TRUE,eval=FALSE,results=!havesmoothed,fig.keep="none"--------------------------------------------------------------------------------
 sigma2 <- mean(aws.gaussian(
    T1[121:220, 121:200, 121:200], 10)@sigma2)

## ----"AWSestimates",echo=TRUE,eval=FALSE,results=!havesmoothed,fig.keep="none"------------------------------------------------------------------------------------
 T1aws6 <- aws(T1, hmax = 1.15, sigma2 = sigma2)@theta
 T1aws12 <- aws(T1, hmax = 1.7, sigma2 = sigma2)@theta
 T1aws18 <- aws(T1, hmax = 2.7, sigma2 = sigma2)@theta
 T1aws24 <- aws(T1, hmax = 4.2, sigma2 = sigma2)@theta
 T1aws30 <- aws(T1, hmax = 6.6, sigma2 = sigma2)@theta
}

## ----label="Figure_A_4", echo= -(1:2), eval=TRUE, fig.cap="Smoothing using AWS: results for a detail of the $T_1$ map.", fig.width=7.5, fig.height=5.33-----------
par(mfrow=c(2,3), mar=c(3,3,3,.1), mgp=c(2,1,0))
ix <- 151:200; iy <- 160; iz <- 121:165
rimage(ix, iz, T1[ix, iy, iz],
       main = "original")
rimage(ix, iz, T1aws6[ix, iy, iz],
       main = "steps = 6")
rimage(ix, iz, T1aws12[ix, iy, iz],
       main = "steps = 12")
rimage(ix, iz, T1aws18[ix, iy, iz],
       main = "steps = 18")
rimage(ix, iz, T1aws24[ix, iy, iz],
       main = "steps = 24")
rimage(ix, iz, T1aws30[ix, iy, iz],
       main = "steps = 30")

## ----"PAWSestimates",echo=-(1:2),eval=FALSE,results=!havesmoothed,fig.keep="none"---------------------------------------------------------------------------------
if(!havesmoothed){
 T1paws6 <- paws(T1, hmax = 1.15, sigma2 = sigma2)@theta
 T1paws12 <- paws(T1, hmax = 1.7, sigma2 = sigma2)@theta
 T1paws18 <- paws(T1, hmax = 2.7, sigma2 = sigma2)@theta
 T1paws24 <- paws(T1, hmax = 4.2, sigma2 = sigma2)@theta
 T1paws30 <- paws(T1, hmax = 6.6, sigma2 = sigma2)@theta
}

## ----label="Figure_A_5", echo= -(1:2), eval=TRUE, fig.cap="Smoothing using patchwise AWS: results for a detail of the $T_1$ map.", fig.width=7.5, fig.height=5.33----
par(mfrow=c(2,3), mar=c(3,3,3,.1), mgp=c(2,1,0))
ix <- 151:200; iy <- 160; iz <- 121:165
rimage(ix, iz, T1[ix, iy, iz],
       main = "original")
rimage(ix, iz, T1paws6[ix, iy, iz],
       main = "steps = 6")
rimage(ix, iz, T1paws12[ix, iy, iz],
       main = "steps = 12")
rimage(ix, iz, T1paws18[ix, iy, iz],
       main = "steps = 18")
rimage(ix, iz, T1paws24[ix, iy, iz],
       main = "steps = 24")
rimage(ix, iz, T1paws30[ix, iy, iz],
       main = "steps = 30")

