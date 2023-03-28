## ----label="initQMRI", echo=FALSE, eval=TRUE-----------
library(knitr)
knitr::opts_chunk$set(warning = FALSE,
                      message = FALSE)
options(width=50, digits=3)


## ----label="packages and filenames", eval=TRUE, echo=-1----------
if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_IRMRI","chapter_IRMRI_init.R"))

## ----label="initIR", echo=FALSE, eval=TRUE, message=FALSE, results=FALSE-----------
rimage.options(zquantiles=c(0.001,0.999), 
               xlab="x", ylab="z", bty="n")
setCores(4) # number of cores for parallel computations


## ----"IRhaves", eval=TRUE, echo=FALSE----------
haveKirbySegm <- file.exists(file.path( 
                IRresDir,"SegmentedKirbyT1.nii.gz"))
haveIRFiles <- file.exists(file.path(IRT1Dir,
                                     "IR16.nii.gz"))
haveIRestimates <- file.exists(file.path(IRresDir,
                                 "IRestimates.rsc"))


## ----"segment kirbyT1", echo=TRUE, eval=!haveKirbySegm----------
if(!haveKirbySegm){
download_t1_data()
brainFile <- file.path(IRresDir, "betFileKirby21T1.nii.gz")
fslbet(get_t1_filenames()[6], outfile=brainFile, 
       retimg=TRUE, reorient=FALSE, betcmd="bet2",
       opt="-m")
fastDir <- file.path(IRresDir,"Kirby21T1Fast")
fast(brainFile, outfile=fastDir)
segimg <- readNIfTI(file.path(IRresDir, 
                              "Kirby21T1Fast_pveseg.nii.gz"),
                    reorient=FALSE)[30:142,49:213,88:184]
writeNIfTI(segimg, file.path(IRresDir,
                             "SegmentedKirbyT1"))
}

## ----"set IR parameters ", echo=TRUE, eval=TRUE----------

Sf <- 900
Rf <- 0.000285
Sgm <- 400
Rgm <- 0.00075
fgm <- .15
Swm <- 370
Rwm <- 0.0011
fwm <- .05
InvTimes <- c(100, 200, 400, 600, 800, 1200, 1600, 
              2000, 2500, 3000, 3500, 4000, 4500, 
              5000, 6000, 15000)




## ----label="Fig:IRcurves", fig.width=12, fig.height=6, fig.cap="Intensities as functions of inversion times and tissue type (black for CSF, red for GM and green for WM)"----
x <- seq(100, 15000, 10)
fintCSF <- qMRI:::IRhomogen(c(Sf, Rf), InvTimes)
fintGM <- qMRI:::IRmix2(c(fgm, Rgm, Sgm), 
                        InvTimes, Sf, Rf)
fintWM <- qMRI:::IRmix2(c(fwm, Rwm, Swm), 
                        InvTimes, Sf, Rf)
plot(InvTimes,fintCSF, xlab="Inversion time", 
                        ylab="Intensity")
points(InvTimes,fintGM,col=2)
points(InvTimes,fintWM,col=3)
lines(x,qMRI:::IRhomogen(c(Sf, Rf), x))
lines(x,qMRI:::IRmix2(c(fgm,Rgm,Sgm), x, Sf, Rf), col=2)
lines(x,qMRI:::IRmix2(c(fwm,Rwm,Swm), x, Sf, Rf), col=3)


## ----"read segimg", eval=TRUE, echo=TRUE-----------
segimg <- readNIfTI( 
  file.path(IRresDir,"SegmentedKirbyT1.nii.gz"),
  reorient=FALSE)

## ----"create IR data", eval=!haveIRFiles, echo=TRUE----------
if(!haveIRFiles){  
set.seed(1)
sigma <- 40
nTimes <- length(InvTimes)
nCSF <- sum(segimg==1)
nGM <- sum(segimg==2)
nWM <- sum(segimg==3)
IRdata <- array(0,c(nTimes,prod(dim(segimg))))
IRdata[,segimg==1] <- sqrt(rnorm(nTimes*nCSF, fintCSF,
                                 sigma)^2 + rnorm(nTimes*nCSF, 0, sigma)^2)
IRdata[,segimg==2] <- sqrt(rnorm(nTimes*nGM, fintGM,
                                 sigma)^2 + rnorm(nTimes*nGM, 0, sigma)^2)
IRdata[,segimg==3] <- sqrt(rnorm(nTimes*nWM, fintWM,
                                 sigma)^2 + rnorm(nTimes*nWM, 0, sigma)^2)
dim(IRdata) <- c(nTimes,dim(segimg))
for(i in 1:nTimes) writeNIfTI(as.nifti(IRdata[i,,,]), 
                              file.path(IRT1Dir,paste0("IR",str_pad(i,2,pad="0"))))
}

## ----"Read IR data", echo=TRUE, eval=TRUE----------
IRFiles <- list.files(IRT1Dir, "*.nii.gz",
                      full.names=TRUE)
segmFile <- file.path(IRresDir,
                      "SegmentedKirbyT1")
IRdata <- readIRData(IRFiles, InvTimes, segmFile, 
        sigma=sigma, L=1, segmCodes=c("CSF","GM","WM"))

## ----"Load IR estimates", echo=FALSE, eval=haveIRestimates----------
if(haveIRestimates){
  load(file.path(IRresDir,"IRestimates.rsc"))
cat("Estimated parameters Sf:", IRfluid$Sf, 
    " Rf:", IRfluid$Rf, "\n")
} else {
  
  ## ----"Estimate Fluid parameters", echo=TRUE, eval=!haveIRestimates ----------
  IRfluid <- estimateIRfluid(IRdata, method="NLR", 
                             verbose=FALSE)
}
  cat("Estimated parameters Sf:", IRfluid$Sf, 
    " Rf:", IRfluid$Rf, "\n")

## ----"Estimate solid parameters", echo=TRUE, eval=!haveIRestimates----------
if(!haveIRestimates){
    IRmix <- estimateIRsolid(IRfluid, verbose=FALSE)
}

## ----"Quality of results", eval=TRUE, echo=TRUE----------
cat("MAE: f(GM)", mean(abs(IRmix$fx[segimg==2]-fgm)),
  "f(WM)", mean(abs(IRmix$fx[segimg==3]-fwm)),"\n",
  "    S(GM)", mean(abs(IRmix$Sx[segimg==2]-Sgm)),
  "S(WM)", mean(abs(IRmix$Sx[segimg==3]-Swm)),"\n",
  "    R(GM)", mean(abs(IRmix$Rx[segimg==2]-Rgm)),
  "R(WM)", mean(abs(IRmix$Rx[segimg==3]-Rwm)))

## ----label="Fig_7_1", echo=-1, eval=TRUE, fig.width=10, fig.height=3.78, fig.cap="Estimated IR parameter maps for central slice"-----------
par(mfrow=c(1,4),mar=c(3,3,3,.5),mgp=c(2,1,0))
rimage(segimg[,,48])
title("Segmentation map")
rimage(IRmix$Sx[,,48],zlim=c(250,500))
title("Solid intensity map")
rimage(IRmix$Rx[,,48],zlim=c(0,.0015))
title("Solid relaxation rate map")
rimage(IRmix$fx[,,48],zlim=c(0,.4))
title("Fluid proportion map")

## ----"Smooth solid parameters", echo=TRUE, eval=!haveIRestimates----------
if(!haveIRestimates){
  alpha <- 1e-4
sIRmix <- smoothIRSolid(IRmix, alpha=alpha, 
                        partial=FALSE, verbose=FALSE)
sIRmix <- estimateIRsolidfixed(sIRmix, verbose=FALSE)
save(IRfluid, IRmix, sIRmix, file=file.path(IRresDir,"IRestimates.rsc"))
}

## ----label="Fig_7_2", echo=-1, eval=TRUE, fig.width=10, fig.height=3.78, fig.cap="Smoothed IR parameter maps for central slice"----------
par(mfrow=c(1,4),mar=c(3,3,3,.5),mgp=c(2,1,0))
rimage(segimg[,,48])
title("Segmentation")
rimage(sIRmix$Sx[,,48],zlim=c(250,500))
title("solid intensity map")
rimage(sIRmix$Rx[,,48],zlim=c(0,.0015))
title("solid relaxation rate map")
rimage(sIRmix$fx[,,48],zlim=c(0,.4))
title("fluid proportion map")

## ----"Quality of smoothed results", eval=TRUE, echo=TRUE-----------
cat("MAE: f(GM)", mean(abs(sIRmix$fx[segimg==2]-fgm)),
  "f(WM)", mean(abs(sIRmix$fx[segimg==3]-fwm)),"\n",
  "    S(GM)", mean(abs(sIRmix$Sx[segimg==2]-Sgm)),
  "S(WM)", mean(abs(sIRmix$Sx[segimg==3]-Swm)),"\n",
  "    R(GM)", mean(abs(sIRmix$Rx[segimg==2]-Rgm)),
  "R(WM)", mean(abs(sIRmix$Rx[segimg==3]-Rwm)))

## ----"all IR in one", eval=FALSE, echo=TRUE----------
## sIRmix <- estimateIR(IRdata, method="QL") ## computationally expensive !

