IRresDir <- "../results/IR"
if(!dir.exists(IRresDir)) dir.create(IRresDir)
IRT1Dir <- file.path(IRresDir,"T1")
if(!exists(IRT1Dir)) dir.create(file.path(IRresDir,"T1"))
library(kirby21.t1)
download_t1_data()
library(fslr)
betFile <- file.path(IRresDir,"betFileKirby21T1.nii.gz")
fslbet(get_t1_filenames()[6],outfile=betFile,retimg=TRUE,reorient=FALSE,betcmd="bet2",opt="-m")
orthographic(readNIfTI("../results/IR/betFileKirby21T1.nii.gz",reorient=FALSE))
fastDir <- file.path(IRresDir,"FileKirby21T1Fast")
fast(betFile, outfile=fastDir)
segimg <- readNIfTI("../results/IR/FileKirby21T1Fast_pveseg.nii.gz",reorient=FALSE)[30:142,49:213,88:184]
writeNIfTI(segimg,file.path(IRresDir,"SegmentedKirbyT1"))

segimg <- readNIfTI(file.path(IRresDir,"SegmentedKirbyT1.nii.gz"),reorient=FALSE)[30:142,49:213,88:184]

Sf <- 900
Rf <- 0.000285
Sgm <- 400
Rgm <- 0.00075
fgm <- .15
Swm <- 370
Rwm <- 0.0011
fwm <- .05
InvTimes <- c(100, 200, 400, 600, 800, 1200, 1600, 2000, 2500, 3000, 
              3500, 4000, 4500, 5000, 6000, Inf)
InvTimes0 <- c(100, 200, 400, 600, 800, 1200, 1600, 2000, 2500, 3000, 
               3500, 4000, 4500, 5000, 6000, 15000)

fintCSF <- qMRI:::IRhomogen(c(Sf,Rf),InvTimes0)
fintGM <- qMRI:::IRmix2(c(fgm,Rgm,Sgm),InvTimes0,Sf,Rf)
fintWM <- qMRI:::IRmix2(c(fwm,Rwm,Swm),InvTimes0,Sf,Rf)

set.seed(1)
sigma <- 40
nTimes <- length(InvTimes0)
nCSF <- sum(segimg==1)
nGM <- sum(segimg==2)
nWM <- sum(segimg==3)
IRdata <- array(0,c(nTimes,prod(dim(segimg))))
IRdata[,segimg==1] <- sqrt(rnorm(nTimes*nCSF,fintCSF,sigma)^2+
                           rnorm(nTimes*nCSF,0,sigma)^2)
IRdata[,segimg==2] <- sqrt(rnorm(nTimes*nGM,fintGM,sigma)^2+
                           rnorm(nTimes*nGM,0,sigma)^2)
IRdata[,segimg==3] <- sqrt(rnorm(nTimes*nWM,fintWM,sigma)^2+
                           rnorm(nTimes*nWM,0,sigma)^2)
dim(IRdata) <- c(nTimes,dim(segimg))
for(i in 1:9) writeNIfTI(as.nifti(IRdata[i,,,]), 
                         file.path(IRT1Dir,paste0("IR0",i)))
for(i in 10:nTimes) writeNIfTI(as.nifti(IRdata[i,,,]), 
                               file.path(IRT1Dir,paste0("IR",i)))

