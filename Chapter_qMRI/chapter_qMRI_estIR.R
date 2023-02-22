IRresDir <- "../results/IR"
if(!dir.exists(IRresDir)) dir.create(IRresDir)
IRT1Dir <- file.path(IRresDir,"T1")
if(!dir.exists(IRT1Dir)) dir.create(
                  file.path(IRresDir,"T1"))

library(kirby21.t1)
library(qMRI)
setCores(4)
t1Files <- list.files(IRT1Dir,"*.nii.gz",full.names=TRUE)
segmFile <- file.path(IRresDir,"SegmentedKirbyT1")
IRdata <- readIRData(t1Files, InvTimes0, segmFile, sigma=sigma,
                     L=1, segmCodes=c("CSF","GM","WM"))
IRfluid <- estimateIRfluid(IRdata, method="NLR", verbose=FALSE)
cat("Estimated parameters Sf:", IRfluid$Sf, 
    " Rf:", IRfluid$Rf, "\n")
IRmix <- estimateIRsolid(IRfluid, verbose=FALSE)
alpha <- 1e-4
sIRmix <- smoothIRSolid(IRmix, alpha=alpha, verbose=FALSE)

sIRmix <- estimateIRsolidfixed(sIRmix, verbose=FALSE)
save(IRfluid, IRmix, sIRmix, file=file.path(IRresDir,"IRestimates.rsc"))


