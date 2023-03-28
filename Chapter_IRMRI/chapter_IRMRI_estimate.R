if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_IRMRI","chapter_IRMRI_init.R"))

library(qMRI)
setCores(16)
t1Files <- list.files(IRT1Dir,"*.nii.gz",full.names=TRUE)
segmFile <- file.path(IRresDir,"SegmentedKirbyT1.nii.gz")
if(!file.exists(segmFile)) source(file.path(baseDir,"MRIwithR","Chapter_IRMRI","chapter_IRMRI_prepare.R"))
InvTimes <- c(100, 200, 400, 600, 800, 1200, 1600, 2000, 2500, 3000, 
               3500, 4000, 4500, 5000, 6000, 15000)
IRdata <- readIRData(t1Files, InvTimes, segmFile, sigma=sigma,
                     L=1, segmCodes=c("CSF","GM","WM"))
IRfluid <- estimateIRfluid(IRdata, method="NLR", verbose=FALSE)
cat("Estimated parameters Sf:", IRfluid$Sf, 
    " Rf:", IRfluid$Rf, "\n")
IRmix <- estimateIRsolid(IRfluid, verbose=FALSE)
alpha <- 1e-4
sIRmix <- smoothIRSolid(IRmix, alpha=alpha, partial=FALSE, verbose=FALSE)

sIRmix <- estimateIRsolidfixed(sIRmix, verbose=FALSE)
save(IRfluid, IRmix, sIRmix, file=file.path(IRresDir,"IRestimates.rsc"))


