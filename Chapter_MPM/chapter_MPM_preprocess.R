if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_MPM","chapter_MPM_init.R"))

fslbet(t1Files[1], outfile = maskFile, retimg = TRUE,
       reorient = FALSE, opts = "-f 0.3 -m",
       betcmd = "bet2", verbose = FALSE)

## ----label="RegisterMPMdata1", echo=TRUE, eval=TRUE---------------------------------------------------------
library(ANTsR)
T1 <- antsImageRead(t1Files[1])
MT <- antsImageRead(mtFiles[1])
PD <- antsImageRead(pdFiles[1])
mask <- antsImageRead(paste0(maskFile, ".nii.gz"))
aregMT2T1 <- antsRegistration(fixed = mask,
        moving = MT, typeofTransform = "DenseRigid",
        verbose = TRUE)


## ----label="RegisterMPMdata2", echo=TRUE, eval=TRUE---------------------------------------------------------
for (id in 1:6){
  img <- antsImageRead(mtFiles[id])
  warpedMT <- antsApplyTransforms(fixed = T1,
               moving = img,
               transformlist = aregMT2T1$fwdtransforms)
  antsCopyImageInfo(img, warpedMT)
  antsImageWrite(warpedMT, rmtFiles[id])
}

## ----label="RegisterMPMdata3", echo=TRUE, eval=TRUE---------------------------------------------------------
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

