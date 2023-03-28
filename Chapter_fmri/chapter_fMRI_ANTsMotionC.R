if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_fmri","chapter_fMRI_init.R"))

if(!haveANTsR) stop("need package ANTsR installed \n precomputed objects  available in results")
## ----label="Use ANTsR for motion correction1", echo=TRUE, eval=FALSE, results=FALSE-------------------------
fmrirun1 <- antsImageRead(f105s1r1)
fmrirun1MC <- antsrMotionCalculation(fmrirun1,
                    typeofTransform = "BOLDRigid")


## ----label="Use ANTsR for motion correction2", echo=TRUE, eval=FALSE, results=FALSE-------------------------
badtimes <- which(fmrirun1MC$fd$MeanDisplacement > .2 )
## motion corrected time series
antsImageWrite(fmrirun1MC$moco_img, f105s1r1mc)
## average corrected image
if(packageVersion("ANTsR") == '0.6.1'){
  avimg <- as.antsImage(fmrirun1MC$moco_avg_img,
      spacing = antsGetSpacing(fmrirun1)[1:3],
      origin = antsGetOrigin(fmrirun1)[1:3],
      direction = antsGetDirection(fmrirun1)[1:3, 1:3])
} else{
  avimg <- fmrirun1MC$moco_avg_img
}
antsImageWrite(avimg, f105s1r1amc)
## save estimated parameters
mocoparams <- fmrirun1MC$moco_params
save(badtimes, mocoparams, file = f105s1r1motion)
