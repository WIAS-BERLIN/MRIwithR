if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_fmri","chapter_fMRI_init.R"))

## ----label="slice timing version 1", echo=TRUE, eval=FALSE, warning=FALSE-----------------------------------
niirun1.o <- readNIfTI(f117s1r1, reorient = FALSE)
niirun1 <- oro2fmri(niirun1.o, setmask = FALSE)
sltimes <- c(seq(1, 33, 2), seq(2, 32, 2))
stcniirun1 <- slicetiming(niirun1, sltimes)
stcdata <- extractData(stcniirun1)
stcdata.nii <- as.nifti(stcdata, niirun1.o)
writeNIfTI(stcdata.nii, file = f117s1r1stc)

