if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_fmri","chapter_fmri_init.R"))

## ----label="Gaussian smoothing", echo=TRUE, eval=TRUE, results=FALSE, cache=TRUE----------------------------
nii <- snii <- readNIfTI(f105s1r1mc, reorient = FALSE)
bw <- 8/nii@pixdim[2:4]
for (i in 1:dim(nii)[4])
   snii@.Data[ , , , i] <- kernsm(nii[, , , i],
                              h = bw, unit="FWHM")@yhat

## ----label="create linear model1", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE-------------------------
nii <- readNIfTI(f105s1r1mc, reorient = FALSE)
TR <- 2.5
scans <- dim(nii)[4]


## ----label="create linear model2", echo=TRUE, eval=TRUE, results=FALSE,cache=FALSE--------------------------
ttt <- read.table(f105s1r1events, header = TRUE)
ntrials <- dim(ttt)[1]
indScissors <- (1:ntrials)[ttt$trial_type == "scissors"]
indFace <- (1:ntrials)[ttt$trial_type == "face"]
indCat <- (1:ntrials)[ttt$trial_type == "cat"]
indShoe <- (1:ntrials)[ttt$trial_type == "shoe"]
indHouse <- (1:ntrials)[ttt$trial_type == "house"]
indScrb <- (1:ntrials)[ttt$trial_type == "scrambledpix"]
indBottle <- (1:ntrials)[ttt$trial_type == "bottle"]
indChair <- (1:ntrials)[ttt$trial_type == "chair"]
onsets <- ttt$onset
duration <- ttt$duration


## ----label="create linear model3", echo=TRUE, eval=TRUE, results=FALSE--------------------------------------
hrfScissors <- fmri.stimulus(scans, onsets[indScissors],
                             duration[indScissors],
                             TR = TR, times = TRUE)
hrfFace <- fmri.stimulus(scans, onsets[indFace],
                         duration[indFace],
                         TR = TR, times = TRUE)
hrfCat <- fmri.stimulus(scans, onsets[indCat],
                        duration[indCat],
                        TR = TR, times=TRUE)
hrfShoe <- fmri.stimulus(scans, onsets[indShoe],
                         duration[indShoe],
                         TR = TR, times = TRUE)
hrfHouse <- fmri.stimulus(scans, onsets[indHouse],
                          duration[indHouse],
                          TR = TR, times = TRUE)
hrfScrambled <- fmri.stimulus(scans, onsets[indScrb],
                              duration[indScrb],
                              TR = TR, times = TRUE)
hrfBottle <- fmri.stimulus(scans, onsets[indBottle],
                           duration[indBottle],
                           TR = TR, times = TRUE)
hrfChair <- fmri.stimulus(scans, onsets[indChair],
                          duration[indChair],
                          TR = TR, times = TRUE)
hrf <- cbind(hrfScissors,hrfFace,hrfCat,hrfShoe,
             hrfHouse,hrfScrambled,hrfBottle,hrfChair)


## ----label="create linear model4", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE-------------------------
load(f105s1r1motion)
xdesign <- fmri.design(hrf, order = 2, cef = mocoparams)

## ----label="evaluate GLM", echo=TRUE, eval=TRUE, results=FALSE----------------------------------------------
nii <- readNIfTI(f105s1r1mc, reorient = FALSE)
mask <- readNIfTI(f105s1boldmask) > 0
ds <- oro2fmri(nii)
ds$mask <- mask
contrast <- rep(1, 8) # Image versus rest
spm <- fmri.lm(ds, xdesign, contrast = contrast,
               mask=mask)

## ----label="Figure_4_9", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Voxelwise signal detection.", fig.width=10,fig.height=3.3----
pvalue.v <- fmri.pvalue(spm, "voxelwise")

## ----label="Figure_4_10", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Voxelwise signal detection with Bonferrroni correction.", fig.width=10,fig.height=3.3----
pvalue.b <- fmri.pvalue(spm, "Bonferroni")


## ----label="smooth and reevaluate GLM", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, warning=FALSE-----
# thats 8mm spatial FWHM bandwidth
ds.s <- oro2fmri(snii)
ds.s$mask <- mask
spm.s <- fmri.lm(ds.s, xdesign, contrast = contrast, mask=mask)


## ----label="Figure_4_11", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE, fig.cap="Signal detection by RFT.", fig.width=10,fig.height=3.3----
pvaluesm <- fmri.pvalue(spm.s)




## ----label="Figure_4_13", echo=TRUE, eval=TRUE, results=FALSE, fig.cap="Voxelwise signal detection by FDR.", fig.width=10,fig.height=3.3----
pvaluefdr <- fmri.pvalue(spm, mode = "FDR", alpha = .01)


## ----label="p-value with cluster thresholding", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE------------
pvaluecl <-fmri.cluster(spm, alpha =.01,
                        ncmin =8, ncmax=15)




## ----label="adaptive smoothing with fmri", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE-----------------
spm.smooth <- fmri.smooth(spm, hmax = 3,
                          adaptation = "awsfull")


## ----label="p-values with aws fmri", echo=TRUE, eval=TRUE, results=FALSE------------------------------------
pvalueaws <- fmri.pvalue(spm.smooth, alpha=.01)

## ----label="segment with fmri", echo=TRUE, eval=TRUE, results=FALSE, cache=FALSE----------------------------
spm.segment <- fmri.smooth(spm, hmax = 3,
                  adaptation = "segment", alpha = 0.01)



## ----label="MVPA with fmri", echo=TRUE, eval=TRUE, results=FALSE--------------------------------------------
nii <- readNIfTI(f105s1r1mc, reorient = FALSE)
ds <- oro2fmri(nii)
ds$mask <- mask
cfaces <- c(-1/7, 1, rep(-1/7,6))
chouses <- c(rep(-1/7, 4), 1, rep(-1/7, 3))
spm.faces <- fmri.lm(ds, xdesign, contrast = cfaces, mask=mask)
spm.houses <- fmri.lm(ds, xdesign, contrast = chouses, mask=mask)
pvaluesl.faces <- fmri.searchlight(spm.faces,
                                   alpha = .01,
                                   radius = 2,
                                   kind = "abs")
pvaluesl.houses <- fmri.searchlight(spm.houses,
                                    alpha = .01,
                                    radius = 2,
                                    kind = "abs")


pattern.faces <- getSearchlightPattern(spm.faces,
                      pvaluesl.faces$pvalue < 0.05,
                      radius = 2)
dim(pattern.faces)

save(nii, snii, scans, xdesign, spm, pattern.faces, pvaluesl.faces, 
     pvaluesl.houses, spm.segment, spm.smooth, pvalueaws, pvaluecl, pvaluefdr, pvaluesm, 
     pvalue.v, pvalue.b, file=file.path(resDir,"fMRI","spms.rsc"))
