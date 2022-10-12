## ----label="ANTsR pre-processing for group ICA1", echo=TRUE, eval=FALSE, results=FALSE----------------------
fmrirest <- antsImageRead(frest)
pfmrirest <- preprocessfMRI(fmrirest)
antsImageWrite(pfmrirest$cleanBoldImage,fprest)
antsImageWrite(pfmrirest$maskImage, fprestmask)
avimg <-getAverageOfTimeSeries(pfmrirest$cleanBoldImage)


## ----label="ANTsR pre-processing for group ICA2", echo=TRUE, eval=FALSE, results=FALSE----------------------
for( i in 1:9){
  fr <- paste0("sub-01_ses-105_task-rest_run-00",
               i, "_bold.nii.gz")
  file <- file.path(dataDir, "MyConnectome",
              "sub-01", "ses-105", "func", fr)
  fmrirest <- antsImageRead(file)
  pfmrirest <- preprocessfMRI(fmrirest,
            meanBoldFixedImageForMotionCorrection=avimg)
  pfile <- file.path(resDir, "MyConnectome",
        "sub-01", "ses-105", "func", paste0("prmc", fr))
  antsImageWrite(pfmrirest$cleanBoldImage, pfile)
  fileprestmask <- file.path(resDir, "MyConnectome",
                    "sub-01", "ses-105", "func",
              paste0("sub-01_ses-105_task-rest_run_00",
                     i, "_mask.nii.gz"))
  antsImageWrite(pfmrirest$maskImage, fileprestmask)
}


## ----label="sICA", echo=TRUE, eval=TRUE, results=FALSE, cache=TRUE------------------------------------------
mask <- as.array(readNIfTI(fprestmask, reorient=FALSE))
nii <- readNIfTI(frest, reorient = FALSE)
ds <- oro2fmri(nii)
ICAresult <- fmri.sICA(ds, mask, ncomp = 40,
          degree = 3, bws = 8, bwt = 4, unit = "FWHM")


## ----label="ICAfingerprint", echo=TRUE, eval=TRUE, results=FALSE, cache=TRUE--------------------------------
ICAresult <- ICAfingerprint(ICAresult)


## ----label="Figure_4_20", echo=TRUE, eval=TRUE, results=FALSE, fig.cap="Diagnostic plot for IC evaluation.", fig.width=10,fig.height=5.,warning=FALSE----
plot(ICAresult, 21)


## ----label="sICA group study", echo=TRUE, eval=TRUE, results=FALSE, cache=TRUE------------------------------
ICAresults <- list(NULL)
for(i in 1:9){
  fr <- paste0("prmcsub-01_ses-105_task-rest_run-00",
               i, "_bold.nii.gz")
  file <- file.path(resDir, "MyConnectome",
                    "sub-01", "ses-105",
                    "func", fr)
  nii <- readNIfTI(file, reorient = FALSE)
  ds <- oro2fmri(nii)
  ICAresults[[i]] <- fmri.sICA(ds, mask,
                               ncomp = 40, degree = 3,
                               bws = 8, bwt = 4,
                               unit = "FWHM",
                               alg.typ = "parallel",
                               alpha = 1.)
}
ICAresults[[10]] <- ICAresult
ICallruns <- fmri.sgroupICA(ICAresults, thresh = .9,
                           minsize = 6)


## ----label="Figure_4_21", echo=TRUE, eval=TRUE, results=FALSE, fig.cap="First independent component from combined IC analysis of 10 runs.", fig.width=9,fig.height=3.8,cache=FALSE,warning=FALSE----
plot(ICallruns, 1)


## ----"get mean BOLD image for registration",eval=TRUE,echo=TRUE,results=FALSE-------------------------------
pmcrest <- readNIfTI(fprest1, reorient = FALSE)
ameanpmcrest <- as.nifti(apply(abs(pmcrest), 1:3, mean),
                        pmcrest)
writeNIfTI(ameanpmcrest, mfprest0)


## ----"Register preprocessed resting data to MNI",eval=FALSE, echo=TRUE--------------------------------------
fmrirest <- antsImageRead(fprest1)
mni <- antsImageRead(fileMNI2)
prest2MNI <- antsRegistration(mni, mfprest, "Affine")
mnirest <- antsApplyTransforms(mni, fmrirest,
                prest2MNI$fwdtransforms, imagetype = 3)
antsImageWrite(mnirest, mnirestfile)


## ----"Access the Talairach atlas",eval=TRUE,echo=TRUE-------------------------------------------------------
atlas <- readNIfTI(fnTalairach)
rsdata <- readNIfTI(mnirestfile)
mask <- apply(abs(rsdata),1:3,max) > 0


## ----"read atlas metadata",echo=TRUE,eval=TRUE,warnings=FALSE-----------------------------------------------
Tmeta <- xmlTreeParse(file.path(fsl_atlas_dir(),
                                "Talairach.xml"))
Tmeta <- xmlToList(Tmeta)
regionNames <- as.character(Tmeta$data[1, ])
ind <- rep(0, length(regionNames))
getregion <- function(regionNames){
  ind <- rep(0, length(regionNames))
  for (i in 1:length(regionNames)){
    z <- strsplit(regionNames[i], "Brodmann area")
    if(length(z[[1]]) == 2){
      ind[i] <- as.numeric(z[[1]][2])
      if(length(grep("Left", z[[1]][1])) > 0)
        ind[i] <- ind[i] + 50
    }
  }
  ind
}
ind <- getregion(regionNames)
regions <- sort(unique(ind[ind > 0]))


## ----"aggregate atlas",eval=TRUE,echo=TRUE------------------------------------------------------------------
batlas <- atlas@.Data
batlas[, , ] <- 0
for(i in 1:97){
  jind <- (1:length(ind))[ind == i]
  batlas[atlas %in% jind] <- i
}
batlas[!mask] <- 0


## ----"create region names", eval=TRUE, echo=TRUE------------------------------------------------------------
ni <- numeric(length(regions))
for(i in 1:length(regions))
  ni[i] <- sum(batlas == regions[i])
regions <- regions[ni > 1]
nregions <- length(regions)
nL <- sum(regions > 50)
regionNames <- paste0(regions %% 50,
                      c(rep("R", nregions-nL),
                        rep("L", nL)))


## ----"create mean time series",eval=TRUE, echo=TRUE---------------------------------------------------------
mts <- matrix(0, 240, nregions)
dim(rsdata) <- c(prod(rsdata@dim_[2:4]),
                 rsdata@dim_[5])
for(i in 1:nregions)
   mts[, i] <- apply(rsdata[batlas == regions[i], ],
                     2, mean)
colnames(mts) <- regionNames


## ----"get correlation matrix", eval=TRUE, echo=TRUE---------------------------------------------------------
rcor <- cor(mts)


## ----"Figure_4_22",echo=TRUE,eval=TRUE,out.width='100%', fig.align="center",fig.cap="Functional  connectivity  between Brodmann areas: Heatmap of correlation matrix from resting state fMRI.",fig.width=10,fig.height=10----
heatmap(rcor, symm = TRUE)


## ----"get sparse partial correlation matrix",eval=TRUE,echo=TRUE--------------------------------------------
smat <- cov(mts)
pmat <- glasso(smat, .002)$wi
dsmath <- diag(1/sqrt(diag(pmat)))
pcor <- -dsmath %*% pmat %*% dsmath
diag(pcor) <- 1
dimnames(pcor) <- dimnames(rcor)


## ----"Construction of graphs",eval=TRUE,echo=TRUE-----------------------------------------------------------
diag(rcor) <- diag(pcor) <- 0
gcor <- graph_from_adjacency_matrix(rcor,
                                    "undirected",
                                    weighted = TRUE)
gpc <- graph_from_adjacency_matrix(pcor,
                                   "undirected",
                                   weighted = TRUE)


## ----"Figure_4_23",eval=TRUE,echo=TRUE,out.width='100%', fig.align="center",fig.cap="Functional connectivity between Brodmann areas: Connectivity graphs obtained from resting state fMRI using correlation (left) and partial correlation matrices (right).",fig.width=14,fig.height=8----
par(mfrow = c(1, 2),
    mar = c(3, 1, 3, .5), mgp = c(2, 1, 0))
colind <- min(255*(rcor+1)/2):max(255*(rcor+1)/2)
plot(permute(gcor,
             canonical_permutation(gcor)$labeling),
     layout = layout_in_circle, vertex.size = 12,
     edge.color = cm.colors(256)[colind])
title("functional connectivity (correlation)")
colind <- min(255*(pcor+1)/2):max(255*(pcor+1)/2)
plot(permute(gpc, canonical_permutation(gpc)$labeling),
     layout = layout_in_circle, vertex.size = 12,
     edge.color = cm.colors(256)[colind])
title("functional connectivity (partial correlation)")

