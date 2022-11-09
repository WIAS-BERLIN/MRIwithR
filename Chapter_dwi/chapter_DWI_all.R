## ----label="init dwi", echo=FALSE, eval = TRUE, results='hide',warning=FALSE,message=FALSE----------------------
opts_knit$set(progress = TRUE, verbose = TRUE, warning=FALSE, message=FALSE)
library(rgl)
knit_hooks$set(rgl = hook_rgl)
options(width=50)


## ----label="init2", echo=TRUE, eval = TRUE, results='hide',warning=FALSE,message=FALSE--------------------------
#opts_chunk$set(warning=FALSE,message=FALSE)
installed <- installed.packages(fields="")
library(dti)
library(adimpro)
library(aws)
library(rgl)
library(igraph)
haveFSLr <- require(fslr)
library(XML)
library(oro.nifti)
library(jsonlite)
haveANTsR <- require(ANTsR)
library(igraph)
library(XML)


## ----label="init3", echo=TRUE, eval=TRUE, message=FALSE---------------------------------------------------------
dataDir <- "../data"
resDir <- "../results"
dwid <- file.path(dataDir, "MyConnectome", "sub-01",
                  "ses-106", "dwi")
rdwipd <- file.path(resDir, "MyConnectome", "sub-01",
                   "ses-106", "dwi-proc")
ldwipd <- file.path(resDir, "MyConnectome", "sub-01",
                   "ses-106", "logdir")
rdwi <- file.path(resDir,"DWI")
dwijson <-  file.path(dataDir, "MyConnectome",
                      "sourcedata", "dicom_headers",
                      "sub-01", "ses-106", "dwi")
tmpdir <- tempdir()
rimage.options(zquantiles=c(0.001,0.98),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")


## ----"determine code to be executed, no need to reproduce results already computed in code snippets", eval=TRUE, echo=FALSE----
preprocessed <- file.exists(file.path(rdwipd,"sub-01_ses-106_dwi_proc.nii.gz"))
# have results from chapter_DWI_preprocess.R
if(!haveFSLr&&!processed) stop("You need to install fslr or directory results")
havedataobj <- file.exists(file.path(rdwi,"dataobj.rsc"))
# data will be saved after first run of the document
havesigma <- file.exists(file.path(rdwipd, "sub-01_ses-106_sigma.nii.gz"))
# have results from chapter_DWI_varest.R
havetensor <- file.exists(file.path(rdwi,"tensorObjects.rsc"))
# have results from chapter_DWI_tensorest.R
havedkitensor <- file.exists(file.path(rdwi,"dkiObjects.rsc"))
# have results from chapter_DWI_dkiest.R
haveQball <- file.exists(file.path(rdwi,"qBallObjects.rsc"))
# have results from chapter_DWI_qBallest.R
havemixtensor <- file.exists(file.path(rdwi,"dmtcomb.rsc"))
# have results from chapter_DWI_mixtens.R
havesmooth <- file.exists(file.path(rdwi,"dtiobjsmooth.rsc"))
# have results from chapter_DWI_smooth.R
havexfibres <- file.exists(file.path(ldwipd,"dyads2.nii.gz"))
# have results from chapter_DWI_xfibres.R
haveprobtrackx <- file.exists(file.path(ldwipd,"pfiber","fdt_parhs.nii.gz"))
# have results from chapter_DWI_probtrackx.R
havefibers <- file.exists(file.path(rdwi,"fibertracks.rsc"))
# have results from chapter_DWI_fibertrack.R
haveconnect <- file.exists(file.path(rdwi,"Connectivity.rsc"))
if(!haveANTsR&&!haveconnect) stop("You need to install ANTsR or directory results")
# have results from chapter_DWI_ANTs_atlasreg.R


## ----label="Figure_5_1", echo=-1, eval = TRUE, results=FALSE, cache=FALSE, message=FALSE, fig.cap="Effect of phase encoding direction on a non-diffusion weighted image. Left: A-P encoding direction, right: P-A encoding direction. The distortions are rather large in, but not restricted to, the anterior brain regions (upper image part). The both images can be combined to a single distortions-corrected one.", fig.width=6, fig.height=3.3, out.width='66%', fig.align='center'----
par(mfrow = c(1, 2),
    mar = c(1, 1, 2, 0.1), mgp = c(2, 1, 0))
run1FileName <- "sub-01_ses-106_run-001_dwi.nii.gz"
run2FileName <- "sub-01_ses-106_run-002_dwi.nii.gz"
run1File <- file.path(dwid, run1FileName)
run2File <- file.path(dwid, run2FileName)
ds1 <- readNIfTI(run1File)
ds2 <- readNIfTI(run2File)
rimage(ds1[ , , 41, 2], main = "A-P phase encoding")
rimage(ds2[ , , 41, 2], main = "P-A phase encoding")


## ----label="initialize dirs", echo=TRUE, eval = TRUE, results=FALSE, message=FALSE------------------------------
rdwipd <- file.path(resDir, "MyConnectome", "sub-01",
                   "ses-106", "dwi-proc")
if (!dir.exists(rdwipd)) dir.create(rdwipd)






















## ----label="Figure_5_2", echo=-1, eval = TRUE, results=FALSE, message=FALSE, fig.cap="Result of the correction for the susceptibility-induced artifacts and the derived brain map on a non-diffusion weighted image.", fig.width=8, fig.height=4.3, out.width='66%', fig.align='center'----
par(mfrow = c(1, 2),
    mar = c(1, 1, 2, 0.1), mgp = c(2, 1, 0))
topio <- file.path(rdwipd, "sub-01_ses-106_topupiout")
fnmask <- file.path(rdwipd, "sub-01_ses-106_brain")
dscorrected <- readNIfTI(topio)
mask <- readNIfTI(paste0(fnmask, "_mask"))
rimage(dscorrected[ , , 41, 1], main = "Corrected image")
rimage(mask[ , , 41], main = "Estimated brain mask")














## ----label = "showSummaryDWOBJ", echo = TRUE, eval = TRUE, message = FALSE--------------------------------------
summary(dwobj)


## ----label = "showDimensionOfObject", echo = TRUE---------------------------------------------------------------
dim(dwobj@si)


## ----label = "Figure_5_3", echo = -2, eval = TRUE, results = 'markup', fig.cap = "An axial slice of the the first 6 volumes of the example diffusion weighted data. The first and second volume are non-diffusion weighted image. Images are scaled for maximum contrast. This demonstrates that the diffusion weighted images have a smaller signal-to-noise ratio.", fig.width = 6, fig.height = 5.7, out.width='75%',fig.align='center'----
rimage.options(swapy=TRUE)
par(mfrow = c(2, 3),
    mar = c(.5, .5, 2, 0.1), mgp = c(1.5, .5, 0))
for (vol in 1:6) {
  rimage(dwobj@si[, , 41, vol],
         main = paste("Volume", vol, ",",
                      "b-value", dwobj@bvalue[vol]))
}


## ----label = "Figure_5_4", echo = -(1:2), eval = TRUE, results = 'markup', fig.cap = "The same axial slice of three volumes acquired with b-values 0, 1500 and 3000 $s/mm^2$. The same intensity range is used for all images.", fig.width = 6, fig.height = 2.92, fig.align='center',out.width='75%'----
dti.options(swapy=TRUE)
par(mfrow = c(1, 3),
    mar=c(.5, .5, 2, 0.1), mgp = c(1.5, .5, 0),
    xaxt = "n", yaxt = "n")
for (i in 16:18) {
  img <- plot(dwobj, slice = 41, gradient = i)
  title(paste("b-value", dwobj@bvalue[i]))
}


## ----label = "writeSlice", echo = TRUE, eval = FALSE------------------------------------------------------------
## library(adimpro)
## write.image(img, file="myslice.png")


## ----label = "Figure_5_5", echo = TRUE, eval = TRUE, results = 'hide', rgl = TRUE, fig.cap = "3D visualization of diffusion weighted data in the object \\code{dwobj} by the function \\code{show3d}. The data values are coded as distance from the center in the direction of the diffusion gradients. Color codes gradient directions.", fig.width = 6, fig.height = 3, fig.align='center'----
show3d(dwobj,
       xind = 41:44, yind = 55:56, zind = 41,
       what = "data",
       zoom = .5, scale= c(0.8, 0.8, 0.8),
       windowRect = c(0, 0, 1000, 400),
       userMatrix = rotationMatrix(0, 1, 0, 1))


## ----label = "Figure_5_6", echo = TRUE, eval = TRUE, results = 'markup', rgl = TRUE, fig.cap = "3D visualization of the regular polyeders \\code{icosa0}, \\code{icosa1}, \\code{icosa2}, \\code{icosa3} and \\code{icosa4} from the dataset \\code{polyeders} in package \\pkg{dti}.", fig.width = 10, fig.height = 1.88, fig.align='center'----
data(polyeders)
icosa <- list(icosa0, icosa1, icosa2, icosa3, icosa4)
rgl.open()
par3d(zoom = 0.15, windowRect = c(0, 0, 1300, 250))
for(i in 1:5){
  ni <- icosa[[i]]$ni * 3
  ind <- icosa[[i]]$indices
  rgl.triangles(icosa[[i]]$vertices[1, ind] + 3*(i-1),
                icosa[[i]]$vertices[2, ind],
                icosa[[i]]$vertices[3, ind],
                color = rep("red", ni))
}


## ----label = "Figure_5_7", echo = TRUE, eval = TRUE, results = 'hide',  message = FALSE, fig.cap = "Threshold definition using the function \\code{sdpar}.", fig.width = 8.7, fig.height = 4, fig.align='center'----
dwobj <- sdpar(dwobj, level = 1.75)


## ----label = "define mask", echo=TRUE, eval = TRUE, results='markup'--------------------------------------------
lmask <- getmask(dwobj)


## ----label = "Figure_5_8", echo = -1, eval = TRUE,  fig.cap = "Mask definitions, central axial slice. Left: Mean non-diffusion weighted image, center: mask defined using threshold and \\code{sdpar}, right: mask obtained from \\code{topup}.", fig.width = 6., fig.height = 2.9, fig.align='center', out.width='75%'----
par(mfrow=c(1,3),mar=c(1,1,2,.1),mgp=c(2,1,0))
rimage(lmask$s0[,,41],col=grey(0:255/255))
title("Mean image")
rimage(lmask$mask[,,41],col=grey(0:255/255))
title("Mask from thresholding ")
rimage(dwobj@mask[,,41],col=grey(0:255/255))
title("Mask from fslr::topup")


## ----"Figure_5_9",eval = TRUE,echo=FALSE,results='hide',fig.cap="Effect of low SNR for Rician distribution. Densities for $\\sigma=1$. Vertical dashed and dotted lines are located at theta and the corresponding expected value, the length of the horizontal lines indicate the gap between both values.", fig.width=8, fig.height=4----
library(dti)
setCores(6)
Lhalf <- function(x,L=1){
   gamma(L+1/2)/gamma(1.5)/gamma(L)*dti:::hg1f1(-0.5,L, -x)
}
drician <- function(y,theta=0,sigma=1){
  y/sigma^2*exp(-(y^2+theta^2)/2/sigma^2)*besselI(y*theta/sigma^2,0)
}
Erician <- function(theta,sigma=1){
   sigma*sqrt(pi/2)*Lhalf(-theta^2/2/sigma^2)
}
par(mar=c(3,3,3,1),mgp=c(2,1,0))
x <- seq(0,5.5,.01)
plot(x,drician(x,0),type="l",ylab="Density value")
lines(x,drician(x,1),col=2)
lines(x,drician(x,2),col=3)
lines(x,drician(x,3),col=4)
lines(x,drician(x,4),col=5)
lines(c(0,0),c(0,.3),lty=2,col=1,lwd=2)
lines(rep(Erician(0),2),c(0,.3),lty=3,col=1)
lines(c(0,Erician(0)),c(.05,.05),col=1,lend=2,lwd=2)
lines(c(1,1),c(0,.3),lty=2,col=2,lwd=2)
lines(rep(Erician(1),2),c(0,.3),lty=3,col=2)
lines(c(1,Erician(1)),c(.1,.1),col=2,lend=2,lwd=2)
lines(c(2,2),c(0,.3),lty=2,col=3,lwd=2)
lines(rep(Erician(2),2),c(0,.3),lty=3,col=3)
lines(c(2,Erician(2)),c(.15,.15),col=3,lend=2,lwd=2)
lines(c(3,3),c(0,.3),lty=2,col=4,lwd=2)
lines(rep(Erician(3),2),c(0,.3),lty=3,col=4)
lines(c(3,Erician(3)),c(.2,.2),col=4,lend=2,lwd=2)
lines(c(4,4),c(0,.3),lty=2,col=5)
lines(rep(Erician(4),2),c(0,.3),lty=3,col=5)
lines(c(4,Erician(4)),c(.25,.25),col=5,lend=2,lwd=2)
legend(4.5,.6,paste0("theta=",0:4),col=1:5,lty=rep(1,5))
title("Rician distributions")












## ----label="Figure_5_10",echo=-1,eval = TRUE,fig.cap="Local noise estimation using \\code{aflsigmc}.",fig.width=8,fig.height=2.2,fig.align='center'----
par(mfrow = c(1, 4),
    mar = c(1, 1, 2, .1), mgp = c(2, 1, 0),
    xaxt = "n", yaxt = "n")
rimage(afsigma1b0$sigma[, , 41])
title("Effective sigma run 1")
rimage(afsigma1b0$Leff[, , 41])
title("Effective L run 1")
rimage(afsigma2b0$sigma[, , 41])
title("Effective sigma run 2")
rimage(afsigma2b0$Leff[, , 41])
title("Effective L run 2")
rm(afsigma1b0,afsigma2b0)








## ----"Figure_5_11", eval = TRUE, echo= -1,fig.cap="Local estimates of sigma and theta for non-diffusion  weighted (A-P and P-A) images. Right: Combined estimate of sigma (top) and densities of the estimated sigma values (bottom).",fig.width=8,fig.height=5.9,out.width='75%',fig.align='center'----
par(mfrow = c(2, 3), mar = c(1.5, 1.5, 2, .1),
    mgp = c(1, .5, 0))
rimage(sb1$sigma[, , 41], main = "sigma bv=0 A-P")
rimage(sb2$sigma[, , 41], main = "sigma bv=0 P-A")
rimage(sigma[, , 41], main = "sigma bv=0 combined")
rimage(sb1$theta[, , 41], main = "theta A-P")
rimage(sb2$theta[, , 41], main = "theta P-A")
d0 <- density(sb1$sigma, to = .6)
d1 <- density(sb2$sigma)
d2 <- density(sigma)
ylim <- range(d0$y, d1$y, d2$y)
plot(d0, ylim = ylim, ylab = "", xlab = "",
     main = "density(sigma)")
lines(d1, col = 2)
lines(d2, col = 3)
legend(0.3, ylim[2], c("A-P", "P-A", "combined"),
       col = 1:3, lty = rep(1,3))
rm(sb1, sb2)


## ----label="Figure_5_12",eval = TRUE, echo=TRUE,fig.cap="ADC for selected voxel.",fig.width=8,fig.height=3,rgl=TRUE,results='hide'----
show3d(dwobj,
       what = "adc",
       xind = 31:33, yind = 53, zind = 41,
       scale = 1, zoom = .4)





## ----"DTI model 1",eval = TRUE----------------------------------------------------------------------------------
dim(dtiobj@D)


## ----"DTI model 2",eval = TRUE----------------------------------------------------------------------------------
signif(dtiobj@D[, 56, 56, 41],2)


## ----"DTI model 3",eval = TRUE----------------------------------------------------------------------------------
ind <- c(1, 2, 3, 2, 4, 5, 3, 5, 6)
DTensor <- matrix(dtiobj@D[ind, 56, 56, 41], 3, 3)
adcDT <- diag(t(dtiobj@gradient) %*%
                DTensor %*% dtiobj@gradient)
adcDT[3:5]




## ----"nl DTI model 2",eval = TRUE-------------------------------------------------------------------------------
signif(matrix(dtiobjnl@D[, 56, 56, 41][ind],
              c(3, 3)), 3)




## ----"qlDTI model res",eval = TRUE, results='markup'------------------------------------------------------------
signif(matrix(dtiobjql@D[, 56, 56, 41][ind],
              c(3, 3)), 3)


## ----"qlDTI model summary",eval = TRUE, results='markup'--------------------------------------------------------
summary(dtiobjql)


## ----"meanDiffusivity",echo=TRUE, eval=FALSE--------------------------------------------------------------------
## md <- dwiMD(dwobj)




## ----label="Figure_5_13",echo=-1,eval = TRUE, fig.cap="Left: FA map of an axial slice of the example dMRI data. Right:  Corresponding map of the geodesic anisotropy. The latter has been up-scaled in the value range for more image contrast. Thus, large GA values are saturated.", fig.width=5.3, fig.height=3.85,out.width='66%',fig.align='center'----
 par(mfrow = c(1, 2),
     mar = c(1, 1, 2, .1), mgp = c(2, 1, 0))
 rimage(dtiindnl@fa[, , 41], main = "FA map")
 rimage(dtiindnl@ga[, , 41], zlim = c(0, 1.5),
        main = "GA map")


## ----"fig:dti3dprep",echo=-c(6,11,17), eval=TRUE,results='hide'-------------------------------------------------
show3d(dtiobjql,
       xind = 30:39, yind = 28:34, zind = 41,
       scale = .7, zoom = .6, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0, 1, 0, 0))
rgl.snapshot("figure/Figure_5_14a.png")
show3d(dtiobjql,
       xind = 37, yind = 32, zind = 41,
       subdivide = 4, zoom = .7, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700))
rgl.snapshot("figure/Figure_5_14b.png")
show3d(dtiindql,
       xind = 30:39, yind = 28:34, zind = 41,
       zoom = .6, lwd = 4, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0, 1, 0, 0))
rgl.snapshot("figure/Figure_5_14c.png")


## ----"Figure_5_15",echo=-1,eval = TRUE, fig.cap="Upper row: FA maps of an axial slice of the example dMRI data obtained using the three available methods for tensor estimation. Lower row:  Corresponding map of the geodesic anisotropy. The latter has been up-scaled in the value range for more image contrast. Thus, large GA values are saturated.",fig.width=6.6,fig.height=6.4, out.width='75%',fig.align='center'----
par(mfrow = c(2, 3),
    mar = c(.1, .1, 2, .1), mgp = c(2, 1, 0),
    xaxt = "n", yaxt = "n")
plot(dtiind, slice = 41)
title("FA linear model")
plot(dtiindnl, slice = 41)
title("FA non-linear model")
plot(dtiindql, slice = 41)
title("FA quasi-likelihood")
plot(dtiind, what = "ga", slice = 41)
title("GA linear model")
plot(dtiindnl, what = "ga", slice = 41)
title("GA non-linear model")
plot(dtiindql, what = "ga", slice = 41)
title("GA quasi-likelihood")
rm(dtiobj,dtiobjnl,dtiobjql,dtiind,dtiindnl)


## ----"generate muenster color scale",eval = FALSE,echo=FALSE,cache=TRUE,results='hide'--------------------------
## showFAColorScale("figure/Figure_5_16.png")


## ----"Figure_5_17",eval = TRUE,echo=-1,fig.cap="Color coded FA using various available color schemes.",fig.width=8,fig.height=2.92----
par(mfrow = c(1, 4), mar = c(.5, .5, 2, .1),
    mgp = c(2, 1, 0))
for (i in c(1, 3, 4, 6)){
  plot(dtiindql, slice = 41, method = i,
       xaxt = "n", yaxt = "n", mar = par("mar"))
  title(paste0("method=", i))
}





## ----"DTensor from dki",echo=TRUE,eval = TRUE-------------------------------------------------------------------
dim(dkiobj@D)
dim(dkiobj@W)


## ----"show3d on dkiobj",echo=c(1:6, 8:13),eval=TRUE,results='hide',cache=TRUE-----------------------------------
show3d(dkiobj,
       xind = 30:39, yind = 28:34, zind = 41,
       what = "DT",
       zoom = .6, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0, 1, 0, 0))
rgl.snapshot("figure/Figure_5_18a.png")
show3d(dkiobj,
       xind = 30:39, yind = 28:34, zind = 41,
       what = "KT",
       zoom = .6, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0, 1, 0, 0))
rgl.snapshot("figure/Figure_5_18b.png")




## ----"Figure_5_19",eval = TRUE,echo=-1,fig.cap="Indices defined for the DKI model. From left to right: mean diffusivity, fractional (tensor) anisotropy, apparent kurtosis, and fractional kurtosis.", fig.width=8,fig.height=2.9----
par(mfrow = c(1, 4),
    mar = c(.5, .5, 2, .1), mgp = c(2, 1, 0),
    xaxt = "n", yaxt = "n")
plot(dkiind, slice = 41, what = "md", main="MD")
plot(dkiind, slice = 41, what = "fa", main="FA")
plot(dkiind, slice = 41, what = "mk2",
     zlim = c(0, .025), main="Kapp")
plot(dkiind, slice = 41, what = "fak", main="FAK")
rm(dkiind,dkiobj)


## ----"arttensor",echo= -c(19,25),eval = TRUE,results='hide'-----------------------------------------------------
data(polyeders)
D <- diag(c(1, .3, .3))[c(1, 2, 3, 5, 6, 9)]
grad <- cbind(c(0, 0, 0), icosa2$vertices)
btb <- dti:::create.designmatrix.dti(grad)
bv <- c(0, rep(1, 162))
si <- exp(-bv*D%*%btb)
si <- as.nifti(array(si, c(1, 1, 1, 163)))
fn <- file.path(tmpdir,"tensor")
writeNIfTI(si, fn)
dstens <- readDWIdata(grad,
                      file.path(tmpdir,"tensor.nii.gz"),
                      "NIFTI",
                      bvalue = bv)
odf0 <- dwiQball(dstens,
                 what = "ODF",
                 mask = array(TRUE, c(1, 1, 1)),
                 order = 6)
show3d(odf0, bgcolor = "white")
rgl.snapshot("figure/Figure_5_20a.png")
odfw <- dwiQball(dstens,
                 what = "wODF",
                 mask = array(TRUE, c(1, 1, 1)),
                 order = 6)
show3d(odfw, bgcolor = "white")
rgl.snapshot("figure/Figure_5_20b.png")





## ----"Figure_5_21",eval = TRUE, echo=TRUE, cache=TRUE, results='hide',fig.cap="Qball-reconstruction order 4 in a region from a central slice.", fig.width=10, fig.height=7, rgl=TRUE----
show3d(qballw4,
       xind = 30:39, yind = 28:34, zind = 41,
       scale = .7, zoom = .6, bgcolor = "white",
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0, 1, 0, 0))



## ----"Qball figures", eval=TRUE, echo=FALSE, results='hide',cache=TRUE------------------------------------------
show3d(qballw4,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
rgl.snapshot("figure/Figure_5_22w4e2.png")
show3d(qballw4a,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
rgl.snapshot("figure/Figure_5_22w4e3.png")
show3d(qballw4c,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
rgl.snapshot("figure/Figure_5_22w4e1.png")
show3d(qballw6a,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
rgl.snapshot("figure/Figure_5_22w6e3.png")
show3d(qballw6b,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
rgl.snapshot("figure/Figure_5_22w6e2.png")
show3d(qballw6c,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
rgl.snapshot("figure/Figure_5_22w6e1.png")
show3d(qballw8a,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
rgl.snapshot("figure/Figure_5_22w8e3.png")
show3d(qballw8b,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
rgl.snapshot("figure/Figure_5_22w8e2.png")
show3d(qballw8c,xind=34,yind=37,zind=41, subdivide=4, windowRect=c(0,0,512,400),zoom=.6, bgcolor="white")
rgl.snapshot("figure/Figure_5_22w8e1.png")


## ----"Figure_5_23",eval = TRUE,echo=TRUE,fig.cap="ODF's overlayed on a tensor FA map.", fig.width=6, fig.height=4.2,results='hide'----
overlay <- function(fn1, fn2, thresh = 100){
  img1 <- read.image(fn1)
  img2 <- read.image(fn2)
  img1 <- shrink.image(img1,
                       xt = img2$dim[1],
                       yt = img2$dim[2],
                       ratio = FALSE)
  over <- function(x1, x2, thresh = 100){
    if(length(dim(x2)) == 2){
       x2[x2 < thresh] <- x1[x2 < thresh]
    } else {
       ind <- apply(x2 < thresh, c(1, 2), "all")
       dim(x2) <- c(prod(dim(x2)[1:2]), dim(x2)[3])
       for(i in 1:dim(x1)[3])
         x2[ind, i] <- x1[, , i][ind]
       dim(x2) <- dim(x1)
    }
    x2
  }
  combine(img1, img2, over, thresh = thresh)
}
img <- plot(dtiindql[30:39, 28:34, 41], show = FALSE)
write.image(img, file="figure/coloredFA.png")
show3d(qballw4,
       xind = 30:39, yind = 28:34, zind = 41,
       scale = 0.7, zoom = 0.59,
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0*pi, 1, 0, 0))
rgl.snapshot("figure/wODFoverlay.png")
par(mar = c(0, 0, 0, 0))
show.image(overlay("figure/coloredFA.png",
                   "figure/wODFoverlay.png", 100),
           xaxt = "n", yaxt = "n")

## ----"rm Qball",eval=TRUE,echo=FALSE----------------------------------------------------------------------------
rm(qballw4,qballw4a,qballw4c,qballw6a,qballw6b,qballw6c,qballw8a,qballw8b,qballw8c)
gc()


## ----"Using fsl/xfibres", eval=FALSE, echo=TRUE-----------------------------------------------------------------
## gradientFile <- file.path(rdwipd,
##                           "sub-01_ses-106_dwi_proc.eddy_rotated_bvecs")
## bvalueFile <- file.path(rdwipd,
##                         "sub-01_ses-106_dwi_proc.bval")
## dataFile <- file.path(rdwipd,
##                       "sub-01_ses-106_dwi_proc")
## fnmask <- file.path(rdwipd,"sub-01_ses-106_brain_mask")
## xfibres(dataFile, gradientFile, bvalueFile, fnmask, 1,
##         verbose=FALSE, opts="--model=2")


## ----label="initadimpro3", echo=FALSE, eval=TRUE, message=FALSE-------------------------------------------------
rimage.options(zquantiles=c(0.001,0.98),xlab="x",ylab="z",bty="n",xaxt="n",yaxt="n")








## ----"Figure_5_25", eval = TRUE, echo=TRUE, results='hide',fig.cap="Weighted ODF for the tensor mixture model described in the text for the voxel selection of Figure~\\ref{fig:Figure_5_23} estimated using  the dwiMixtensor function of package \\pkg{dti}.", fig.width=8,fig.height=5,rgl=TRUE----
show3d(dmtcomb,
       xind = 30:39, yind = 28:34, zind = 41,
       subdivide = 3, scale = .4, zoom = .6,
       windowRect = c(0, 0, 1000, 700),
       userMatrix = rotationMatrix(0*pi, 1, 0, 0))




## ----"Figure_5_26", eval = TRUE, echo=-1, results='hide',fig.cap="Tensor mixture model of maximal order five: Isotropic compartment size, fractional anisotropy (FA), order of mixture and effective order (EO) of mixture (from left to right) for the central axial slice.",fig.width=10.,fig.height=3.6----
par(mar = c(.5, .5, 2, .1), mgp = c(2, 1, 0))
plot(dmtcomb, what = c("w0", "fa",
                       "order",
                       "eorder"),
     slice = 41, xaxt = "n", yaxt = "n")












## ----"Load msPOAS results",eval = FALSE,echo=FALSE, cache=TRUE--------------------------------------------------
## load("tmp/DWI/dwobjpoas.rsc")




## ----"Figure_5_27",eval = TRUE,echo=-1,fig.cap="Color coded FA maps for tensor estimates using the original data (left), data smoothed by POAS (center) and Gaussian smoothing (right).",fig.width=8.,fig.height=3.83, out.width='80%',fig.align='center'----
par(mfrow = c(1, 3),
    mar = c(.5, .5, 2, .1), mgp = c(2, 1, 0))
plot(dtiindql, slice = 41, xaxt = "n", yaxt = "n")
plot(dtiindql.poas, slice = 41, xaxt = "n", yaxt = "n")
plot(dtiindql.gauss, slice = 41, xaxt = "n", yaxt = "n")
rm(sigma,dtiindql.gauss)


## ----"tensor again",echo=TRUE,eval=FALSE------------------------------------------------------------------------
## dtiobj <- dtiTensor(dwobj)


## ----"tensor indices again",echo=TRUE,eval=FALSE----------------------------------------------------------------
## dtiind <- dtiIndices(dtiobj)
## dtiind@fa[!mask] <- 0






## ----"rm dmtcomb", eval=TRUE, echo=FALSE------------------------------------------------------------------------
rm(dmtcomb)


## ----"Tractography illustrations",eval = FALSE,echo=c(1:3,5:6,8:9)----------------------------------------------
## wRec <- c(0, 0, 1000, 1000)
## show3d(trxql,
##        zoom = .5, windowRect = wRec, bgcolor = "white")
## rgl.snapshot("figure/Figure_5_29a.png")
## show3d(trxql.poas,
##        zoom = .5, windowRect = wRec, bgcolor = "white")
## rgl.snapshot("figure/Figure_5_29b.png")
## show3d(trxcomb5,
##        zoom = .5, windowRect = wRec, bgcolor = "white")
## rgl.snapshot("figure/Figure_5_29c.png")


## ----"Probablilistic fiber tracking",eval=FALSE,echo=TRUE-------------------------------------------------------
## probtrackx(samples = "logdir",
##            mask = fnmask, seed = fnmask)













## ----"rm tracks", eval=TRUE, echo=FALSE-------------------------------------------------------------------------
rm(trxql, trxql.poas, trxcomb5)


## ----"Figure_5_31", eval=TRUE, echo=TRUE,out.width='100%', fig.align="center",fig.cap="Anatomic connectivity for the HarvardOxford cortical atlas: Adjacency matrices obtained from fiber tracking results for the tensor mixture model.",fig.width=10,fig.height=10----
heatmap(zMT, symm = TRUE, col = grey(0:255/255))
title("Tensor Mixture Model")




## ----"Figure_5_32", eval=TRUE, echo=-1,out.width='100%', fig.align="center",fig.cap="Anatomic connectivity for the HarvardOxford cortical atlas: Visualization of connectivity graphs obtained from fiber tracking results for the tensor mixture model.",fig.width=11,fig.height=6----
par(mfrow = c(1, 2),
    mar = c(3, 3, 3, 1), mgp = c(2, 1, 0))
plot(permute(gMT,canonical_permutation(gMT)$labeling),
     layout = layout_on_sphere,
     edge.color = grey(33:100/100))
plot(permute(gMT,canonical_permutation(gMT)$labeling),
     layout = layout_in_circle,
     edge.color = grey(33:100/100))
title("Connectivity graph from MT")

