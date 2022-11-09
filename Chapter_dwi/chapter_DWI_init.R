## load packages needed in DWI chapter
dataDir <- "../data"
resDir <- "../results"
tmpdir <- tempdir()
dwid <- file.path(dataDir, "MyConnectome", "sub-01",
                  "ses-106", "dwi")
dwijson <-  file.path(dataDir, "MyConnectome",
                      "sourcedata", "dicom_headers",
                      "sub-01", "ses-106", "dwi")
rdwipd <- file.path(resDir, "MyConnectome", "sub-01",
                   "ses-106", "dwi-proc")
if (!dir.exists(rdwipd)) dir.create(rdwipd)
rdwi <- file.path(resDir,"DWI")
if (!dir.exists(rdwipd)) dir.create(rdwi)
ldwipd <- file.path(resDir, "MyConnectome", "sub-01",
                   "ses-106", "logdir")

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
setCores(8)

rimage.options(zquantiles=c(0.001,0.98),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
