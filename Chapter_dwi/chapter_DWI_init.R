if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","book_init.R"))
## load packages needed in DWI chapter
tmpdir <- tempdir()

dwid <- file.path(dataDir, "MyConnectome", "sub-01",
                  "ses-2015", "dwi")
rdwipd <- file.path(resDir, "MyConnectome")
if (!dir.exists(rdwipd)) dir.create(rdwipd)
rdwipd <- file.path(resDir, "MyConnectome", "sub-01")
if (!dir.exists(rdwipd)) dir.create(rdwipd)
rdwipd <- file.path(resDir, "MyConnectome", "sub-01",
                    "ses-2015")
if (!dir.exists(rdwipd)) dir.create(rdwipd)
rdwipd <- file.path(resDir, "MyConnectome", "sub-01",
                    "ses-2015", "dwi-proc")
if (!dir.exists(rdwipd)) dir.create(rdwipd)
rdwi <- file.path(resDir,"DWI")
if (!dir.exists(rdwi)) dir.create(rdwi)
ldwipd <- file.path(resDir, "MyConnectome", "sub-01",
                   "ses-2015", "logdir")
if (!dir.exists(ldwipd)) dir.create(ldwipd)
ptdwipd <- file.path(resDir, "MyConnectome", "sub-01",
                    "ses-2015", "pfiber")
if (!dir.exists(ldwipd)) dir.create(ldwipd)

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
library(XML)
setCores(16)

rimage.options(zquantiles=c(0.001,0.98),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
