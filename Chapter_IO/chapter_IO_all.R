## ----label="init-io", echo=FALSE, eval = TRUE, results='hide',warnings=FALSE,message=FALSE----------
opts_knit$set(progress = TRUE, verbose = TRUE, warning=FALSE, message=FALSE)
library(oro.dicom)
library(oro.nifti)
options(width=50)


## ----label="initadimpro2", echo=FALSE, eval=TRUE, message=FALSE----------
library(adimpro)
rimage.options(zquantiles=c(0.001,0.98),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")


## ----"DICOM data",echo=TRUE,eval=TRUE----------
if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","book_init.R"))

files <- dir(file.path(dataDir,
                       "example-dicom", "dicoms"),
             all.files = TRUE, full.names = TRUE,
             no.. = TRUE)


## ----label="Load a single DICOM file", echo=TRUE, results='markup'----------
ds <- readDICOMFile(files[250])


## ----label="Show dimension of DICOM image data", echo=TRUE, results='markup'----------
dim(ds$img)


## ----label="Figure_3_1", eval=TRUE, echo=2:3, results='markup', out.width="66%", fig.align="center", fig.cap="Slice number 250  of the example T1 dataset.", fig.width=8, fig.label="Show image of DICOM data"----
par(mar = c(3, 3, 3, .1), mgp = c(2, 1, 0))
rimage(1:384, 1:274, ds$img)


## ----label="Show header element of DICOM data", echo=TRUE, results='markup'----------
ds$hdr[c(36, 37, 54), ]


## ----label="Extract header elements of DICOM data", echo=TRUE, results='markup'----------
extractHeader(ds$hdr, "InstanceNumber")


## ----label="Read DICOM data from directory general", eval = TRUE, echo=TRUE, results='markup'----------
dicomdir <- file.path(dataDir, "example-dicom",
                      "dicoms")
dsseries <- readDICOM(dicomdir)


## ----label="Reproduce figure", eval=FALSE, echo=TRUE, results='markup'----------
## rimage(dsseries$img[[250]])


## ----label="Find number of files", echo=TRUE, results='markup'----------
length(dsseries$img)


## ----label="Extract metadata", echo=TRUE, results='markup'----------
head(extractHeader(dsseries$hdr, "InstanceNumber"))


## ----"Convert to NIfTI", echo=TRUE, eval=TRUE, warning=FALSE, message=FALSE----------
niftiobj <- dicom2nifti(dsseries)
analyzeobj <- dicom2analyze(dsseries)


## ----"Write NIfTI/ANALYZE", echo=TRUE, eval=TRUE----------
writeNIfTI(niftiobj, file.path(dataDir,
                               "example-dicom", "T1"))
writeANALYZE(analyzeobj, file.path(dataDir,
                               "example-dicom", "T1"),
             gzipped = FALSE)


## ----label="Read ANALYZE and NIfTI files general", eval=TRUE, echo=TRUE, results='markup'----------
ds1 <- readANALYZE(file.path(dataDir,
                             "example-dicom", "T1"))
ds2 <- readNIfTI(file.path(dataDir,
                           "example-dicom", "T1"))


## ----label="Read ANALYZE and NIfTI files", echo=TRUE, results='markup'----------
slotNames(ds1)
ds1@pixdim


## ----label="Figure_3_2", echo=TRUE, results='markup', out.width='75%', fig.align="center", fig.cap="Orthographic view of the T1 image."----------
orthographic(ds1, zlim = c(0, 100))


## ----label="simple mask",eval=TRUE,echo=TRUE----------
mask <- ds2 > 20


## ----label="create NIfTI object", echo=TRUE, results='markup'----------
mask <- as.nifti(mask, ds2)


## ----label="write NIfTI object", eval=FALSE, echo=TRUE, results='markup'----------
## write(mask, '<filename>')


## ----label="AFNI processing", echo=TRUE, eval=FALSE, results='markup'----------
## ds <- readAFNI(<filename>)
## writeAFNI(ds, <filename>)


## ----label="BIDS", echo=TRUE, eval=TRUE, results='markup'----------
library(jsonlite)
fname <- file.path(dataDir, "MyConnectome",
                   "sub-01", "ses-2015", "anat",
                   "sub-01_ses-2015_T1w.json")
ttt <- read_json(fname)
head(names(ttt),10)
ttt$RepetitionTime
ttt$EchoTime
ttt$FlipAngle

