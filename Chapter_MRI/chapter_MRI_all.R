## ----label="init mri", echo=FALSE, eval = TRUE, results='hide',warnings=FALSE,message=FALSE------------------------------------------------------
opts_knit$set(progress = TRUE, verbose = TRUE, warning=FALSE, message=FALSE,self.contained=FALSE,cache=FALSE)
options(width=50)
library(oro.nifti)


## ----"generate orthographic views", echo = c(1:8,10,13,16), eval=TRUE, results='hide',fig.keep='none',massage=FALSE------------------------------
if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","book_init.R"))

library(oro.nifti)
fn1 <- file.path(dataDir, "Kirby21", "T1",
                 "visit_1", "113",
                 "113-01-T1.nii.gz")
fn2 <- file.path(dataDir, "Kirby21", "T2",
                 "visit_1", "113",
                 "113-01-T2.nii.gz")
fn3 <- file.path(dataDir, "Kirby21", "FLAIR",
                 "visit_1", "113",
                 "113-01-FLAIR.nii.gz")
T1 <- readNIfTI(fn1)
T2 <- readNIfTI(fn2)
FLAIR <- readNIfTI(fn3)
png("figure/Figure_2_3a.png",
    width = 1000, height = 1250)
orthographic(T1,
        xyz = c(85.5, 128.5, 160), zlim = c(0, 1e6))
dev.off()
png("figure/Figure_2_3b.png",
    width = 1000, height = 1250)
orthographic(T2,
        xyz = c(90.5, 128.5, 160), zlim = c(0, 6e6))
dev.off()
png("figure/Figure_2_3c.png",
    width = 1000, height = 1250)
orthographic(FLAIR,
        xyz = c(164, 288.5, 360), zlim = c(0, 2e5))
dev.off()

