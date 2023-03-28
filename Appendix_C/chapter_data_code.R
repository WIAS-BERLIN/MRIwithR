if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","book_init.R"))


## ----"read code", echo=TRUE, eval=FALSE----------------------------------------------------------------
setwd(baseDir)
srcurl <- file.path("https://archive.wias-berlin.de",
                     "servlets","MCRFileNodeServlet",
                     "wias_derivate_00003583",
                     "MRIwithRcode.tgz")
 download.file(srcurl, file.path(baseDir,
                    "MRIwithRcode.tgz"), method="curl")
 untar(file.path(baseDir, "MRIwithRcode.tgz"))
 file.remove(file.path(baseDir, "MRIwithRcode.tgz"))
