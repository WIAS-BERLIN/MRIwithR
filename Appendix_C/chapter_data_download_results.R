if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","book_init.R"))

## ----"read intermediate results", echo=TRUE, eval=FALSE------------------------------------------------
 setwd(baseDir)
 srcurl <- file.path("https://archive.wias-berlin.de",
                     "servlets","MCRFileNodeServlet",
                     "wias_derivate_00003583",
                     "MRIwithRresults.tgz")
 download.file(srcurl, file.path(baseDir,
               "MRIwithRresults.tgz"), method="curl")
 untar(file.path(baseDir, "MRIwithRresults.tgz"))
 file.remove(file.path(baseDir, "MRIwithRresults.tgz"))
