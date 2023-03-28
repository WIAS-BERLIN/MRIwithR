if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_dwi","chapter_DWI_init.R"))

load(file.path(rdwi,"dataobj.rsc"))
## ----"dwiQball1",eval = TRUE, echo=TRUE, cache=TRUE, results='hide'-----------------------------------------
qballw4 <- dwiQball(dwobj,
                    what = "wODF",
                    order = 4, lambda = 1e-2)
qballw4a <- dwiQball(dwobj, what="wODF", order=4, lambda=1e-3)
qballw4c <- dwiQball(dwobj, what="wODF", order=4, lambda=1e-1)
qballw6a <- dwiQball(dwobj, what="wODF", order=6, lambda=1e-3)
qballw6b <- dwiQball(dwobj, what="wODF", order=6, lambda=1e-2)
qballw6c <- dwiQball(dwobj, what="wODF", order=6, lambda=1e-1)
qballw8a <- dwiQball(dwobj, what="wODF", order=8, lambda=1e-3)
qballw8b <- dwiQball(dwobj, what="wODF", order=8, lambda=1e-2)
qballw8c <- dwiQball(dwobj, what="wODF", order=8, lambda=1e-1)

save(qballw4, qballw4a, qballw4c, qballw6a, qballw6b, qballw6c, 
     qballw8a, qballw8b, qballw8c, file=file.path(rdwi,"qBallObjects.rsc"))

