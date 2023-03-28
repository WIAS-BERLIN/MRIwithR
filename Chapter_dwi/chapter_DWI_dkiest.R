if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_dwi","chapter_DWI_init.R"))

load(file.path(rdwi,"dataobj.rsc"))
## ----"Compute dkitensor",eval = FALSE,echo=TRUE, results='hide'---------------------------------------------
dkiobj <- dkiTensor(dwobj,
                    method = "CLLS-QP")

## ----"Compute dkiIndices",eval = TRUE, echo=TRUE, cache=TRUE, results='hide',warning=FALSE,message=FALSE----
dkiind <- dkiIndices(dkiobj)

save(dkiobj, dkiind, file=file.path(rdwi,"dkiObjects.rsc"))

