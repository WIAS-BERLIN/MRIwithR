source("../CodeSecondEdition/chapter_DWI_init.R")
load(file.path(rdwi,"dataobj.rsc"))
setCores(64) 

## ----"Tensor mixtures", eval=FALSE, echo=TRUE, results='hide',cache=TRUE------------------------------------
dmtobj5 <- dwiMixtensor(dwobj, maxcomp = 5, model = "MTiso")
dmtobj4 <- dwiMixtensor(dwobj, maxcomp = 4, model = "MTiso")
dmtobj3 <- dwiMixtensor(dwobj, maxcomp = 3, model = "MTiso")
dmtobj2 <- dwiMixtensor(dwobj, maxcomp = 2, model = "MTiso")
dmtobj1 <- dwiMixtensor(dwobj, maxcomp = 1, model = "MTiso")
dmtcomb <- dwiMtCombine(dmtobj5, dmtobj4)
dmtcomb <- dwiMtCombine(dmtcomb, dmtobj3)
dmtcomb <- dwiMtCombine(dmtcomb, dmtobj2)
dmtcomb <- dwiMtCombine(dmtcomb, dmtobj1)

mtindices <- extract(dmtcomb,
                    what = c("w0", "fa", "eorder", "order"))

save(dmtcomb, file=file.path(rdwi,"dmtcomb.rsc"))
rm(dmtobj1,dmtobj2,dmtobj3,dmtobj4,dmtobj5,dmtcomb,
     mtindices)


