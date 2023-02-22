source("../CodeSecondEdition/chapter_DWI_init.R")
load(file.path(rdwi,"dataobj.rsc"))

setCores(4)
dtiobj.sm <- dti.smooth(dwobj,hmax=3)
dtiind.sm <- dtiIndices(dtiobj.sm)

for(i in 1:dwobj@ngrad){
  dwobj@si[, , , i] <-
    kernsm(dwobj@si[, , , i], h = 1)@yhat
}

fnsigma <- file.path(rdwipd, "sub-01_ses-106_sigma")
sigma <- as.array(readNIfTI(fnsigma))
setCores(6)
dtiobjql.gauss <- dtiTensor(dwobj,
                            method = "quasi-likelihood",
                            sigma = sigma)
dtiindql.gauss <- dtiIndices(dtiobjql.gauss)
save(dtiobjql.gauss, dtiindql.gauss, dtiobj.sm, dtiind.sm, 
     file=file.path(rdwi,"dtiobjsmooth.rsc"))
rm(dtiobjql.gauss, dtiindql.gauss, dtiobj.sm, dtiind.sm)
gc()
gc()

load(file.path(rdwi,"dataobj.rsc"))

sigmap <- awslsigmc(dwobj@si[ , , , 1],
                    steps = 16, mask=dwobj@mask)$sigma
dwobj <- dwi.smooth.ms(dwobj, kstar = 12,
                            sigma = sigmap)

dtiobjql.poas <- dtiTensor(dwobj,
                           method = "quasi-likelihood",
                           sigma = sigmap)
dtiindql.poas <- dtiIndices(dtiobjql.poas)

save(dtiindql.gauss, dtiindql.poas, file=file.path(rdwi,"dtiobjsmooth.rsc"))
