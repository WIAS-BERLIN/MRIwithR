if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_dwi","chapter_DWI_init.R"))

load(file.path(rdwi,"dataobj.rsc"))

dtiobj.sm <- dti.smooth(dwobj,hmax=3)
dtiind.sm <- dtiIndices(dtiobj.sm)

for(i in 1:dwobj@ngrad){
  dwobj@si[, , , i] <-
    kernsm(dwobj@si[, , , i], h = 1)@yhat
}

fnsigma <- file.path(rdwipd, "sub-01_ses-2015_sigma")
sigma <- as.array(readNIfTI(fnsigma))
dtiobjql.gauss <- dtiTensor(dwobj,
                            method = "quasi-likelihood",
                            sigma = sigma)
dtiindql.gauss <- dtiIndices(dtiobjql.gauss)

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
