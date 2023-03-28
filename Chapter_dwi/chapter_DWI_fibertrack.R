if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_dwi","chapter_DWI_init.R"))

load(file=file.path(rdwi,"tensorIndices.rsc"))
load(file=file.path(rdwi,"dmtcomb.rsc"))
load(file=file.path(rdwi,"dtiobjsmooth.rsc"))

trxql <- tracking(dtiindql)
trxql <- selectFibers(trxql, minlength = 20)
trxql <- reduceFibers(trxql)

trxql.poas <- tracking(dtiindql.poas)
trxql.poas <- selectFibers(trxql.poas, minlength = 20)
trxql.poas <- reduceFibers(trxql.poas)

trxql.gauss <- tracking(dtiindql.gauss)
trxql.gauss <- selectFibers(trxql.gauss, minlength = 20)
trxql.gauss <- reduceFibers(trxql.gauss)


trxcomb5 <- tracking(dmtcomb, mincompartsize = .05)
trxcomb5 <- trxcomb50 <- reduceFibers(trxcomb5)
trxcomb5 <- selectFibers(trxcomb5, minlength = 20)

save(trxql, trxql.poas, trxcomb5, trxcomb50, file=file.path(rdwi,"fibertracks.rsc"))

rm(trxql, trxql.poas, trxcomb5,  trxcomb50, trxql.gauss)
