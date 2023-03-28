if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_dwi","chapter_DWI_init.R"))

load(file.path(rdwi,"dataobj.rsc"))
fnsigma <- file.path(rdwipd, "sub-01_ses-2015_sigma")
sigma <- as.array(readNIfTI(fnsigma))
sigma <- sigma[dwobj@xind,dwobj@yind,]
## ----"DTI model",eval = FALSE,echo=TRUE,results='hide'------------------------------------------------------
dtiobj <- dtiTensor(dwobj,
                    method = "linear")


## ----"DTI model 1",eval = TRUE------------------------------------------------------------------------------
dim(dtiobj@D)


## ----"DTI model 2",eval = TRUE------------------------------------------------------------------------------
signif(dtiobj@D[, 56, 56, 41],2)


## ----"DTI model 3",eval = TRUE------------------------------------------------------------------------------
ind <- c(1, 2, 3, 2, 4, 5, 3, 5, 6)
DTensor <- matrix(dtiobj@D[ind, 56, 56, 41], 3, 3)
adcDT <- diag(t(dtiobj@gradient) %*%
                DTensor %*% dtiobj@gradient)
adcDT[3:5]


## ----"nlDTI model",eval = FALSE,echo=TRUE,results='hide'----------------------------------------------------
dtiobjnl <- dtiTensor(dwobj,
                      method = "nonlinear")


## ----"nl DTI model 2",eval = TRUE---------------------------------------------------------------------------
signif(matrix(dtiobjnl@D[, 56, 56, 41][ind],
              c(3, 3)), 3)


## ----"qlDTI model",eval=FALSE,echo=2:9,results='hide'-------------------------------------------------------
## ## sigma contains negative values introduced in topup
thresh <- quantile(sigma[dwobj@mask], .1)
sigma[sigma < thresh] <- thresh
dtiobjql <- dtiTensor(dwobj,
                      method = "quasi-likelihood",
                      sigma = sigma)
# this runs  >90 min


## ----"qlDTI model res",eval = TRUE, results='markup'--------------------------------------------------------
signif(matrix(dtiobjql@D[, 56, 56, 41][ind],
              c(3, 3)), 3)


## ----"qlDTI model summary",eval = TRUE, results='markup'----------------------------------------------------
summary(dtiobjql)

if(!dir.exists(rdwi)) dir.create(rdwi)
save(dtiobj,dtiobjnl,dtiobjql, file=file.path(rdwi,"tensorObjects.rsc"))




## ----"dtiIndices", echo=TRUE, eval = FALSE, results='hide'--------------------------------------------------
dtiind <- dtiIndices(dtiobj)
dtiindnl <- dtiIndices(dtiobjnl)
dtiindql <- dtiIndices(dtiobjql)

save(dtiind,dtiindnl,dtiindql, file=file.path(rdwi,"tensorIndices.rsc"))
rm(dtiobj,dtiobjnl,dtiobjql,dtiind,dtiindnl,dtiindql)
