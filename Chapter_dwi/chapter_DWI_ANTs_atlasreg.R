if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_dwi","chapter_DWI_init.R"))

load(file.path(rdwi,"fibertracks.rsc"))
pho <- readNIfTI(file.path(fsl_atlas_dir(),
                           "HarvardOxford",
                 "HarvardOxford-cort-prob-1mm.nii.gz"))
img <- as.nifti(apply(pho, 1:3, sum) > 0, pho)
fnhomask <- file.path(tmpdir, "HOmask")
writeNIfTI(img, fnhomask)

fnmask <- file.path(rdwipd,
                    "sub-01_ses-2015_brain_mask.nii.gz")
mask <- antsImageRead(fnmask)
homask <- antsImageRead(file.path(tmpdir,
                                  "HOmask.nii.gz"))
ho2dwi <- antsRegistration(fixed = mask,
                           moving = homask,
                           typeofTransform = "Affine")

fnatlas <- file.path(fsl_atlas_dir(), "HarvardOxford",
        "HarvardOxford-cort-maxprob-thr25-1mm.nii.gz")
atlas <- antsImageRead(fnatlas)
hoindwi <- antsApplyTransforms(mask, atlas,
  ho2dwi$fwdtransforms, interpolator = "genericLabel")
fnHOatlas <- file.path(rdwipd,
                       "sub-01_ses-2015_HOatlas.nii.gz")
antsImageWrite(hoindwi, fnHOatlas)

HOatlas <- readNIfTI(fnHOatlas)
zMT <- AdjacencyMatrix(trxcomb50, HOatlas)


set.seed(1)
gMT <- graph_from_adjacency_matrix(zMT,
                                   "undirected",
                                   weighted = TRUE)

HOmeta <- xmlTreeParse(file.path(fsl_atlas_dir(),
                       "HarvardOxford-Cortical.xml"))
HOmeta <- xmlToList(HOmeta)
regionNames <- as.character(HOmeta$data[1, ])
coords <- matrix(as.numeric(unlist(HOmeta$data[2, ])),
                 48, 4)[, 2:4]

save(regionNames,coords,gMT,zMT,file=file.path(rdwi,"Connectivity.rsc"))

