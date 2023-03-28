##
##   the call to fslr function probtrackx currently fails
##   need to use FSL command probtracks2 directly 
##   requires sample files xx<i>samples.nii.gz to be renamed merged_xx<i>samples.nii.gz
##
##   note that this is computationally expensive (2 days) and requires excessive memory (125 GByte)
##
if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Chapter_dwi","chapter_DWI_init.R"))


fast(file.path(rdwipd,"sub-01_ses-106_brain.nii.gz"),outfile=rdwipd)
img2 <- readNIfTI(file.path(rdwipd,"_pve_2.nii.gz"))
img1 <- readNIfTI(file.path(rdwipd,"_pve_1.nii.gz"))
writeNIfTI(img1*img2>0,file.path(rdwipd,"seeds"))

sdwipd <- file.path(ldwipd,"merged")
fnmask <- file.path(rdwipd, "sub-01_ses-106_brain.nii.gz")
seeds <- file.path(rdwipd,"seeds.nii.gz")

opt <- paste("--forcedir -V 2 -l --onewaycondition -c 0.2 -S 2000 
       --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=25.0 
       --sampvox=0.0 --opd --omatrix1 --dir=",ptdwipd)
cmd <- paste("probtrackx2 -s",sdwipd,"-m",fnmask,"-x",seeds,opt)
system(cmd)

# probrackx(samples = sdwipd, mask = fnmask, seed = fnmask, opt=opt)
