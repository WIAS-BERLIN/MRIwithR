if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Appendix_A","chapter_nonpara_init.R"))

z1 <- awstestprop(c(200, 200, 200), hmax = 6)
z2 <- awstestprop(c(200, 200, 200), hmax = 6, ladjust = 1.5)

save(z1, z2, file=file.path(npresDir,"propagation.rsc"))
