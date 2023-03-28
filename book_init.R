if(!exists("baseDir")) baseDir <- dirname(getwd())
codeDir <- file.path(baseDir,"MRIwithR") 
if(!dir.exists(codeDir)) stop(paste("incorrect baseDir:",baseDir, "directory MRIwithR not found")) 
dataDir <- file.path(baseDir,"MRIwithRdata")
if(!dir.exists(dataDir)) stop(paste("no data in directory baseDir:",baseDir, "directory MRIwithRdata not found")) 
resDir <- file.path(baseDir,"MRIwithRresults")
if(!dir.exists(resDir)){ 
   warning(paste("no intermediate results in directory baseDir:",baseDir, "directory MRIwithRresults is created")) 
   dir.create(resDir)
   }
