if(!exists("baseDir")) baseDir <- dirname(dirname(getwd()))
source(file.path(baseDir,"MRIwithR","Appendix_A","chapter_nonpara_init.R"))

sigma2 <- mean(aws.gaussian(
   T1[121:220, 121:200, 121:200], 10)@sigma2)

T1aws6 <- aws(T1, hmax = 1.15, sigma2 = sigma2)@theta
T1aws12 <- aws(T1, hmax = 1.7, sigma2 = sigma2)@theta
T1aws18 <- aws(T1, hmax = 2.7, sigma2 = sigma2)@theta
T1aws24 <- aws(T1, hmax = 4.2, sigma2 = sigma2)@theta
T1aws30 <- aws(T1, hmax = 6.6, sigma2 = sigma2)@theta

T1paws6 <- paws(T1, hmax = 1.15, sigma2 = sigma2)@theta
T1paws12 <- paws(T1, hmax = 1.7, sigma2 = sigma2)@theta
T1paws18 <- paws(T1, hmax = 2.7, sigma2 = sigma2)@theta
T1paws24 <- paws(T1, hmax = 4.2, sigma2 = sigma2)@theta
T1paws30 <- paws(T1, hmax = 6.6, sigma2 = sigma2)@theta

save(sigma2, T1aws6, T1aws12, T1aws18, T1aws24, T1aws30,
        T1paws6, T1paws12, T1paws18, T1paws24, T1paws30, 
        file=file.path(npresDir,"awsestimate.rsc"))
