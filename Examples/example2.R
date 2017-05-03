## This example simulates the same design multiple times in parallel.

library(FAIRsimulator)
library(foreach)
library(doParallel)
library(purrr)

set.seed(32423)


probTemperation <- function(probs) {
  probs <- sqrt(probs)/sum(sqrt(probs))
  return(probs)
}

StudyObjIni <- createStudy(latestTimeForNewBirthCohorts=18*30,studyStopTime = 32*30,
                        nSubjects = c(320,320,320),
                        randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
                        #minAllocationProbabilities = list(c(0,rep(0,4)),c(0,rep(0,4)),c(0,rep(0,4))),
                        minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
                        treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
                        effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
                        Recruitmentfunction=function(...) {return(5000)},
                        minSubjects = 10,
                        probTemperationFunction = probTemperation)



registerDoParallel(cores=6)

myRes <- foreach(i=1:12) %dopar% AdaptiveStudy(StudyObjIni)

tmp <- map(myRes,getProbData)
bind_rows(tmp)

## Extract the probabilities used in the plotProbs plot above.



tmp <- getProbData(StudyObj)
