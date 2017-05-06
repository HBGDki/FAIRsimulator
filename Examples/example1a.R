library(FAIRsimulator)

set.seed(32423)

## Square root prob temperation
probTemperation <- function(probs) {
  probs <- sqrt(probs)/sum(sqrt(probs))
  return(probs)
}

## Symmetric prob temperation
probTemperation <- function(probs) {
  tmp<-function(x) {return(return(3.5*(x-0.5)^3+1/16*x+0.5))} 
  return(tmp(probs)/sum(tmp(probs)))
}

StudyObj <- createStudy(latestTimeForNewBirthCohorts=18*30,studyStopTime = 32*30,
                        nSubjects = c(320,320,320),
                        randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
                        #minAllocationProbabilities = list(c(0,rep(0,4)),c(0,rep(0,4)),c(0,rep(0,4))),
                        minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
                        treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
                        effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
                        Recruitmentfunction=function(...) {return(5000)},
                        probTemperationFunction = probTemperation,
                        Futilityfunction=futilityFunction,
                        checkFutility="before")



StudyObj<-AdaptiveStudy(StudyObj)


plotProbs(StudyObj)
plotProbs(StudyObj,strProb="UnWeightedRandomizationProbabilities")
