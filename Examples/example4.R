## This example simulates a 4 cohort design with 4 different age cohorts and 5 treatments, allowing for new cohorts at birth
## A recruitment rate of 32 subjects per day is assumed, futility function is not present

library(FAIRsimulator)
set.seed(324124)
 
StudyObj <- createStudy(latestTimeForNewBirthCohorts=12,studyStopTime = 13*30+1,
                        nCohorts = 4,
                        nSubjects = c(320,320,320,320),cohortStartTimes = c(0,0,0,0),
                        samplingDesign = list(c(0,1,2,3,4,5,6)*30,c(0,1,2,3,4,5,6)*30,c(0,1,2,3,4,5,6)*30,c(0,1,2,3,4,5,6)*30),
                        randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5),rep(0.20,5)),
                        minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
                        treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12"),c("SoC-4","TRT-13","TRT-14","TRT-15","TRT-16")),
                        effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
                        Recruitmentfunction=function(...) {return(32)},
                        minSubjects = 10,newCohortLink = list(2,3,4,NULL),
                        recruitmentAges = list(c(0,1)*30,c(6,7)*30,c(12,13)*30,c(18,19)*30),
                        dropoutRates = rep(0.2/(6*30),4),
                        Futilityfunction = function(probs,...){return(probs)})

StudyObj<-AdaptiveStudy(StudyObj)
 
## Plot the design
plotStudyCohorts(StudyObj)
 
## plot the number of active subjects per treatment cycle and cohort
plotActiveSubjects(StudyObj)
 
## Plot the HAZ profiles versus age.
plotHAZ(StudyObj)
