##This example simulates a staggered design with 3 cohort ages and an interim analysis

library(FAIRsimulator)
# 
# set.seed(324124)
# 
# StudyObj <- createStudy(latestTimeForNewBirthCohorts=0,studyStopTime = 6*30+1,
#                         nSubjects = c(1000,1000,1000),cohortStartTimes = c(0,0,0),
#                         samplingDesign = list(c(0,1,2,3,4,5,6)*30,c(0,1,2,3,4,5,6)*30,c(0,1,2,3,4,5,6)*30),
#                         randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
#                         minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
#                         treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
#                         effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
#                         Recruitmentfunction=function(...) {return(5000)},
#                         minSubjects = 10,newCohortLink = list(NULL,NULL,NULL),
#                         dropoutRates = rep(0.2/(6*30),3))

# #Set residual error to 0
# StudyObj$InitEvent <- function(StudyObj) {
#   StudyObj<-FAIRsimulator::InitEvent(StudyObj)
#   StudyObj$sig<-0 #Set residual error to 0
#   return(StudyObj)
# }
# 
# 
# StudyObj<-AdaptiveStudy(StudyObj)
# 
# ## Plot the design
# plotStudyCohorts(StudyObj)
# 
# ## plot the number of active subjects per treatment cycle and cohort
# plotActiveSubjects(StudyObj)
# 
# ## Plot the HAZ profiles versus age.
# plotHAZ(StudyObj)
# 
# ## Plot the HAZ data and the treatment effects
# plotHAZTreatmentEff(StudyObj)
# 
# # Plot the randomization probabilities at end of study
# plotProbs(StudyObj,strProb = "UnWeightedUpdateProbabilities")
# print("LMER Coefficients for each cohort:")
# StudyObj$CohortList[[1]]$UpdateCoefficients
# StudyObj$CohortList[[2]]$UpdateCoefficients
# StudyObj$CohortList[[3]]$UpdateCoefficients
