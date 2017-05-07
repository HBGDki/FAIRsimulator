#This example have 3 cohort's but a 12-24 cohort instead of a 12-18 month cohort and an interim-analysis
#The design also includes a weight function for the probabilities that down weight both low and high probabilities to be more robust versus futility treatments
#but with the drawback at converging at a slower rate
#The design also have a recruitment rate which gives ~ 20 days for recruitment

library(FAIRsimulator)
 
set.seed(32156)
 
StudyObj <- createStudy(latestTimeForNewBirthCohorts=10*30,studyStopTime = 12*30+20,
                         nSubjects = c(500,500,500),cohortStartTimes = c(0,0,0),
                         samplingDesign = list(c(0,1,2,3,4,5,6)*30,c(0,1,2,3,4,5,6)*30,seq(0,12)*30),
                         randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
                         minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
                         treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
                         effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
                         Recruitmentfunction=function(...) {return(25)},
                        newCohortLink = list(2,NULL,NULL),
                         dropoutRates = rep(0.2/(6*30),3),
                         Futilityfunction = function(probs,...){return(probs)},
                         probTemperationFunction =  function(probs) {
                          tmp<-function(x) {return(return(3.5*(x-0.5)^3+1/16*x+0.5))} 
                           return(tmp(probs)/sum(tmp(probs)))
                         },recruitmentAges = list(c(0,1)*30,c(6,7)*30,c(12,13)*30))

 
StudyObj<-AdaptiveStudy(StudyObj)
 
## Plot the design
plotStudyCohorts(StudyObj)
 
## plot the number of active subjects per treatment cycle and cohort
plotActiveSubjects(StudyObj)
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
