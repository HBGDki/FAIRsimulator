## This example simulates one cohort with a relatively slow recruitment and a large number of subjects.

## It doesn't work right now, but the aim is to have an interim analysis after 6 months and then recruit all subsequent subjects acording to these probabilities.

library(FAIRsimulator)

set.seed(324124)

StudyObj <- createStudy(
  nCohorts = 2,
  latestTimeForNewBirthCohorts=7*30,
  studyStopTime = 24*30,
  nSubjects = c(320,320),
  recruitmentAges = list(c(0,1),c(0,1)),
  samplingDesign = list(c(0,1,2,3,4,5,6)*30,c(0,1,2,3,4,5,6)*30),
  randomizationProbabilities = list(rep(0.20,5),rep(0.20,5)),
  minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4))),
  treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8")),
  effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
  Recruitmentfunction=function(...) {return(5000)})

StudyObj <- createStudy(
  nCohorts = 1,
  latestTimeForNewBirthCohorts = 0,
  studyStopTime = 12*30+1,
  nSubjects = 1000,
  cohortStartTimes = 0,
  samplingDesign = list(c(0,1,2,3,4,5,6)*30),
  randomizationProbabilities = list(rep(0.20,5)),
  minAllocationProbabilities = list(c(0.2,rep(0,4))),
  treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4")),
  effSizes = list(c(0,0.05,0.1,0.15,0.25)),
  Recruitmentfunction=function(...) {return(10)},
  newCohortLink = NULL,
  dropoutRates = rep(0.2/(6*30),3))

# #Set residual error to 0
# StudyObj$InitEvent <- function(StudyObj) {
#   StudyObj<-FAIRsimulator::InitEvent(StudyObj)
#   StudyObj$sig<-0 #Set residual error to 0
#   return(StudyObj)
# }


StudyObj<-AdaptiveStudy(StudyObj)

## Plot the design
plotStudyCohorts(StudyObj)

## plot the number of active subjects per treatment cycle and cohort
plotActiveSubjects(StudyObj)

## Plot the HAZ profiles versus age.
plotHAZ(StudyObj)

## Plot the HAZ data and the treatment effects
plotHAZTreatmentEff(StudyObj)

# Plot the randomization probabilities at end of study
#plotProbs(StudyObj,strProb = "UnWeightedUpdateProbabilities")

tmp <- getProbData(StudyObj,strProb = "UnWeightedUpdateProbabilities")
ggplot(tmp,aes(x=TreatmentName,y=Prob,fill=TreatmentName)) +
  geom_bar(stat="identity") +
  geom_text(aes(x=TreatmentName,y=Prob+0.02,label=paste("P =",round(Prob,2)))) +
  ylim(0,1) +
  ylab("Probability") +
  xlab("Treatment") +
  labs(fill=NULL,linetype=NULL) + 
  theme(legend.position="top") 
