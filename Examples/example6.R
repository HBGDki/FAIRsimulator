## This example simulates a design with 3 cohort ages, a recruitment rate of 20 choldren/day and an interim analysis after 5 months

library(FAIRsimulator)
set.seed(8585)

InterimAnalyzesTime<-function(Cohort,StudyObj) {
   TimeToPerformInterim <-FALSE
  
  if (Cohort$CurrentTime %in% c(5*30)) {
    TimeToPerformInterim<-TRUE
  }
   
  return(TimeToPerformInterim)
}


StudyObjIni <- createStudy(
  cohortStartTimes = c(0,0,0),
  newCohortLink = list(2, 3, NULL),
  recruitmentAges = list(c(0,1)*30,c(6,7)*30,c(12,13)*30),
  nSubjects = c(300,300,300),
  Recruitmentfunction=function(...) {return(20)},
  samplingDesign = list(0:6*30, seq(0,6,by=2)*30, seq(0,6,by=2)*30),
  studyStopTime = 25*30+3,
  latestTimeForNewBirthCohorts=0*30,
  treatments =list(
    c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),
    c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),
    c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
  effSizes = list(c(0,0.05,0.1,0.15,0.25),
                  c(0,0.05,0.1,0.15,0.25),
                  c(0,0.05,0.1,0.15,0.25)),
  randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
  strCovariates = c("BIRTHWT","MAGE", "MHTCM", "SEXN", "SANITATN"),
  minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4)))
)

# Run the study
StudyObj <- AdaptiveStudy(StudyObjIni)

## Plot the design and interim analysis time points
plotStudyCohorts(StudyObj,plotAnaTimes = T)

## Plot the active subjects
plotActiveSubjects(StudyObj)

# Check the randomization probabilities based on the interim analysis in cohort 4
StudyObj$CohortList[[4]]$PreviousRandomizationProbabilities[[1]]$RandomizationProbabilities

# Check the randomization probabilities based on the final analysis in cohort 4
StudyObj$CohortList[[4]]$RandomizationProbabilities
