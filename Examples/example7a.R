## Simulates 3 cohorts with the same recruitment age but with start times 0 and 6 months. The second and third cohort have its randomization probabilities udpated based on the first/second cohort.

## It doesn't work right now.

library(FAIRsimulator)

set.seed(324124)

StudyObjIni <- createStudy(
  nCohorts                     = 1,
  recruitmentAges              = list(c(6,7)*30),
  nSubjects                    = c(500),
  cohortStartTimes             = c(0*30),
  newCohortLink                = list(NULL),
  Recruitmentfunction          = function(...) {return(5000)},
  samplingDesign               = list(seq(0,12,by=2)*30),
  studyStopTime                = 36*30+3,
  latestTimeForNewBirthCohorts = 24*30+1,
  treatments                   = list(c("SoC-1","Cell 1","Cell 2"," Cell 3"," Cell 4")),
  effSizes                     = list(c(0,0.0633,0.1037,0.1574,0.1687)),
  randomizationProbabilities   = list(rep(0.20,5)),
  minAllocationProbabilities   = list(c(0.2,rep(0,4))),
  AddNewBirthCohortEventFunction = AddNewSixMonthCohortEvent
)

StudyObj <- AdaptiveStudy(StudyObjIni)

plotStudyCohorts(StudyObj,plotAnaTimes = T)

StudyObj$CohortList %listmap% "RandomizationProbabilities"
StudyObj$CohortList %listmap% "UpdateProbabilities"
StudyObj$CohortList %listmap% "UnWeightedUpdateProbabilities"


## One adapted cohort
AddNewSixMonthCohortEvent<-function(StudyObj) {
  #Note - assuming birth cohort is first element in StudyDesignSettings
  
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      Cohort<-StudyObj$CohortList[[i]]
      if (Cohort$Active==FALSE && Cohort$RandomizationAgeRange[1]==6*30 && Cohort$CycleNum==1 && Cohort$CohortStartTime+Cohort$CurrentTime==StudyObj$CurrentTime && StudyObj$CurrentTime<=StudyObj$StudyDesignSettings$LatestTimeForNewBirthCohorts) #If this cohort just ended and its time to add a new cohort
      {
        BirthCohortIndex<-which(unlist(lapply(StudyObj$StudyDesignSettings$CohortAgeRange,function(x) {x[1]==6*30}))==TRUE) #Get the birth cohort index (i.e. lower AgeRangeAtRandomization=0)
        NewChildCohort<-NewCohort(StudyObj,CohortNum=BirthCohortIndex)
        NewChildCohort$StartNum <- max(StudyObj$CohortList %listmap% "StartNum")+1
        NewChildCohort$Name<-paste0("C-",NewChildCohort$StartNum,"-",1," [",Cohort$RandomizationAgeRange[1]/30,"-",Cohort$RandomizationAgeRange[2]/30,"m @ rand]")
        DebugPrint(paste0("Create new 6 month cohort ",NewChildCohort$Name," based on probabilities in ",Cohort$Name," at time: ",StudyObj$CurrentTime),1,StudyObj)
        NewChildCohort$RandomizationProbabilities<-Cohort$UpdateProbabilities #Use the latest probabilities from previous birth cohort
        NewChildCohort$UnWeightedRandomizationProbabilities<-Cohort$UnWeightedUpdateProbabilities
        StudyObj$CohortList[[(length(StudyObj$CohortList)+1)]]<-NewChildCohort
        break
      }
    }
  }
  return(StudyObj)
}


probTemperation <- function(probs) {
  probs <- sqrt(probs)/sum(sqrt(probs))
  return(probs)
}

StudyObjIni <- createStudy(
  recruitmentAges = list(c(0,1)*30,c(6,7)*30,c(12,13)*30),
  nSubjects = c(300,300,300),
  Recruitmentfunction=function(...) {return(5000)},
  samplingDesign = list(0:6*30, c(0,3,6)*30, c(0,3,6)*30),
  studyStopTime = 32*30,
  latestTimeForNewBirthCohorts=18*30,
  treatments =list(
    c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),
    c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),
    c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
  effSizes = list(c(0,0.05,0.1,0.15,0.25),
                  c(0,0.05,0.1,0.15,0.25),
                  c(0,0.05,0.1,0.15,0.25)),
  randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
  strCovariates = c("BIRTHWT","MAGE", "MHTCM", "SEXN", "SANITATN"),
  minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
  accumulatedData = TRUE,updateProbabilitiesFunction = UpdateProbabilities
)

StudyObj <- AdaptiveStudy(StudyObjIni)

plotStudyCohorts(StudyObj,plotAnaTimes = TRUE)

tmp1 <- getItemsFromSubjects(StudyObj$CohortList[[7]],prevTreatmen=TRUE)
tmp2 <- getItemsFromSubjects(StudyObj$CohortList[[6]],prevTreatmen=TRUE)

tmp3 <- full_join(tmp1,tmp2)
names(tmp2)
