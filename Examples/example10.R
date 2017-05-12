### This example tries to mimic the standard design for the India meeting,
### In this case a 12-18 month cohort (with 3 batches of children)
### And a 6-12 month cohort (with 3 batches of children)

library(FAIRsimulator)

set.seed(324124)


## One adapted cohort
AddNewSixMonthCohortEvent<-function(StudyObj) {
  #Note - assuming birth cohort is first element in StudyDesignSettings
  
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      Cohort<-StudyObj$CohortList[[i]]
      if ((Cohort$Active==FALSE && Cohort$RandomizationAgeRange[1]==6*30 && Cohort$CycleNum==1 && Cohort$CohortStartTime+Cohort$CurrentTime==StudyObj$CurrentTime && StudyObj$CurrentTime<=StudyObj$StudyDesignSettings$LatestTimeForNewBirthCohorts) || #If this cohort just ended and its time to add a new cohort
         (Cohort$Active==TRUE && Cohort$RandomizationAgeRange[1]==6*30 && Cohort$CycleNum==1 && Cohort$CohortStartTime+6*30==StudyObj$CurrentTime && StudyObj$CurrentTime<=StudyObj$StudyDesignSettings$LatestTimeForNewBirthCohorts))
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

## One adapted cohort
AddNewTwelveMonthCohortEvent<-function(StudyObj) {
  #Note - assuming birth cohort is first element in StudyDesignSettings
  
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      Cohort<-StudyObj$CohortList[[i]]
      if (Cohort$Active==FALSE && Cohort$RandomizationAgeRange[1]==12*30 && Cohort$CycleNum==1 && Cohort$CohortStartTime+Cohort$CurrentTime==StudyObj$CurrentTime && StudyObj$CurrentTime<=StudyObj$StudyDesignSettings$LatestTimeForNewBirthCohorts) #if this cohort just ended and its time to add a new cohort
      {
        BirthCohortIndex<-which(unlist(lapply(StudyObj$StudyDesignSettings$CohortAgeRange,function(x) {x[1]==12*30}))==TRUE) #Get the birth cohort index (i.e. lower AgeRangeAtRandomization=0)
        NewChildCohort<-NewCohort(StudyObj,CohortNum=BirthCohortIndex)
        NewChildCohort$StartNum <- max(StudyObj$CohortList %listmap% "StartNum")+1
        NewChildCohort$Name<-paste0("C-",NewChildCohort$StartNum,"-",1," [",Cohort$RandomizationAgeRange[1]/30,"-",Cohort$RandomizationAgeRange[2]/30,"m @ rand]")
        DebugPrint(paste0("Create new 12 month cohort ",NewChildCohort$Name," based on probabilities in ",Cohort$Name," at time: ",StudyObj$CurrentTime),1,StudyObj)
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

InterimAnalyzesTime<-function(Cohort,StudyObj) {
  DebugPrint(paste0("Check if cohort ",Cohort$Name," is about to have an interim analysis",StudyObj$CurrentTime),3,StudyObj)
  TimeToPerformInterim<-FALSE
  tmp<-sum(unlist(lapply(Cohort$SubjectList,function(Subject,LastSample){
    return(as.numeric(Subject$Status==0 || (LastSample %in% Subject$SubjectSampleTime)))},max(Cohort$SamplingDesign))))/Cohort$MaxNumberOfSubjects
  if (length(tmp)==0) tmp<-0
  
  if (Cohort$RandomizationAgeRange[1]==6*30 && StudyObj$CurrentTime==6*30-1) {
    DebugPrint(paste0("Time to do an interim analyses for cohort ",Cohort$Name," at time: ",StudyObj$CurrentTime," (",round(tmp*100,1)," % subjects completed)"),1,StudyObj)
    TimeToPerformInterim<-TRUE
  }
  return(TimeToPerformInterim)
}



StudyObjIni <- createStudy(
  nCohorts                       = 2,
  recruitmentAges                = list(c(6,7)*30,c(12,13)*30),
  nSubjects                      = c(500,500),
  cohortStartTimes               = c(0*30,0*30),
  newCohortLink                  = list(NULL,NULL),
  Recruitmentfunction            = function(...) {return(5000)},
  samplingDesign                 = list(seq(0,12,by=2)*30,seq(0,6,by=1)*30),
  studyStopTime                  = 18*30+3,
  latestTimeForNewBirthCohorts   = 14*30+1,
  treatments                     = list(c("SoC-1","Cell 1","Cell 2"," Cell 3"," Cell 4"),c("SoC-1","Cell 1","Cell 2"," Cell 3"," Cell 4")),
  effSizes                       = list(c(0,0.0633,0.1037,0.1574,0.1687),c(0,0.0633,0.1037,0.1574,0.1687)),
  randomizationProbabilities     = list(rep(0.20,5),rep(0.20,5)),
  minAllocationProbabilities     = list(c(0.2,rep(0,4)),c(0.2,rep(0,4))),
  probTemperationFunction        = probTemperation,
  AddNewBirthCohortEventFunction = AddNewSixMonthCohortEvent,
  interimAnalyzesTimeFunction    = InterimAnalyzesTime
)

StudyObjIni$EventList[[length(StudyObjIni$EventList)+1]]<-AddNewTwelveMonthCohortEvent


StudyObj <- AdaptiveStudy(StudyObjIni)



plotStudyCohorts(StudyObj,plotAnaTimes = T,shiftByLevel = 0.2)

StudyObj$CohortList %listmap% "RandomizationProbabilities"
StudyObj$CohortList %listmap% "UpdateProbabilities"
StudyObj$CohortList %listmap% "UnWeightedUpdateProbabilities"

