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
  DebugPrint(paste0("Check if cohort ",Cohort$Name," is about to have an interim analysis ",StudyObj$CurrentTime),3,StudyObj)
  TimeToPerformInterim<-FALSE
  tmp<-sum(unlist(lapply(Cohort$SubjectList,function(Subject,LastSample){
    return(as.numeric(Subject$Status==0 || (LastSample %in% Subject$SubjectSampleTime)))},max(Cohort$SamplingDesign))))/Cohort$MaxNumberOfSubjects
  if (length(tmp)==0) tmp<-0
  
  if (Cohort$RandomizationAgeRange[1]==6*30 && (StudyObj$CurrentTime==6*30-1 || StudyObj$CurrentTime==18*30-1 )) {
    DebugPrint(paste0("Time to do an interim analyses for cohort ",Cohort$Name," at time: ",StudyObj$CurrentTime," (",round(tmp*100,1)," % subjects completed)"),1,StudyObj)
    TimeToPerformInterim<-TRUE
  }
  return(TimeToPerformInterim)
}

FinalAnalysesEvent<-function(StudyObj) {
  DebugPrint(paste0("Check if time to do a final analyses at study time: ",StudyObj$CurrentTime),3,StudyObj)
  if (StudyObj$CurrentTime==StudyObj$StudyDesignSettings$StudyStopTime-1) {
    DebugPrint(paste0("Time to do a final analyses with all data at time: ",StudyObj$CurrentTime),1,StudyObj)
    df<-getAllSubjectData(StudyObj,prevTreatment = TRUE) #Get all data that is neede for final analyses
    ## Rename, select and order columns
    if(length(names(df)[grep(names(df),pattern = "PTRT")]) >0) {
      df <- df %>% 
        dplyr::rename(ID=StudyID,DATA=HAZ,AGE=Age,TRT=TreatmentIndex,TRTS=Treatment) %>% 
        select(ID,DATA,AGE,TRT,TRTS,one_of(StudyObj$StudyDesignSettings$Covariates),Level,CohortName,matches("PTRT"))
      for (cstr in names(df)[grep(names(df),pattern = "PTRT")]) { #Set no previous treatment to "0"
        df[,cstr]<-ifelse(is.na(df[,cstr]),0,df[,cstr])
      }
    } else {
      df <- df %>% 
        dplyr::rename(ID=StudyID,DATA=HAZ,AGE=Age,TRT=TreatmentIndex,TRTS=Treatment) %>% 
        select(ID,DATA,AGE,TRT,TRTS,one_of(StudyObj$StudyDesignSettings$Covariates),Level,CohortName)
    }
    
    df[df==-99]<-NA #Set -99 to missing
    df<-ImputeCovariates(df,StudyObj,method=StudyObj$StudyDesignSettings$ImpMethod) #Impute missing covariates
    
    #Sort to get correct TRT order
    df <- df[order(df["TRT"],df["ID"], df["AGE"]),] 
    
    ##### Make some covariate factors
    df$TRT<-as.factor(df$TRT) 
    df$SEXN<-as.factor(df$SEXN)
    df$SANITATN<-as.factor(df$SANITATN)
    df$AGE<-df$AGE/(12*30) #Rescale time to years
    
    ### If we have previous treatments
    ptrti <- grep("PTRT",names(df))
    if(length(ptrti)>0) {
      
      for(i in 1:length(ptrti)) {
        df[,ptrti[i]] <- as.factor(df[,ptrti[i]])
      }
      myPTRTs <- names(df)[ptrti]
    }
    
    cohortlevels<-unique(StudyObj$CohortList %listmap% "Level")
    browser()
    for (l in 1:length(cohortlevels)) {
      dflevel<-subset(df,Level==cohortlevels[l])
      
      #### Perform LME estimation based on some covariates and treatment effects for each level
      if(length(ptrti)>0) {
        
        lmeFormula <- paste0("DATA~1 + AGE + AGE:TRT + (AGE|ID) +",paste0(c(StudyObj$StudyDesignSettings$Covariates,myPTRTs),collapse = " + "))
      } else {
        lmeFormula <- paste0("DATA~1 + AGE + AGE:TRT + (AGE|ID) +",paste0(StudyObj$StudyDesignSettings$Covariates,collapse = " + "))
      }
      
      lmefit <- lmer(lmeFormula,data=dflevel,REML=FALSE) # IIV on baseline only
      
      ##### Calculate new probabilites based on another cohort LME results
      lmecoef<-summary(lmefit)$coefficients[,1] #Get coefficicents from LME
      lmese<-summary(lmefit)$coefficients[,2] #Get SE from LME
      lmecoef<-lmecoef[regexpr('AGE:TRT.*',names(lmecoef))==1]
      lmese<-lmese[regexpr('AGE:TRT.*',names(lmese))==1]
      print(lmefit)
      ### Get first cohort with level = l
      if (!is.null(StudyObj$CohortList)) {
        for (i in 1:length(StudyObj$CohortList)) {
          if (StudyObj$CohortList[[i]]$Level==l) {
            Cohort<-StudyObj$CohortList[[i]]
            break
          }
        }

        if ((length(lmecoef)+1)!=length(Cohort$RandomizationProbabilities)) {
          lmecoefnew<-rep(0,length(Cohort$RandomizationProbabilities)-1)
          lmesenew<-rep(0,length(Cohort$RandomizationProbabilities)-1)
          for (i in 1:(length(Cohort$RandomizationProbabilities)-1)) {
            iIndex<-which(names(lmecoef)==paste0("AGE:TRT",i+1))
            if (length(iIndex)!=0) {
              lmecoefnew[i]<-lmecoef[iIndex]
              lmesenew[i]<-lmese[iIndex]
            }
          }
          names(lmecoefnew)<-paste0("AGE:TRT",2:length(Cohort$RandomizationProbabilities))
          names(lmesenew)<-paste0("AGE:TRT",2:length(Cohort$RandomizationProbabilities))
          lmecoef<-lmecoefnew
          lmese<-lmesenew
        }
      }
      #Get probability of beeing best
      probs <- GetNewRandomizationProbabilities(trtcoeff=lmecoef,trtse=lmese,
                                                StudyObj$StudyDesignSettings$iNumPosteriorSamples) #Calculate randomization probs based on posterior distribution
      
      ## Apply any probability temperation function, e.g. sqrt
      probstemp <- StudyObj$StudyDesignSettings$probTemperation(probs)
    
      ## Save finala analyses output
      FinalAnalyses<-list()
      
      FinalAnalyses$Level = l
      FinalAnalyses$Time = StudyObj$CurrentTime
      FinalAnalyses$LMEFit<-lmefit
      FinalAnalyses$UnWeightedUpdateProbabilities<-probs
      FinalAnalyses$UnWeightedTemperatedUpdateProbabilities<-probstemp
      FinalAnalyses$LMECoeff<-lmecoef
      FinalAnalyses$LMESE<-lmese
      StudyObj$FinalAnalysesList[[length(StudyObj$FinalAnalysesList)+1]]<-FinalAnalyses
    }
  }
  return(StudyObj)
}


StudyObjIni <- createStudy(
  nCohorts                       = 2,
  recruitmentAges                = list(c(6,7)*30,c(12,13)*30),
  nSubjects                      = c(500,500),
  cohortStartTimes               = c(0*30,0*30),
  newCohortLink                  = list(NULL,NULL),
  Recruitmentfunction            = function(...) {return(5000)},
  samplingDesign                 = list(seq(0,12,by=2)*30,seq(0,6,by=1)*30),
  studyStopTime                  = 18*30+5,
  latestTimeForNewBirthCohorts   = 14*30+1,
  treatments                     = list(c("SoC-1","Cell 1","Cell 2"," Cell 3"," Cell 4"),c("SoC-1","Cell 1","Cell 2"," Cell 3"," Cell 4")),
  effSizes                       = list(c(0,0.0633,0.1037,0.1574,0.1687),c(0,0.0633,0.1037,0.1574,0.1687)),
  randomizationProbabilities     = list(rep(0.20,5),rep(0.20,5)),
  minAllocationProbabilities     = list(c(0.2,rep(0,4)),c(0.2,rep(0,4))),
  probTemperationFunction        = probTemperation,
  AddNewBirthCohortEventFunction = AddNewSixMonthCohortEvent,
  interimAnalyzesTimeFunction    = InterimAnalyzesTime,
  accumulatedData = TRUE
)

StudyObjIni$EventList[[length(StudyObjIni$EventList)+1]]<-AddNewTwelveMonthCohortEvent
StudyObjIni$EventList[[length(StudyObjIni$EventList)+1]]<-FinalAnalysesEvent

StudyObj <- AdaptiveStudy(StudyObjIni)



plotStudyCohorts(StudyObj,plotAnaTimes = T,shiftWithinLevel = 1)

plotProbs(StudyObj,cohortAgeNames=c("6-12 months","12-18 months"),strProb="UpdateProbabilities")

StudyObj$CohortList %listmap% "RandomizationProbabilities"
StudyObj$CohortList %listmap% "UpdateProbabilities"
StudyObj$CohortList %listmap% "UnWeightedUpdateProbabilities"

