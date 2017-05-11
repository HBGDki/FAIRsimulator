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


## Only Cohort data

##

myCohorts <- cohorts(StudyObj)

GetCohortData(myCohorts[[7]],StudyObj,accumulatedData=TRUE) %>% View

GetCohortData<-function(Cohort,StudyObj, accumulatedData = FALSE) {
  DebugPrint(paste0("Assembly cohort ",Cohort$Name," data available at time ",StudyObj$CurrentTime),2,StudyObj)

  if(!accumulatedData) {
  ## Get the data from the current cohort
  myDf <- getItemsFromSubjects(Cohort,
                               scalarItems=c("StudyID","TreatmentIndex","Treatment"),
                               longitudinalItems = c("Data","SampleAge"))
  } else {
    
    # Loop over the other cohorts and get the data from those of the same level
    myDf <- getItemsFromSubjects(Cohort,
                                 scalarItems=c("StudyID","TreatmentIndex","Treatment"),
                                 longitudinalItems = c("Data","SampleAge"),
                                 prevTreatment = TRUE)

    myCohorts <- cohorts(StudyObj)
    
    for(i in 1:length(myCohorts)) {
      if(Cohort$Name == myCohorts[[i]]$Name) next   ## Don't include the same cohort again
      if(Cohort$Level != myCohorts[[i]]$Level) next ## Only include cohorts from the same level
      
      newDf <- getItemsFromSubjects(myCohorts[[i]],scalarItems=c("StudyID","TreatmentIndex","Treatment"),
                                    longitudinalItems = c("Data","SampleAge"),prevTreatment = TRUE)
  
      myDf <- merge(myDf,newDf, all=TRUE)
    }
  }
  
  
  ## Rename, select and order columns
  if(length(names(myDf)[grep(names(myDf),pattern = "PTRT")]) >0) {
    myDf <- myDf %>% 
      rename(ID=StudyID,DATA=Data,AGE=SampleAge,TRT=TreatmentIndex,TRTS=Treatment) %>% 
      select(ID,DATA,AGE,TRT,TRTS,one_of(StudyObj$StudyDesignSettings$Covariates),Level,CohortName,matches("PTRT"))
  } else {
    myDf <- myDf %>% 
      rename(ID=StudyID,DATA=Data,AGE=SampleAge,TRT=TreatmentIndex,TRTS=Treatment) %>% 
      select(ID,DATA,AGE,TRT,TRTS,one_of(StudyObj$StudyDesignSettings$Covariates),Level,CohortName)
  }

  return (myDf)  
}

##

UpdateProbabilities<-function(Cohort,StudyObj,cohortindex=NULL) {
  DebugPrint(paste0("Doing an analysis based on data in cohort ",Cohort$Name," at time ",StudyObj$CurrentTime),1,StudyObj)
  
  df<-GetCohortData(Cohort,StudyObj,accumulatedData = StudyObj$StudyDesignSettings$AccumulatedData) #Get Cohorts data up to this point in time
  df[df==-99]<-NA #Set -99 to missing
  df<-ImputeCovariates(df,StudyObj,method=StudyObj$StudyDesignSettings$ImpMethod) #Impute missgin covariates
  

  
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
  
  #### Perform LME estimation based on some covariates and treatment effects for each cohort
  if(length(ptrti)>0) {
    
    lmeFormula <- paste0("DATA~1 + AGE + AGE:TRT + (AGE|ID) +",paste0(c(StudyObj$StudyDesignSettings$Covariates,myPTRTs),collapse = " + "))
  } else {
    lmeFormula <- paste0("DATA~1 + AGE + AGE:TRT + (AGE|ID) +",paste0(StudyObj$StudyDesignSettings$Covariates,collapse = " + "))
  }
   
  lmefit <- lmer(lmeFormula,data=df,REML=FALSE) # IIV on baseline only
  
  #lmefit <- lmer(paste0("DATA~1 + AGE + AGE:TRT + (AGE|ID) +",paste0(StudyObj$StudyDesignSettings$Covariates,collapse = " + ")),data=df,REML=FALSE) # IIV on baseline only
  
  #lmefit <- lmer(paste0("DATA~1 + AGE + AGE:TRT + (1+AGE|ID) +",paste0(StudyObj$StudyDesignSettings$Covariates,collapse = " + ")),data=df,REML=FALSE)
  ##### Calculate new probabilites based on another cohort LME results
  lmecoef<-summary(lmefit)$coefficients[,1] #Get coefficicents from LME
  lmese<-summary(lmefit)$coefficients[,2] #Get SE from LME
  lmecoef<-lmecoef[regexpr('AGE:TRT.*',names(lmecoef))==1]
  lmese<-lmese[regexpr('AGE:TRT.*',names(lmese))==1]
  print(lmefit)
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
  #DebugPrint(paste0("Estimated treatment effect in cohort ",Cohort$Name," at time ",StudyObj$CurrentTime),1,StudyObj)
  #DebugPrint(lmecoef,1,StudyObj)
  
  probs <- GetNewRandomizationProbabilities(trtcoeff=lmecoef,trtse=lmese,
                                            StudyObj$StudyDesignSettings$iNumPosteriorSamples) #Calculate randomization probs based on posterior distribution
  
  
  ## Apply any probabiity temperation function, e.g. sqrt
  probs <- StudyObj$StudyDesignSettings$probTemperation(probs)
  
  nonupdateprobs<-probs #Save the not updated probs for statistics
  
  ### Futility - returns prob 0 for futile treatments
  if(StudyObj$StudyDesignSettings$CheckFutility == "before") {
    probs <- StudyObj$Futilityfunction(probs,Cohort,StudyObj)
  }
  
  ## Adjust the probabilities so that minimum allocation is honored
  probs <- updateProbs(StudyObj,probs,Cohort)
  
  ### Futility - returns prob 0 for futile treatments
  if(StudyObj$StudyDesignSettings$CheckFutility == "after") {
    probs <- StudyObj$Futilityfunction(probs,Cohort,StudyObj)
  }
  
  Cohort$UpdateProbabilities<-probs #The latest probability updates
  Cohort$UnWeightedUpdateProbabilities<-nonupdateprobs #The latest probability updates
  
  Cohort$AnalysisTime <- c(Cohort$AnalysisTime,StudyObj$CurrentTime) # Save the study time of analysis
  
  Cohort$UpdateCoefficients<-lmecoef #The latest coefficients
  Cohort$UpdateSE<-lmese #The latest coefficients standard errors
  StudyObj$CohortList[[cohortindex]]<-Cohort #Save the updated cohort
  
  for (j in 1:length(StudyObj$CohortList)) {#Update all dependent cohorts
    if (!is.null(StudyObj$CohortList[[j]]$ProbabilityCohort) && cohortindex!=j && cohortindex==StudyObj$CohortList[[j]]$ProbabilityCohort) {#If cohort j should be updated based on prob in cohort i
      DebugPrint(paste0("Updating probabilities in cohort ",StudyObj$CohortList[[j]]$Name," based on analysis in cohort ",Cohort$Name," at time ",StudyObj$CurrentTime),1,StudyObj)
      RandProbs<-list() #Save previous randomization probabilities on cohort
      RandProbs$CohortTime<-StudyObj$CohortList[[j]]$CurrentTime #The time until the probability was valid
      RandProbs$StudyTime<-StudyObj$CurrentTime
      RandProbs$FromCohort<-cohortindex
      RandProbs$RandomizationProbabilities<-StudyObj$CohortList[[j]]$RandomizationProbabilities
      RandProbs$UnWeightedRandomizationProbabilities<-StudyObj$CohortList[[j]]$UnWeightedRandomizationProbabilities
      RandProbs$Treatments<-StudyObj$CohortList[[j]]$Treatments
      StudyObj$CohortList[[j]]$PreviousRandomizationProbabilities[[length(StudyObj$CohortList[[j]]$PreviousRandomizationProbabilities)+1]]<-RandProbs
      StudyObj$CohortList[[j]]$RandomizationProbabilities<-probs #Update the probabilities for the child cohort
      StudyObj$CohortList[[j]]$UnWeightedRandomizationProbabilities<-nonupdateprobs #The unqeighted update probabilities for the cohort
    }
  }
  return(StudyObj)
}

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
  accumulatedData = TRUE
)

StudyObj <- AdaptiveStudy(StudyObjIni)

plotStudyCohorts(StudyObj)

tmp1 <- getItemsFromSubjects(StudyObj$CohortList[[7]],prevTreatmen=TRUE)
tmp2 <- getItemsFromSubjects(StudyObj$CohortList[[6]],prevTreatmen=TRUE)

tmp3 <- full_join(tmp1,tmp2)
names(tmp2)
