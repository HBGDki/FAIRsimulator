#' DebugPrint
#' @description Function to print debug information
#'
#' @param str A string.
#' @param iLevel The level of debug messages.
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A debug message.
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' }
DebugPrint <- function(str,iLevel,StudyObj) {
  if (iLevel<=StudyObj$DebugLevel) cat(paste0(paste0(rep("\t",iLevel-1),collapse=""), str, "\n"))
}

#' StudyIncrementEvent
#' @description Increments the time in a FAIRsimulator \code{study} object.
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{study} object
#' @export
#'
#' @examples
#' \dontrun{}
StudyIncrementEvent<-function(StudyObj){
  DebugPrint(paste0("Increasing time to: ",StudyObj$CurrentTime+1),3,StudyObj)
  StudyObj$CurrentTime<-StudyObj$CurrentTime+1 #Define the increment function
  StudyObj$CurrentDate<-StudyObj$CurrentDate+1 #The current date
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      if (StudyObj$CohortList[[i]]$Active) { #Only update cohort information if active
        StudyObj$CohortList[[i]]$CurrentTime<-StudyObj$CohortList[[i]]$CurrentTime+1
        ### Update time
        StudyObj$CohortList[[i]]$SubjectList<-lapply(StudyObj$CohortList[[i]]$SubjectList,
                                                     FUN=function(Subject){
                                                       if (Subject$Status==1) {#If an active subject
                                                         Subject$CurrentAge<-Subject$CurrentAge+1
                                                         Subject$CurrentCohortTime<-Subject$CurrentCohortTime+1
                                                         Subject$SubjectCohortTime<-Subject$SubjectCohortTime+1
                                                       }
                                                       return(Subject)
                                                     })
      }
    }
  }
  return(StudyObj)
}

#' CohortCompleted
#' @description Function to check if a cohort is completed.
#' @param Cohort A FAIRsimulator \code{cohort} object.
#' @param StudyObj A FAIRsimulator \code{study} object.
#'
#' @return TRUE or FALSE
#' @export
#'
#' @examples
#' \dontrun{}
CohortCompleted<-function(Cohort,StudyObj) { 
  DebugPrint(paste0("Check if cohort ",Cohort$Name," is completed at time: ",StudyObj$CurrentTime),3,StudyObj)
  Completed<-FALSE
  
  ###All subjects are recruited (now assuming that this will never happen with 2nd statement) && all subjects have dropped out or have had the last samplingtime
  ###Cohort$MaxNumberOfSubjects==Cohort$NumberOfRecruitedSubjects &&

  if (length(Cohort$SubjectList)!=0 && all(lapply(Cohort$SubjectList,function(Subject,LastSample){
    return(Subject$Status==0 || (LastSample %in% Subject$SubjectSampleTime))
  },max(Cohort$SamplingDesign))==TRUE)) {
    DebugPrint(paste0("Cohort ",Cohort$Name," is completed at time: ",StudyObj$CurrentTime),1,StudyObj)
    Completed<-TRUE
  }
  return(Completed)
}



#' InterimAnalyzesTime
#' @description Check if it istime to do an interim analysis for a specific cohrt
#' @param Cohort A FAIRsimulator \code{cohort} object 
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return TRUE or FALSE
#' @export
#'
#' @examples
#' \dontrun{}
InterimAnalyzesTime<-function(Cohort,StudyObj) {
  DebugPrint(paste0("Check if cohort ",Cohort$Name," is about to have an interim analysis",StudyObj$CurrentTime),3,StudyObj)
  TimeToPerformInterim<-FALSE
  tmp<-sum(unlist(lapply(Cohort$SubjectList,function(Subject,LastSample){
    return(as.numeric(Subject$Status==0 || (LastSample %in% Subject$SubjectSampleTime)))},max(Cohort$SamplingDesign))))/Cohort$MaxNumberOfSubjects
  #print(paste0("Percent completed in Cohort ",Cohort$Name,": ",round(tmp*100,1)))
  if (length(tmp)==0) tmp<-0
  
  # if (Cohort$RandomizationAgeRange[1]==6*30 && StudyObj$CurrentTime==6*30-1) {
  #    DebugPrint(paste0("Time to do an interim analyses for cohort ",Cohort$Name," at time: ",StudyObj$CurrentTime," (",round(tmp*100,1)," % subjects completed)"),1,StudyObj)
  #    TimeToPerformInterim<-TRUE
  # }
  
  # if (tmp>3/5 || (Cohort$CurrentTime+round(Cohort$RandomizationAgeRange[1]))==20*30) {
  #  # DebugPrint(paste0("Time to do an interim analyses for cohort ",Cohort$Name," at time: ",StudyObj$CurrentTime," (",round(tmp*100,1)," % subjects completed)"),1,StudyObj)
  #   TimeToPerformInterim<-TRUE
  # }
  # if (Cohort$CurrentTime==(6*30+1)) {
  #   DebugPrint(paste0("Time to do an interim analyses for cohort ",Cohort$Name," at time: ",StudyObj$CurrentTime," (",round(tmp*100,1)," % subjects completed)"),1,StudyObj)
  #   TimeToPerformInterim<-TRUE
  # }
  return(TimeToPerformInterim)
}


#' ImputeCovariates
#' @description Single imputation of covariates by two methods: predictive mean matching and median imputation. Predictive mean matching (default) aggregates the variable to the child level, imputes the covariates by predictive mean matching, and expands the result to all observations per individual (\code{strIDVariable}). This method accounts for the correlations between the variables. Median imputation imputes the median imputation per individual of data frame (without correlation) of covariates in covariates vector.
#'
#' @param df A \code{data.frame} in which to impute missing covariate values.
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param strIDVariable A string with the name of the columns in \code{df} containing the subject identifiers.
#' @param method A string identifying the imputation method. Current implemented methods are 
#' \code{method = "pmm"} (default) and \code{method = "median"}.
#' 
#' @return A \code{data.frame} with all missing covariates imputed.
#' @export
#'
#' @examples
#' \dontrun{}
ImputeCovariates <- function(df, StudyObj, strIDVariable = "ID", 
                             method = c("pmm", "median")) {

  method <- match.arg(method)
  
  if (method == "median") {
    for (i in 1:length(StudyObj$StudyDesignSettings$Covariates)) {
      strcov<-StudyObj$StudyDesignSettings$Covariates[i]
      df[is.na(df[[strcov]]),strcov]<-median(df[!duplicated(strIDVariable),strcov],na.rm=TRUE)
    }
  }
  
  if (method == "pmm") {
    ini <- mice::mice(df, maxit = 0, m = 1, printFlag = FALSE)
    pred <- ini$predictorMatrix
    strIDVariable <- strIDVariable[1]
    pred[pred[, strIDVariable] == 1, strIDVariable] <- (-2)
    meth <- ini$method
    meth[ini$nmis > 0] <- "2lonly.pmm"
    imp <- mice::mice(df, pred = pred, meth = meth, m = 1, printFlag = FALSE)
    df <- mice::complete(imp, 1)
  }
  
  return(df)
}


#' GetNewRandomizationProbabilities
#' @description Function for calculating best treatment probabilitites based on sampling from estimated treatments effect with SEs. 
#' Can also calculate best treatment based on random_sampling (randsam) variable with dim of NumTreatments x iNumPosteriorSamples.
#' @param trtcoeff Treatment coefficients
#' @param trtse Standard errors of the treatment coefficients
#' @param iNumPosteriorSamples Number of samples to use for computing the posterior probability
#' @param randsam An object
#'
#' @return A vector of the same length as the \code{trtcoeff} with the updated posterior probabilities for each treatment to be the best treatment.
#' @export
#'
#' @examples
#' \dontrun{}
GetNewRandomizationProbabilities<- function(trtcoeff,trtse,iNumPosteriorSamples,randsam=NULL) {
  
  sums<-rep(0,length(trtcoeff)+1) #Vector where each treatment is sampled to be best (including SoC)
  
  if (is.null(randsam)) {
    randsam<-matrix(0,length(trtcoeff),iNumPosteriorSamples)
    for (i in 1:length(trtcoeff)) {
      randsam[i,]<-rnorm(iNumPosteriorSamples,mean=trtcoeff[i],sd=trtse[i])
    }
  }
  
  for (i in 1:iNumPosteriorSamples){
    sums[1]<-sums[1]+sum(all(randsam[,i]<=0))
    for (j in 1:length(trtcoeff)) {
      sums[j+1]<-sums[j+1]+sum(max(randsam[,i])==randsam[j,i] & randsam[j,i]>0)
    }
  }
  
  return(sums/iNumPosteriorSamples)
}

#' GetCohortData
#'
#' @description Get dataset for a Cohort object
#'
#' @param Cohort A FAIRsimulator \code{cohort} object
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param accumulatedData If using one cohort (FALSE) or all cohorts with the same level (age) (TRUE)
#'
#' @return A \code{data.frame} with the Cohort level subject specific data (and covariates)
#' @export
#'
#' @examples
#' \dontrun{}
GetCohortData<-function(Cohort,StudyObj, accumulatedData = FALSE) {
  DebugPrint(paste0("Assembly cohort ",Cohort$Name," data available at time ",StudyObj$CurrentTime),2,StudyObj)
  
  if(!accumulatedData) {
    ## Get the data from the current cohort
    myDf <- getItemsFromSubjects(Cohort,
                                 scalarItems=c("StudyID","TreatmentIndex","Treatment"),
                                 longitudinalItems = c("Data","SampleAge"))

    
    if(is.null(myDf)) return(NULL) # No subjects in cohort

  } else {
    
    # Loop over the other cohorts and get the data from those of the same level
    myDf <- getItemsFromSubjects(Cohort,
                                 scalarItems=c("StudyID","TreatmentIndex","Treatment"),
                                 longitudinalItems = c("Data","SampleAge"),
                                 prevTreatment = TRUE)
    

    if(is.null(myDf)) return(NULL) # No subjects in cohort
   
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
      dplyr::rename(ID=StudyID,DATA=Data,AGE=SampleAge,TRT=TreatmentIndex,TRTS=Treatment) %>% 
      select(ID,DATA,AGE,TRT,TRTS,one_of(StudyObj$StudyDesignSettings$Covariates),Level,CohortName,matches("PTRT"))
    for (cstr in names(myDf)[grep(names(myDf),pattern = "PTRT")]) { #Set no previous treatment to "0"
      myDf[,cstr]<-ifelse(is.na(myDf[,cstr]),0,myDf[,cstr])
    }
  } else {
    myDf <- myDf %>% 
      dplyr::rename(ID=StudyID,DATA=Data,AGE=SampleAge,TRT=TreatmentIndex,TRTS=Treatment) %>% 
      select(ID,DATA,AGE,TRT,TRTS,one_of(StudyObj$StudyDesignSettings$Covariates),Level,CohortName)
  }
  DebugPrint(paste0("Cohort data available in ",Cohort$Name, " at time ",StudyObj$CurrentTime," with ",nrow(myDf[!duplicated(myDf$ID),])," IDs."),2,StudyObj)
  return (myDf)  
}



#' UpdateProbabilities
#'
#' @description Update probabilities based on data analysis and probability of being best
#' @param Cohort A FAIRsimulator \code{cohort} object
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param cohortindex Cohort index.
#'
#' @return A FAIRsimulator \code{cohort} object
#' @export
#'
#' @examples
#' \dontrun{}
UpdateProbabilities<-function(Cohort,StudyObj,cohortindex=NULL) {
  DebugPrint(paste0("Doing an analysis based on data in cohort ",Cohort$Name," at time ",StudyObj$CurrentTime),1,StudyObj)
  
  df<-GetCohortData(Cohort,StudyObj,accumulatedData = StudyObj$StudyDesignSettings$AccumulatedData) #Get Cohorts data up to this point in time
  
  if(is.null(df) || nrow(df)==0) return(StudyObj) # No subjects in cohort

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
  
  #### Perform LME estimation based on some covariates and treatment effects for each cohort
  if(length(ptrti)>0) {
    
    lmeFormula <- paste0("DATA~1 + AGE + AGE:TRT + (AGE|ID) +",paste0(c(StudyObj$StudyDesignSettings$Covariates,myPTRTs),collapse = " + "))
    lmFormula <- paste0("DATA~1 + AGE + AGE:TRT + ",paste0(c(StudyObj$StudyDesignSettings$Covariates,myPTRTs),collapse = " + "))
  } else {
    lmeFormula <- paste0("DATA~1 + AGE + AGE:TRT + (AGE|ID) +",paste0(StudyObj$StudyDesignSettings$Covariates,collapse = " + "))
    lmFormula <- paste0("DATA~1 + AGE + AGE:TRT +",paste0(StudyObj$StudyDesignSettings$Covariates,collapse = " + "))
  }
  
  errFun <- function(err) {
    DebugPrint(paste0("Cannot do an LME analysis, doing an LM analysis instead at time ",StudyObj$CurrentTime),1,StudyObj)
    DebugPrint(paste0("Error from LME analysis: ",err),1,StudyObj)
    lmfit <- lm(lmFormula,data=df)
    return(lmfit)
  }
  
  lmefit <- tryCatch({lmer(lmeFormula,data=df,REML=FALSE)},error = errFun)
  
  
   # IIV on baseline only
  
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
    tmpprobs<-probs
    tmp<-StudyObj$Futilityfunction(probs,Cohort,StudyObj)
    probs <-tmp[[1]]
    if (tmp[[2]]==TRUE) {
      Futility<-list()
      Futility$PreviousProbability<-tmpprobs
      Futility$FutileProbability<-probs
      Futility$StudyTime<-StudyObj$CurrentTime
      Futility$CohortTime<-Cohort$CurrentTime
      Futility$CohortName<-Cohort$Name
      StudyObj$FutilityList[[length(StudyObj$FutilityList)+1]]<-Futility
    }
  }
  
  ## Adjust the probabilities so that minimum allocation is honored
  probs <- updateProbs(StudyObj,probs,Cohort)
  
  ### Futility - returns prob 0 for futile treatments
  if(StudyObj$StudyDesignSettings$CheckFutility == "after") {
    tmpprobs<-probs
    tmp<-StudyObj$Futilityfunction(probs,Cohort,StudyObj)
    probs <-tmp[[1]]
    if (tmp[[2]]==TRUE) {
      Futility<-list()
      Futility$PreviousProbability<-tmpprobs
      Futility$FutileProbability<-probs
      Futility$StudyTime<-StudyObj$CurrentTime
      Futility$CohortTime<-Cohort$CurrentTime
      Futility$CohortName<-Cohort$Name
      StudyObj$FutilityList[[length(StudyObj$FutilityList)+1]]<-Futility
    }
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

#' AnalyzeDataEvent
#' @description An event to check if the cohort is completed and also set that subjects are completed according to the criteria set.
#' @param StudyObj A FAIRsimulator \code{cohort} object
#'
#' @return A FAIRsimulator \code{cohort} object
#' @export
#'
#' @examples
#' \dontrun{}
AnalyzeDataEvent<-function(StudyObj) { 
  DebugPrint(paste0("Check if we want to analyze data at study time: ",StudyObj$CurrentTime),3,StudyObj)
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      Cohort<-StudyObj$CohortList[[i]]
      if (Cohort$Active) { #if this is an active cohort
        cohortCompleted<-CohortCompleted(Cohort,StudyObj)
        cohortInterimAnalysis<-StudyObj$StudyDesignSettings$InterimAnalyzesTime(Cohort,StudyObj)
        if (cohortCompleted || cohortInterimAnalysis) { #If we should update probabilities of child cohorts
          StudyObj<-StudyObj$StudyDesignSettings$UpdateProbabilities(Cohort,StudyObj,i)
          #StudyObj<-UpdateProbabilities(Cohort,StudyObj,i)
        }
        if (cohortCompleted) { #If this is a completed cohort
          StudyObj$CohortList[[i]]$Active<-FALSE
        }
      }
    }
  }
  return(StudyObj)  
}


#' StopEvent
#' @description The stop criteria  for the whole study
#' @param StudyObj A FAIRsimulator \code{cohort} object
#'
#' @return The stop criteria for the study.
#' @export
#'
#' @examples
#' \dontrun{}
StopEvent<-function(StudyObj) {
  Crit<-StudyObj$CurrentTime>=StudyObj$StudyDesignSettings$StudyStopTime
  if (Crit) DebugPrint(paste0("Ending study simulations at study time ",StudyObj$CurrentTime),1,StudyObj)
  return(Crit)
}

#' RecruitmentRatefunction
#'
#' @description Function to give back the recruitment rate from a FAIRsimulator \code{study} object.
#' @param StudyObj A FAIRsimulator \code{cohort} object
#' @param Cohort 
#'
#' @return The target number of subjects to recruit per day.
#' @export
#'
#' @examples
#' \dontrun{}
RecruitmentRatefunction<-function(StudyObj,Cohort) {
  #  return(5000) #Instantaneous randomization
  return(20) #20 randomized subjects per time unit
}

#' NewCohort
#'
#' @description Create a newFAIRsimulator \code{cohort} object in a FAIRsimulator \code{study} object, or update an existing one.
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param CohortNum The cohort index number to update.
#'
#' @return A FAIRsimulator \code{cohort} object
#' @export
#'
#' @examples
#' \dontrun{}
NewCohort<-function(StudyObj,CohortNum=NULL) { 
  Cohort<-list()
  Cohort$NumberOfRecruitedSubjects<-0 #The number of already recruited subjects
  Cohort$Active<-TRUE #If the cohort is still active
  Cohort$CohortStartTime<-StudyObj$CurrentTime #The Cohort start time in relation to the study time
  Cohort$CohortStartDate<-StudyObj$CurrentDate #The Cohort start date
  Cohort$CurrentTime<-0 #The Cohort time
  Cohort$AnalysisTime <- NULL # Study time of analyses
  Cohort$RecruitmentRateFunction<-StudyObj$Recruitmentfunction # The recruitment rate function

  if (!is.null(CohortNum)) {
    Cohort$MaxNumberOfSubjects        <- StudyObj$StudyDesignSettings$MaxNumberofSubjects[[CohortNum]] #The maximum number of subject in this cohort
    Cohort$SamplingDesign             <- StudyObj$StudyDesignSettings$SamplingDesigns[[CohortNum]] #The sampling design for this cohort
    Cohort$RandomizationProbabilities <- StudyObj$StudyDesignSettings$RandomizationProbabilities[[CohortNum]] #Get the randomization probabilities
    Cohort$UnWeightedRandomizationProbabilities<-Cohort$RandomizationProbabilities
    Cohort$UpdateProbabilities        <- Cohort$RandomizationProbabilities #The Current updated probabilities
    Cohort$MinAllocationProbabilities <- StudyObj$StudyDesignSettings$MinAllocationProbabilities[[CohortNum]] #Get the minimum randomization probabilities
    Cohort$Treatments                 <- StudyObj$StudyDesignSettings$Treatments[[CohortNum]] #Get the treatments (treatment codes)
    Cohort$EffSizes                   <- StudyObj$StudyDesignSettings$EffSizes[[CohortNum]] #Get the effect sizes for each treatment
    Cohort$RandomizationAgeRange      <- StudyObj$StudyDesignSettings$CohortAgeRange[[CohortNum]] #The age range at randomization for this cohort
    Cohort$DropoutRate                <- StudyObj$StudyDesignSettings$CohortDropoutRate[CohortNum] #The dropout rate for this cohort
    Cohort$NewCohortLink              <- StudyObj$StudyDesignSettings$NewCohortLink[[CohortNum]] #The link to another cohort to get the randomization probabilities
    Cohort$Name                       <- paste0("C-",CohortNum,"-1 [",Cohort$RandomizationAgeRange[1]/30,"-",Cohort$RandomizationAgeRange[2]/30,"m @ rand]")
    Cohort$StartNum                   <- CohortNum #The cohort start number
    Cohort$CycleNum                   <- 1 #The cycle number
    Cohort$Level                      <- which(unlist(StudyObj$StudyDesignSettings$CohortAgeRange)[seq(1,length(StudyObj$StudyDesignSettings$CohortAgeRange)*2,by=2)]==Cohort$RandomizationAgeRange[1])#The relative age level of this cohort, youngest = 1, oldest = length(StudyDesignSettings$CohortAgeRange)
  }
  class(Cohort)<-"cohort"
  return(Cohort)  
}

#' GetRecruitedSubjects
#' 
#' @description Get the number of recruited subjects at current time.
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param Cohort A FAIRsimulator \code{cohort} object
#'
#' @return The number of subjects
#' @export
#'
#' @examples
GetRecruitedSubjects<-function(StudyObj,Cohort) { 
  pcount <- rpois(n=1, lambda=Cohort$RecruitmentRateFunction(StudyObj,Cohort)) #Recruit according to a Poisson distribution
  pcount<-min(pcount,Cohort$MaxNumberOfSubject-Cohort$NumberOfRecruitedSubjects) #Make sure we're not recruting too many subjects
  return(pcount)
}

#' AddCohortEvent
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{study} object
#' @export
#'
#' @examples
AddCohortEvent<-function(StudyObj) { #Initialize a new cohort and add it to the cohort list
  DebugPrint(paste0("Check for new cohorts at study time: ",StudyObj$CurrentTime),3,StudyObj)
  for (i in 1:length(StudyObj$StudyDesignSettings$CohortStartTimes)) {
    if (StudyObj$StudyDesignSettings$CohortStartTimes[i]==StudyObj$CurrentTime) {
      StudyObj$CohortList[[length(StudyObj$CohortList)+1]]<-NewCohort(StudyObj,i) #Add new cohort
      DebugPrint(paste0("Add cohort ",StudyObj$CohortList[[length(StudyObj$CohortList)]]$Name," at study time: ",StudyObj$CurrentTime),1,StudyObj)
    }
  }
  return(StudyObj)
}

#' GetTreatmentRandomizations
#'
#' @description Get the treatment randomizations
#'
#' @param Num The treatment randomization probabilities.
#' @param Cohort A FAIRsimulator \code{cohort} object
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return The treatment index, i.e. the treatment the subjects are randomized to.
#' @export
#'
#' @examples
GetTreatmentRandomizations<-function(Num,Cohort,StudyObj) { #Simple uniform random allocation between treatments based on randomization probabilities
  DebugPrint(paste0("Allocated ",Num," new subjects to treatments"),3,StudyObj)
  rnd<-runif(Num)
  treatmentIndex<-rep(NA,Num)
  for (i in 1:Num){
    lowprob<-0
    highprob<-0
    if (any(is.nan(Cohort$RandomizationProbabilities))) {

      print(paste0("Cohort ",Cohort$Name))
      print(Cohort$RandomizationProbabilities)
    }
    for (j in 1:length(Cohort$RandomizationProbabilities)) {
      highprob<-highprob+Cohort$RandomizationProbabilities[[j]]
      if (rnd[i]>=lowprob && rnd[i]<highprob) {
        treatmentIndex[i]<-j
        break
      }
      lowprob<-lowprob+Cohort$RandomizationProbabilities[[j]]
    }
  }
  return(treatmentIndex)
}


#' GetAgesAtRandomization
#' 
#' @description Returns the uniform age at randomization given an age range
#'
#' @param NumRecruited The number of ages to return
#' @param Cohort A FAIRsimulator \code{cohort} object
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A vector of uniform ages
#' @export
#'
#' @examples
GetAgesAtRandomization<-function(NumRecruited,Cohort,StudyObj) {
  DebugPrint(paste0("Simulate age at randomization for ",NumRecruited," subjects uniformly within [",Cohort$RandomizationAgeRange[1]/30,"-",Cohort$RandomizationAgeRange[2]/30,"] month"),3,StudyObj)
  return(runif(n=NumRecruited,min = Cohort$RandomizationAgeRange[1],max=Cohort$RandomizationAgeRange[2]))
}

#' DropoutTime
#'
#' @description Constant hazard dropout function
#' @param lambda The lambda value
#' @param randnum An optional random number. Will be generated by the fuction if NULL.
#'
#' @return The dropout day.
#' @export
#'
#' @examples
DropoutTime<-function(lambda,randnum=NULL){
  if (is.null(randnum)) randnum<-runif(1)
  return(-log(randnum)/lambda)
}

#' GetSubject
#' 
#' @description Creates a FAIRsimulator \code{individual} object.
#'
#' @param ID The subject identification number
#' @param TRTIndex The treatment index
#' @param AgeAtRand The age at randomization
#' @param Cohort A FAIRsimulator \code{cohort} object
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{individual} object.
#' @export
#'
#' @examples
GetSubject<-function(ID,TRTIndex,AgeAtRand,Cohort,StudyObj) {
  
  Subject<-list()
  Subject$StudyID<-ID #The study specific ID
  Subject$AgeAtRand<-AgeAtRand #The Age at randomization or at cohort start if evolved cohort
  Subject$DateAtRand<-StudyObj$CurrentDate #The Date at randomization at cohort start if evolved cohort
  Subject$RandStudyTime<-StudyObj$CurrentTime #The Randomization/evolved time time (in study time)
  Subject$RandCohortTime<-Cohort$CurrentTime # The Randomization time/evolved time (in cohort time)
  
  Subject$CurrentAge<-AgeAtRand #The current age at a particular study time
  Subject$CurrentCohortTime<-Subject$RandCohortTime #The current cohort time at a particular study time
  Subject$SubjectCohortTime<-0 #The subject specific cohort time, i.e. how long the subject has been in the cohort (starting at 0)
  
  
  Subject$TreatmentIndex<-TRTIndex #Save the randomized treatment
  Subject$Treatment<-Cohort$Treatments[[TRTIndex]]
  Subject$TreatmentEff<-Cohort$EffSizes[[TRTIndex]]
  
  Subject$RandNum<-runif(1) #Individual uniform random number that could be used for e.g. dropout assessment
  
  
  irow<-sample(x=1:nrow(StudyObj$dfFFEMPool),1,replace=TRUE) #Sample row from pool of individuals
  
  #Add FFEM Coefficients
  dfPool<-StudyObj$dfFFEMPool[irow,]
  Subject$FREMCoeffs<-as.numeric(dfPool[1,names(dfPool)[regexpr('COEFF.*',names(dfPool))==1]]) #Get individual coefficicents
  
  #Add covariates
  dfCovs<-StudyObj$dfSubjPool[irow,]
  Subject$Covariates<-dfCovs[,StudyObj$StudyDesignSettings$Covariates] #Store the covariates
  
  #Add IIV Samples
  iNumOMEGA<-length(names(dfPool)[regexpr('VAR.*',names(dfPool))==1])
  iNumOMDIM<--1/2+sqrt(1/4+2*iNumOMEGA)
  
  OMEGAS<-as.numeric(dfPool[1,names(dfPool)[regexpr('VAR.*',names(dfPool))==1]])
  OM                             <- matrix(0, nrow=iNumOMDIM, ncol=iNumOMDIM) #Define an empty matrix
  OM[lower.tri(OM,diag = TRUE)]  <- OMEGAS #Assign lower triangular + diag
  tOM                            <- t(OM) #Get a transposed matrix
  OM[upper.tri(OM,diag = FALSE)] <- tOM[upper.tri(tOM,diag = FALSE)] #Assign the upper triangular except diag
  Subject$IndSamples<-mvrnorm(n = 1, rep(0,iNumOMDIM), OM, tol = 1e-6, empirical = FALSE, EISPACK = FALSE) #Simulate individual from var-cov matrix
  
  Subject$Status<-1 #0=Dropout, 1=Active, 2=Completed
  Subject$DropoutCohortTime<-NA #Set the dropout time with cohort as reference
  Subject$DropoutSubjectTime<-NA #Set the dropout time with subject cohort time as reference
  Subject$DropoutStudyTime<-NA #Set the dropout time with study time as reference
  
  
  #Subject$PreviousTreatmentList<-
  class(Subject) <- "individual"
  return(Subject)
  
}

#' GetSubjects
#'
#' @description Create a list of subjects
#' @param NumRecruited The number of subjects to recruit.
#' @param TreatmentIndex The treatment index
#' @param RandAges The randomization ages
#' @param Cohort A FAIRsimulator \code{cohort} object
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A list of FAIRsimulator \code{individual} objects.
#' @export
#'
#' @examples
GetSubjects<-function(NumRecruited,TreatmentIndex,RandAges,Cohort,StudyObj) {
  SubjectList<-list()
  ID<-StudyObj$StudyDesignSettings$CurrentID
  if (!is.null(Cohort$SubjectList)) SubjectList<-Cohort$SubjectList
  if (NumRecruited>0) {
    for (i in 1:NumRecruited) {
      SubjectList[[length(SubjectList)+1]]<-GetSubject(ID,TreatmentIndex[i],RandAges[i],Cohort,StudyObj)
      ID<-ID+1
    }
  }
  return(SubjectList)
}


#' RecruitmentEvent
#' 
#' @description Event that recruites subjects to cohort
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{study} object
#' @export
#'
#' @examples
RecruitmentEvent<-function(StudyObj) { #´
  DebugPrint(paste0("Check for new recruitments at study time: ",StudyObj$CurrentTime),3,StudyObj)
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      Cohort<-StudyObj$CohortList[[i]]
      if (is.null(Cohort$ParentCohort) && Cohort$Active &&  Cohort$NumberOfRecruitedSubjects<Cohort$MaxNumberOfSubjects) { #if this is an active cohort who still needs to recruit subjects
        NumRecruited<-GetRecruitedSubjects(StudyObj,Cohort)
        if (NumRecruited>0) {
          DebugPrint(paste0("Recruited ",NumRecruited," subjects in cohort ",Cohort$Name," at study time ",StudyObj$CurrentTime),2,StudyObj)
          Cohort$NumberOfRecruitedSubjects<-Cohort$NumberOfRecruitedSubjects+NumRecruited
          TreatmentIndex<-GetTreatmentRandomizations(NumRecruited,Cohort,StudyObj) #Get the treatment codes for the recruited subjects
          RandAges<-GetAgesAtRandomization(NumRecruited,Cohort,StudyObj) #Get the randomization age range of the subjects to recruit
          Cohort$SubjectList<-GetSubjects(NumRecruited,TreatmentIndex,RandAges,Cohort,StudyObj) #Get the subjects
          StudyObj$StudyDesignSettings$CurrentID <- StudyObj$StudyDesignSettings$CurrentID + NumRecruited #Updated the global ID counter with the numbers recruited
          StudyObj$CohortList[[i]]<-Cohort
        }
        if (Cohort$NumberOfRecruitedSubjects==Cohort$MaxNumberOfSubjects) DebugPrint(paste0("All (",Cohort$MaxNumberOfSubjects,") subjects are recruited for cohort ",Cohort$Name," at study time ",StudyObj$CurrentTime),1,StudyObj)
      }
    }
  }
  return(StudyObj)  
}

#' AddEffect
#' 
#' @description Add a treatment effect to a HAZ value
#'
#' @param y The HAZ value to add a treatment effect to
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param Subject A FAIRsimulator \code{individual} object.
#' @param time The time of the HAZ observation
#' @param effecttime The effecttime
#' @param cumeffect Should the effect be cumulative to previous cohorts.
#'
#' @return
#' @export
#'
#' @examples
AddEffect<-function(y,StudyObj,Subject,time,effecttime=6*30,cumeffect=NULL) {
  DebugPrint(paste0("Adding HAZ effect for subject ",Subject$StudyID," at sample time = ",time,""),4,StudyObj)
  cumeff<-0
  #If we have a cumulative effect from previous treatments, i.e. at time=0
  if (!is.null(Subject$CumulativeEffect)) cumeff<-Subject$CumulativeEffect
  if (!is.null(cumeffect)) cumeff<-cumeffect #Use this as previous effect instead
  return(y+Subject$TreatmentEff/effecttime*time + cumeff) #Assume linear effect (and additive)
}

#' SimulateHAZ
#'
#' @description Simulate HAZ values
#' 
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param Subject A FAIRsimulator \code{individual} object.
#' @param time The time of the new HAZ value
#' @param age The age of the subject
#'
#' @return The simulated HAZ value
#' @export
#'
#' @examples
SimulateHAZ<-function(StudyObj,Subject,time,age) {
  DebugPrint(paste0("Simulate HAZ for subject ",Subject$StudyID," at study time: ",StudyObj$CurrentTime," (sample time = ",time,")"),4,StudyObj)
  res_err<-mvrnorm(n = length(age), 0, StudyObj$sig , tol = 1e-6, empirical = FALSE, EISPACK = FALSE) #Simulate residual error on HAZ
  y <-  StudyObj$calcHAZVector(age = age/(30*12),basethetas = StudyObj$thbasevector,covthetas = Subject$FREMCoeff,etas = Subject$IndSamples) #Calculate Y without residual error
  yres<-y+res_err #Add residual error
  yres<-AddEffect(yres,StudyObj,Subject,time) #Add HAZ Effect
  DebugPrint(paste0("HAZ observation ",yres," simulated for subject with age ",age),4,StudyObj)
  return(yres)
}


#' SimulateDataEvent
#' 
#' @description Simulate data for all subjects which should have a sample at this date
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{study} object
#' @export
#'
#' @examples
SimulateDataEvent<-function(StudyObj) { #
  DebugPrint(paste0("Check for data to simulate at study time: ",StudyObj$CurrentTime),3,StudyObj)
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      Cohort<-StudyObj$CohortList[[i]]
      if (Cohort$Active) { #if this is an active cohort
        Cohort$SubjectList<-lapply(Cohort$SubjectList,FUN=function(Subject,Cohort,StudyObj){
          SamplingDesign<-Cohort$SamplingDesign #Might switch to individual sample time
          if (Subject$Status==1 && any(SamplingDesign==Subject$SubjectCohortTime)) {
            #Simulate data 
            Subject$Data<-c(Subject$Data,SimulateHAZ(StudyObj,Subject,
                                                     SamplingDesign[SamplingDesign %in% Subject$SubjectCohortTime],
                                                     Subject$CurrentAge))
            #Add the sample ages
            Subject$SampleAge<-c(Subject$SampleAge,Subject$CurrentAge)
            #Add the Subject specific Sample Time
            Subject$SubjectSampleTime<-c(Subject$SubjectSampleTime,Subject$SubjectCohortTime)
            #Add the Cohort specific sample Time
            Subject$CohortSampleTime<-c(Subject$CohortSampleTime,Subject$CurrentCohortTime)
            #Add the Study specific Sample Time
            Subject$StudySampleTime<-c(Subject$StudySampleTime,StudyObj$CurrentTime)
            if (max(SamplingDesign)==Subject$SubjectCohortTime) {
              DebugPrint(paste0("Subject ",Subject$StudyID," in cohort ",Cohort$Name," is has completed the cohort at time: ",StudyObj$CurrentTime),3,StudyObj)
              Subject$Status<-2 #Completed Subject, should be rerandomized
            }
          }
          return(Subject)  
        },Cohort,StudyObj)
        StudyObj$CohortList[[i]]<-Cohort  
      }
    }
  }
  return(StudyObj)  
}

#' DropoutEvent
#' 
#' @description Event that check if subject have dropped out
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{study} object
#' @export
#'
#' @examples
DropoutEvent<-function(StudyObj) { #
  DebugPrint(paste0("Check for subject dropout at study time: ",StudyObj$CurrentTime),3,StudyObj)
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      Cohort<-StudyObj$CohortList[[i]]
      if (Cohort$Active) { #if this is an active cohort
        Cohort$SubjectList<-lapply(Cohort$SubjectList,FUN=function(Subject,Cohort,StudyObj){
          if (Subject$Status==1) {## If Active subject
            dropouttime<-DropoutTime(Cohort$DropoutRate,randnum=Subject$RandNum)
            if (Subject$SubjectCohortTime>=dropouttime) {#If this Subject is dropping out 
              DebugPrint(paste0("Subject ",Subject$StudyID," in cohort ",Cohort$Name," is dropping out at time: ",StudyObj$CurrentTime),2,StudyObj)
              Subject$Status<-0 #Set to dropout
              Subject$DropoutCohortTime<-Subject$CurrentCohortTime #Set the dropout time with cohort as reference
              Subject$DropoutSubjectTime<-Subject$SubjectCohortTime #Set the dropout time with subject cohort time as reference
              Subject$DropoutStudyTime<-StudyObj$CurrentTime #Set the dropout time with study time as reference
            }
          }
          return(Subject)  
        },Cohort,StudyObj)
        StudyObj$CohortList[[i]]<-Cohort
      }
    }
  }
  return(StudyObj)
}

#' MoveSubjects
#' 
#' @description Move subjects from FromCohort to ToCohort which are completed (Status==2), Re-randomize treatments based on ToCohort rand probabilities
#'
#' @param FromCohort A FAIRsimulator \code{cohort} object for the cohort to move subjects from.
#' @param ToCohort A FAIRsimulator \code{cohort} object for the cohort to move subjects to.
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{cohort} object (the ToCohort)
#' @export
#'
#' @examples
MoveSubjects<-function(FromCohort,ToCohort,StudyObj) { #
  DebugPrint(paste0("Check for subjects to move from cohort ",FromCohort$Name," to cohort ",ToCohort$Name," at study time: ",StudyObj$CurrentTime),3,StudyObj)
  fromids<-NULL
  fromstatus<-NULL
  toids<-NULL
  if (!is.null(FromCohort$SubjectList)) {
    fromids<-FromCohort$SubjectList %listmap% "StudyID"
    fromstatus<-FromCohort$SubjectList %listmap% "Status"
  }
  if (!is.null(ToCohort$SubjectList)){
    toids<-ToCohort$SubjectList %listmap% "StudyID"
  }
  if (!is.null(fromids)) {
    for (i in 1:length(fromids)) {
      if (fromstatus[i]==2 && !(fromids[i] %in% toids)) { #If this subject is completed and not already moved
        DebugPrint(paste0("Moving subject ",fromids[i]," from cohort ",FromCohort$Name," to cohort ",ToCohort$Name," at study time: ",StudyObj$CurrentTime),3,StudyObj)
        Subject<-FromCohort$SubjectList[[i]] #Get the subject
        #### Update the subject with the ToChort information......
        Subject$AgeAtRand<-Subject$CurrentAge #The Age at randomization or at cohort start if evolved cohort
        Subject$DateAtRand<-StudyObj$CurrentDate #The Date at randomization at cohort start if evolved cohort
        Subject$RandStudyTime<-StudyObj$CurrentTime #The Randomization/evolved time time (in study time)
        Subject$RandCohortTime<-ToCohort$CurrentTime # The Randomization time/evolved time (in new cohort time)
        
        #Get cumulative effect (dHAZ) at this time (i.e. time when new cohort starts)
        Subject$CumulativeEffect<-AddEffect(0,StudyObj,Subject,Subject$SubjectCohortTime)
        EffectAtAge<-AddEffect(0,StudyObj,Subject,Subject$SubjectCohortTime,cumeffect = 0)
        
        Subject$Status<-1 #Active
        Subject$CurrentCohortTime<-Subject$RandCohortTime #The current cohort time at a particular study time
        Subject$SubjectCohortTime<-0 #The subject specific cohort time, i.e. how long the subject has been in the cohort (starting at 0)
        
        TreatmentIndex<-GetTreatmentRandomizations(1,ToCohort,StudyObj) #Get the new treatment code for the recruited subjects
        if (is.null(Subject$PreviousTreatmentIndex)) {
          Subject$PreviousTreatmentIndex<-list()
          Subject$PreviousTreatment<-list()
          Subject$PreviousTreatmentEff<-list()
          Subject$PreviousTreatmentAge<-list()
          Subject$PreviousTreatmentEffectAtAge<-list()
        }
        
        Subject$PreviousTreatmentIndex[[length(Subject$PreviousTreatmentIndex)+1]]<-Subject$TreatmentIndex #Previous treatment index
        Subject$PreviousTreatment[[length(Subject$PreviousTreatment)+1]]<-Subject$Treatment #Previous treatment
        Subject$PreviousTreatmentEff[[length(Subject$PreviousTreatmentEff)+1]]<-Subject$TreatmentEff #Previous treatment effect
        Subject$PreviousTreatmentAge[[length(Subject$PreviousTreatmentAge)+1]]<-Subject$CurrentAge #Previous treatment age
        Subject$PreviousTreatmentEffectAtAge[[length(Subject$PreviousTreatmentEffectAtAge)+1]]<-EffectAtAge #Previous treatment effect @ age 
        
        Subject$TreatmentIndex<-TreatmentIndex #Save the randomized treatment
        Subject$Treatment<-ToCohort$Treatments[[TreatmentIndex]]
        Subject$TreatmentEff<-ToCohort$EffSizes[[TreatmentIndex]]
        
        
        if (StudyObj$StudyDesignSettings$MoveLastSampleToNewCohort) {#If reuse the last sample as baseline sample in new cohort
          Data<-Subject$Data
          SampleAge<-Subject$SampleAge
          SubjectSampleTime<-Subject$SubjectSampleTime
          CohortSampleTime<-Subject$CohortSampleTime
          StudySampleTime<-Subject$StudySampleTime
        } 
        #Reset all previous data and sample times
        Subject$Data<-NULL
        Subject$SampleAge<-NULL
        Subject$SubjectSampleTime<-NULL
        Subject$CohortSampleTime<-NULL
        Subject$StudySampleTime<-NULL
        
        if (StudyObj$StudyDesignSettings$MoveLastSampleToNewCohort) {#If resuse the last sample as baseline sample in new cohort
          if (length(Data)!=0) {
            Subject$Data<-Data[length(Data)] #Take last sample
            Subject$SampleAge<-SampleAge[length(SampleAge)]
            Subject$SubjectSampleTime<-0
            Subject$CohortSampleTime<-ToCohort$CurrentTime
            Subject$StudySampleTime<-StudyObj$CurrentTime
          }
        }
        
        
        Subject$RandNum<-runif(1) #Individual uniform random number that could be used for e.g. dropout assessment
        ToCohort$NumberOfRecruitedSubjects<-ToCohort$NumberOfRecruitedSubjects+1 #Increase the number of subjects in cohort
        if (is.null(ToCohort$SubjectList)) ToCohort$SubjectList<-list()
        ToCohort$SubjectList[[length(ToCohort$SubjectList)+1]]<-Subject #Add the new subject to the new cohort list
      }
    }
  }
  return(ToCohort)
}   



#' GetStartNumIndex
#'
#' @description Return the cohort list index of a certain cohort start number
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param StartNum The cohort start number
#'
#' @return The cohort list index
#' @export
#'
#' @examples
GetStartNumIndex <- function(StudyObj,StartNum) {
  for (i in 1:length(StudyObj$CohortList)) {
    if (StartNum==StudyObj$CohortList[[i]]$StartNum) return(i)
  }
  return (NA)
}

#' MoveCompletedSubjects
#'
#' @description Event that checks if completed subjects should be moved to new Cohort. Only subjects that are in a cohort that is linked to another cohort is affected.
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{study} object
#' @export
#'
#' @examples
MoveCompletedSubjects<-function(StudyObj) { 
  DebugPrint(paste0("Check for completed subjects at study time: ",StudyObj$CurrentTime),3,StudyObj)
  NewChildCohortList<-list()
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      ChildCohort<-NULL
      Cohort<-StudyObj$CohortList[[i]]
      if (Cohort$Active && !is.null(Cohort$NewCohortLink) && any(lapply(Cohort$SubjectList,function(Subject){return(Subject$Status==2)})==TRUE))  { 
        #if this is an active cohort with completed subjects and it has a linked cohort
        if (is.null(Cohort$ChildCohort)) { #If we need to start a new cohort first
          DebugPrint(paste0("A new child cohort is created based on completed subjects from cohort ",Cohort$Name," at study time ",StudyObj$CurrentTime),1,StudyObj)
          NewChildCohort<-NewCohort(StudyObj,CohortNum=NULL)
          NewCohortLinkIndex<-GetStartNumIndex(StudyObj,Cohort$NewCohortLink)
          NewChildCohort$MaxNumberOfSubjects<-Cohort$MaxNumberOfSubjects #The maximum number of subject in this cohort
          NewChildCohort$SamplingDesign<-StudyObj$CohortList[[NewCohortLinkIndex]]$SamplingDesign #The sampling design for this cohort
          NewChildCohort$RandomizationProbabilities<-StudyObj$CohortList[[NewCohortLinkIndex]]$UpdateProbabilities #Get the randomization probabilities
          NewChildCohort$UnWeightedRandomizationProbabilities<-StudyObj$CohortList[[NewCohortLinkIndex]]$UnWeightedUpdateProbabilities #Get the randomization probabilities
          NewChildCohort$UpdateProbabilities<-NewChildCohort$RandomizationProbabilities
          NewChildCohort$MinAllocationProbabilities<-StudyObj$CohortList[[NewCohortLinkIndex]]$MinAllocationProbabilities #Get the minimum randomization probabilities
          NewChildCohort$Treatments<-StudyObj$CohortList[[NewCohortLinkIndex]]$Treatments #Get the treatments (treatment codes)
          NewChildCohort$EffSizes<-StudyObj$CohortList[[NewCohortLinkIndex]]$EffSizes #Get the effect sizes for each treatment
          NewChildCohort$RandomizationAgeRange<-Cohort$RandomizationAgeRange #The age range at randomization for this cohort
          NewChildCohort$DropoutRate<-StudyObj$CohortList[[NewCohortLinkIndex]]$DropoutRate #The dropout rate for this cohort
          NewChildCohort$NewCohortLink<-StudyObj$CohortList[[NewCohortLinkIndex]]$NewCohortLink #The link to another cohort to get the randomization probabilities
          NewChildCohort$Level<-Cohort$Level+1 #This cohort is moving up one level
          NewChildCohort$Name<-paste0("C-",Cohort$StartNum,"-",Cohort$CycleNum+1," [",Cohort$RandomizationAgeRange[1]/30,"-",Cohort$RandomizationAgeRange[2]/30,"m @ rand]")
          NewChildCohort$StartNum<-Cohort$StartNum #The Cohort starting number
          NewChildCohort$CycleNum<-Cohort$CycleNum+1 #The Cohort cycle number
          StudyObj$CohortList[[i]]$ChildCohort<-length(StudyObj$CohortList)+length(NewChildCohortList)+1 #Add reference to the new cohort as a child cohort
          NewChildCohort$ParentCohort<-i #Add a reference to old Cohort as a parent cohort
          NewChildCohort$ProbabilityCohort<-NewCohortLinkIndex #The cohort where to update the probabilities from

          NewChildCohortList[[length(NewChildCohortList)+1]]<-NewChildCohort
        } else {
          # StudyObj$CohortList[[Cohort$ChildCohort]]<-MoveSubjects(Cohort,StudyObj$CohortList[[Cohort$ChildCohort]],StudyObj)          
          # browser()
        }       
      }
    }
  }
  
  ### Add NewChildCohorts to Cohort list
  if (length(NewChildCohortList)!=0) {
    for (i in 1:length(NewChildCohortList)) {
      StudyObj$CohortList[[(length(StudyObj$CohortList)+1)]]<-NewChildCohortList[[i]]
    }
  }
  
  ### Update probabilities from parallell cohorts first then move subjects
  StudyObj<-UpdateProbabilitiesEvent(StudyObj)
  
  #Move subjects
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      ChildCohort<-NULL
      Cohort<-StudyObj$CohortList[[i]]
      if (Cohort$Active && !is.null(Cohort$NewCohortLink) && any(lapply(Cohort$SubjectList,function(Subject){return(Subject$Status==2)})==TRUE))  {
        #if this is an active cohort with completed subjects and it has a linked cohort
        StudyObj$CohortList[[Cohort$ChildCohort]]<-MoveSubjects(Cohort,StudyObj$CohortList[[Cohort$ChildCohort]],StudyObj)
      }
    }
  }
  return(StudyObj)
}
  
#' UpdateProbabilitiesEvent
#'
#' @description Update probabilities of non-connected cohorts.
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{study} object
#' @export
#'
#' @examples
UpdateProbabilitiesEvent<-function(StudyObj) {
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      linkedcohorts<-c()
      Cohorti<-StudyObj$CohortList[[i]]
      for (j in 1:length(StudyObj$CohortList)) {
        Cohortj<-StudyObj$CohortList[[j]]
        if (i!=j && !is.null(Cohorti$ProbabilityCohort) && !is.null(Cohortj$ProbabilityCohort) && Cohorti$ProbabilityCohort==Cohortj$ProbabilityCohort) { #These two cohorts should be updated based on the same probability
         # if (any(Cohorti$RandomizationProbabilities!=Cohortj$UpdateProbabilities) && Cohorti$CohortStartTime>Cohortj$CohortStartTime && Cohorti$Active==TRUE) {
          if (Cohorti$CohortStartTime>Cohortj$CohortStartTime && Cohorti$Active==TRUE) {
              linkedcohorts<-c(linkedcohorts,j) #Add a linked to cohort to the index list
          }
        }
      }  
      ### Check which one is the latest of the linked cohorts and use that for updating
      if (!is.null(linkedcohorts)) {
        j<-max(linkedcohorts)
        Cohortj<-StudyObj$CohortList[[j]]
        if (!is.null(Cohortj$UpdateProbabilities) && any(Cohorti$RandomizationProbabilities!=Cohortj$UpdateProbabilities)) {
          DebugPrint(paste0("Updating probabilities in ",Cohorti$Name," based on probabilities in ",Cohortj$Name," at time: ",StudyObj$CurrentTime),1,StudyObj)
          RandProbs<-list() #Save previous randomization probabilities on cohort
          RandProbs$CohortTime<-StudyObj$CohortList[[i]]$CurrentTime #The time until the probability was valid
          RandProbs$StudyTime<-StudyObj$CurrentTime
          RandProbs$FromCohort<-j
          RandProbs$RandomizationProbabilities<-StudyObj$CohortList[[i]]$RandomizationProbabilities
          RandProbs$UnWeightedRandomizationProbabilities<-StudyObj$CohortList[[i]]$UnWeightedRandomizationProbabilities
          RandProbs$Treatments<-StudyObj$CohortList[[i]]$Treatments
          StudyObj$CohortList[[i]]$PreviousRandomizationProbabilities[[length(StudyObj$CohortList[[i]]$PreviousRandomizationProbabilities)+1]]<-RandProbs
          StudyObj$CohortList[[i]]$RandomizationProbabilities<-Cohortj$UpdateProbabilities #Update the probabilities for the parallell cohort
          StudyObj$CohortList[[i]]$UnWeightedRandomizationProbabilities<-Cohortj$UnWeightedUpdateProbabilities
        }
      }
    }
  }
  return(StudyObj)
}

#' AddNewBirthCohortEvent
#'
#' @description Function for adding new cohorts at birth based on previous ended birth cohorts
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{study} object
#' @export
#'
#' @examples
AddNewBirthCohortEvent<-function(StudyObj) {
  #Note - assuming birth cohort is first element in StudyDesignSettings
  
  if (!is.null(StudyObj$CohortList)) {
    for (i in 1:length(StudyObj$CohortList)) {
      Cohort<-StudyObj$CohortList[[i]]
      if (Cohort$Active==FALSE && Cohort$RandomizationAgeRange[1]==0 && Cohort$CycleNum==1 && Cohort$CohortStartTime+Cohort$CurrentTime==StudyObj$CurrentTime && StudyObj$CurrentTime<=StudyObj$StudyDesignSettings$LatestTimeForNewBirthCohorts) #If this cohort just ended and its time to add a new cohort
      {
        BirthCohortIndex<-which(unlist(lapply(StudyObj$StudyDesignSettings$CohortAgeRange,function(x) {x[1]==0}))==TRUE) #Get the birth cohort index (i.e. lower AgeRangeAtRandomization=0)
        NewChildCohort<-NewCohort(StudyObj,CohortNum=BirthCohortIndex)
        NewChildCohort$StartNum <- max(StudyObj$CohortList %listmap% "StartNum")+1
        NewChildCohort$Name<-paste0("C-",NewChildCohort$StartNum,"-",1," [",Cohort$RandomizationAgeRange[1]/30,"-",Cohort$RandomizationAgeRange[2]/30,"m @ rand]")
        DebugPrint(paste0("Create new birth cohort ",NewChildCohort$Name," based on probabilities in ",Cohort$Name," at time: ",StudyObj$CurrentTime),1,StudyObj)
        NewChildCohort$RandomizationProbabilities<-Cohort$UpdateProbabilities #Use the latest probabilities from previous birth cohort
        NewChildCohort$UnWeightedRandomizationProbabilities<-Cohort$UnWeightedUpdateProbabilities
        StudyObj$CohortList[[(length(StudyObj$CohortList)+1)]]<-NewChildCohort
        break
      }
    }
  }
  return(StudyObj)
}


#' InitEvent
#'
#' @description Initialization event function. Reads in the necessary FREM stuff for the simulations.
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A FAIRsimulator \code{study} object
#' @export
#'
#' @examples
InitEvent <- function(StudyObj) {

  #Read in all the FREM stuff
  runno <-'71' # The updated FREM model
  myExt <- system.file("extdata",paste0("run",runno,".ext"),package="FAIRsimulator")
  dfext <- subset(getExt(extFile = myExt),ITERATION=="-1000000000") #REad in parameter values

  StudyObj$dfFFEMPool<-read.csv(file=system.file("extdata","dfFFEMPool-India-run71.csv",package="FAIRsimulator")) #Read in models that we can use
  StudyObj$dfSubjPool<-read.csv(file=system.file("extdata","dfSubj-India-run71.csv",package="FAIRsimulator")) #Read in covariates that we can use

  StudyObj$calcHAZVector<-calcHAZVector #The function for simulating HAZ observations

  noBaseThetas <- 6
  noCovThetas  <- 33

  StudyObj$thbasevector <- as.numeric(dfext[2:(noBaseThetas+1)])
  StudyObj$sig     <- as.numeric(dfext[noBaseThetas+2+noCovThetas]) #Get the additive residual error (variance)
  
  return(StudyObj)
}

#' futilityFunction
#' 
#' @description Function to determine if treatments are futile.
#'
#' @param probs The current randomization probabilities after having taken the minimum randomisation probabilities into account.
#' @param Cohort A FAIRsimulator \code{cohort} object
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param minSubj The cutoff for the number of subjects expected to be randomised to a treatment.
#' 
#' @details The function checks the expected number of subjects to be randomised to a particular treatment given the probabilities in \code{probs}. Any treatment that has an expected number of subjects 
#' lower than or equal to \code{minSubj} will have their randomisation probabilitiy set to 0. The remaining probabilities will be adjusted to add up to 1.
#'
#' @return A List with vector of updated probabilities and futility used (TRUE/FALSE)
#' @export
#'
#' @examples
futilityFunction <- function(probs,Cohort,StudyObj,minFutProb=StudyObj$StudyDesignSettings$MinimumFutilityProbability) {
  anyfutile<-any(probs<=minFutProb)
  probs <- updateProbs(StudyObj,ifelse(probs<=minFutProb,0,probs),Cohort)
  return(list(probs),anyfutile)
}

#' updateProbs
#' 
#' @description Updates a vector of probabilties so that one or more of the probabilities are of at least a certain size.
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param probs A vector of probabilities to be adjusted. The vector sum does not have to be 1.
#' @param Cohort A FAIRsimulator \code{cohort} object
#' @param minProb A vector of the same length as \code{probs} with the minimum values.
#'
#' @details Typically used to ensure that standard of care always have a certain randomisation probability.
#' @return A vector of probabilties in which the minimum probabilties are honored and the remaining probabilities are adjusted so that the vector sum is 1.
#' @export
#'
#' @examples
updateProbs <- function(StudyObj,probs,Cohort,minProb = Cohort$MinAllocationProbabilities) {
  
  probs[minProb!=0] <- minProb[minProb!=0]
  probspresum       <- sum(probs[minProb!=0])
  probssum          <- sum(probs[minProb==0])
  probs[minProb==0] <- (1-probspresum)*probs[minProb==0]/probssum
  DebugPrint(paste0("Recalculated randomization probabilities in ",Cohort$Name," at time ",StudyObj$CurrentTime),1,StudyObj)
  
  return(probs)
}

#' getCohortAgeRange
#' 
#' Creates a lookup data frame with CohortName, Recruitment age range and a Label for the age group. 
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param cohortAgeNames An optional vector of labels for the age groups
#'
#' @return A data frame with columns CohortName, AgeRange and CohortAge
#' @export
#'
#' @examples
#' \dontrun{
#' getCohortAgeRange(StudyObj)
#' }
getCohortAgeRange <- function(StudyObj,cohortAgeNames = c("0-6 months","6-12 months","12-18 months")) {
  
  myCohorts <- cohorts(StudyObj)
  
  resDf1 <- data.frame(CohortName = names(myCohorts %listmap% "Level"),
                       Level = myCohorts %listmap% "Level")
  
  resDf1$CohortAge <- cohortAgeNames[resDf1$Level] 
  resDf1$CohortName <- as.character(resDf1$CohortName)
  resDf1$CohortAge <- factor(resDf1$CohortAge,levels=cohortAgeNames)
  
  return(resDf1)
}

#' getProbData
#' 
#' Extracts the randomization probabilities from all cohorts over time.
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param strProb A string indicating the probabilities to extract. Possible values are "UnWeightedRandomizationProbabilities" and "RandomizationProbabilities". The former 
#' are the probabilities after adjustment for minimum allocation probabilities and the latter are the unadjusted probabilities.
#' @param ... Other parameters passed to \code{getCohortAgeRange}
#'
#' @return A data.frame with columns CohortName (cohort name), Prob (probability),  TreatmentName, CohortAge (cohort age group label) and  RandStudyTime 
#' (the time the randomisation probabilities were set).
#' @export
#'
#' @examples
#' \dontrun{
#' probData <- getProbData(StudyObj)
#' }

getProbData <- function(StudyObj,strProb="UnWeightedRandomizationProbabilities",...) {
  
  allData <- getAllSubjectData(StudyObj)
  
  resDf   <- getCohortAgeRange(StudyObj,...)  
  allData <- left_join(allData,resDf,by="CohortName")
  
  myCohorts <- cohorts(StudyObj)
  
  cohrtNams    <- dimnames(myCohorts %listmap% strProb)[[2]]
  probs        <- data.frame(myCohorts %listmap% strProb) 
  names(probs) <- cohrtNams
  
  treatnames        <- data.frame(myCohorts %listmap% "Treatments")
  names(treatnames) <- cohrtNams
  

  probs      <- probs      %>%   gather(key="CohortName",value="Prob") 
  treatnames <- treatnames %>%   gather(key="CohortName",value="TreatmentName") 
  
  probs$TreatmentName <- treatnames$TreatmentName
  
  randProbs <- allData %>% distinct(CohortAge,CohortName,.keep_all=TRUE) %>% dplyr::select(CohortAge,CohortName,CohortStartTime) %>% left_join(probs,.,by="CohortName")
  
  randProbs$TreatmentName <-factor(randProbs$TreatmentName,levels=unique(unlist(StudyObj$StudyDesignSettings$Treatments)))
  
  return(randProbs)
}

#' plotProbs
#' 
#' Plots the randomization probabilities over study time stratified by age treatment group.
#'
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param ... Other parameters passed to \code{getCohortAgeRange}
#'
#' @return An invisible list of ggplot objects, one for each age tretament group.
#' @export
#'
#' @examples
#' \dontrun{
#' plotProbs(StudyObj)
#' }

plotProbs <- function(StudyObj,strProb="UnWeightedRandomizationProbabilities",cohortAgeNames=NULL,...) {
  
  if(is.null(cohortAgeNames)) stop("Need to specify the cohortAgeNames.")
  probData <- getProbData(StudyObj,strProb=strProb,cohortAgeNames=cohortAgeNames)

  xs <- split(probData,f = probData$CohortAge)
  
  plotList <- list()
  

  plotList[[1]] <- ggplot(xs[[1]],aes(CohortStartTime/30,Prob,group=TreatmentName,color=TreatmentName)) +
    geom_point() +
    geom_line() +
    theme(legend.position="top") +
    labs(color=NULL) +
    xlab("Study time (months)") +
    ylab("Randomization probability") +
    facet_wrap(~CohortAge)
  
  if(length(xs) > 1) {
    for(i in 2:length(xs)) {
      plotList[[i]] <- plotList[[1]] %+% xs[[i]]
    }
  }
  
  retVal <- plotList 
  
  ## print the plots  
  plotList$ncol <- length(xs)
  do.call(grid.arrange,plotList)
  
  # Return the plot list without printing
  return(invisible(retVal))
}

#' probTemperation
#'
#' This funciton is called after each update of the probabilities. It is a placeholder for a function that makes the randomization updates less drastic. 
#' The default function does nothing to the probabilities.
#' @param probs The current set of probabilities.
#'
#' @return A set of probabilities
#' @export
#'
#' @examples
#' \dontrun{
#' probTemperation(probs) # Default. Will give back the same probabilities.
#' probTemperation <- function(probs) {
#'      probs <- sqrt(probs)/sum(sqrt(probs))
#'      return(probs)
#'}
#' probTemperation(probs) # Will now give back the less dramatic square rooted probabilities.
#' }
probTemperation <- function(probs) {
  return(probs)
}

#' runMultiSim
#' 
#' Simulate multiple instances of a study, optionally in a parallel fashion.
#'
#' @param StudyOnjIni An initial FAIRsimulator \code{study} object created by \code{createStudy}.
#' @param extractProbs Whether to also extract randomization probabilities (in parallel).
#' @param iter The number of studies to simulate.
#' @param ncores The number of cores to use for parallel execution. 
#' @param ... Other parameters. Typically \code{strProb} to specify "UnWeightedRandomizationProbabilities" instead of the default "RandomizationProbabilities".
#'
#' @return A list with a list if the simulated FAIRsimulator \code{study} objects and, optionally, a data.frame with the randomizaion probabilities.
#' @export
#'
#' @examples
#' \dontrun{
#' StudyObjIni <- createStudy(latestTimeForNewBirthCohorts=18*30,studyStopTime = 32*30,
#' nSubjects = c(320,320,320),
#' randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
#' #minAllocationProbabilities = list(c(0,rep(0,4)),c(0,rep(0,4)),c(0,rep(0,4))),
#' minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
#' treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
#' effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
#' Recruitmentfunction=function(...) {return(5000)},
#' minSubjects = 10)
#'
#' myMultStud <- runMultiSim(StudyOnjIni,iter=10,ncores=2)
#' }

runMultiSim <- function(StudyOnjIni,extractProbs=TRUE,iter=1,ncores=1,strProb="UnWeightedRandomizationProbabilities",cohortAgeNames=NULL,clean=TRUE) {
  
  if(extractProbs & is.null(cohortAgeNames)) stop("Need to specify cohortAgeNames.")

  ## Start the parallell engine if ncores > 1
  if(ncores>1) {registerDoParallel(cores=ncores)}
  
  ## Run the simulations
  myRes <- foreach(i=1:iter) %dopar% AdaptiveStudy(StudyObjIni)
  

  cleanFun <- function(myList) {
    myList$dfSubjPool <- NULL
    myList$dfFFEMPool <- NULL
    if (!is.null(myList))  {
      for (i in 1:length(myList$CohortList)) {
        if (!is.null(myList$CohortList[[i]]$SubjectList)) {
          for (j in 1:length(myList$CohortList[[i]]$SubjectList)) {
            myList$CohortList[[i]]$SubjectList[[j]]$FREMCoeffs<-NULL
            myList$CohortList[[i]]$SubjectList[[j]]$Covariates<-NULL
            myList$CohortList[[i]]$SubjectList[[j]]$IndSamples<-NULL
            myList$CohortList[[i]]$SubjectList[[j]]$SampleAge<-NULL
            myList$CohortList[[i]]$SubjectList[[j]]$Data<-NULL
            myList$CohortList[[i]]$SubjectList[[j]]$StudySampleTime<-NULL
            myList$CohortList[[i]]$SubjectList[[j]]$CohortSampleTime<-NULL
            myList$CohortList[[i]]$SubjectList[[j]]$SubjectSampleTime<-NULL
          }
        }
      }
    }
    return(myList)
  }
  
  ## Extract the probabilities if requested
  if(extractProbs) {
    probDf  <- getMultiProbList(myRes,ncores=ncores,cohortAgeNames=cohortAgeNames,strProb=strProb)
      
    if(clean) myRes <- lapply(myRes,cleanFun)
    retList <- list(studList = myRes,probDf = probDf)
  } else {
    
    if(clean) myRes <- lapply(myRes,cleanFun)

    retList <- list(studList = myRes)
  }
  
  return(retList)
  
}

#' getMultiProbList
#' 
#' Extracts the probability lists from a list of FAIRsimulator \code{study} objects. It will optionally be done in parallel.
#'
#' @param multiStudObj A list of FAIRsimulator \code{study} objects
#' @param ncores The number of cores to use. 
#' @param ... Arguments to \code{getProbData}, in particular \code{strProb}
#'
#' @return A data frame with probabilities. Includes a column with Iteration number.
#' @export
#'
#' @examples
#' \dontrun{
#' getMultiProbList(multiStudObj,ncores=2,strProb="UnWeightedRandomizationProbabilities") 
#' }

getMultiProbList <- function(multiStudObj,ncores=1,strProb="UnWeightedRandomizationProbabilities",cohortAgeNames=NULL) {
  
  if(is.null(cohortAgeNames)) stop("Need to specify cohortAgeNames.")

  
  ## Start the parallell engine if ncores > 1
  if(ncores>1) {registerDoParallel(cores=ncores)}
  
  ## get a list of all probability data

  probList <- foreach(i=1:length(multiStudObj)) %dopar% getProbData(multiStudObj[[i]],cohortAgeNames=cohortAgeNames,strProb=strProb)
  
  ## Create one data frame of all the probability data frames
  probDf <- bind_rows(probList)
  
  ## Add an iteration variable (Iter)
  probDf$Iter <- rep(1:length(multiStudObj),each=nrow(probDf)/length(multiStudObj))
  
  return(probDf)  
}

#' plotMultiProb
#' 
#' Plots the updated probabilities over study time, per cohort, with prediction interval.
#'
#' @param probDf A data.frame with probabilities for multiple FAIRsimulator \code{study} objects.
#' @param ylb y-axis label for 
#' @param pup Upper percentile for the prediction interval.
#' @param pdo Lower percentile for the prediction interval.
#'
#' @return An invisible list of ggplot objects as long as the number of cohorts.
#' @export
#'
#' @examples
#' \dontrun{
#' plotMultiProb(probDf)
#' plotMultiProb(probDfUnweighted,ylb="Unweighted randomization probabilities")
#' }
plotMultiProb <- function(probDf,ylb="Randomization probability",pup = 0.95,pdo=0.05) {
  
  ## Define functions to compute the upper and lower limits of the prediction intervals of the probabilities
  up95 <- function(x,PRup=pup) quantile(x,p=PRup)
  do05 <- function(x,PRdo=pdo) quantile(x,p=PRdo)
  
  ## Compute the plot data
  sumProbData <- probDf %>% 

    group_by(CohortAge,TreatmentName,CohortStartTime) %>% 

    summarise(Mean=mean(Prob),Up=up95(Prob),Down=do05(Prob)) 
  
  ## To get legends for each age cohort we need to create separate graphs for the. We'll split the data and create the graphs and then print them 
  ## with grid.arrange.
  
  xs <- split(sumProbData,f = sumProbData$CohortAge)
  
  plotList <- list()
  
  plotList[[1]] <- ggplot(xs[[1]],aes(CohortStartTime/30,Mean,group=TreatmentName,color=TreatmentName)) +
    geom_point() +
    geom_line() +
    geom_ribbon(aes(ymin=Down,ymax=Up,fill=TreatmentName),alpha=0.5) +
    theme(legend.position="top") +
    labs(color=NULL,fill=NULL) +
    xlab("Study time (months)") +
    ylab(ylb) +
    scale_x_continuous(breaks=c(0,6,12,18)) +
    facet_grid(CohortAge~TreatmentName)
  
  for(i in 2:length(xs)) {
    plotList[[i]] <- plotList[[1]] %+% xs[[i]]
  }
  
  retVal <- plotList 
  
  ## print the plots  
  plotList$nrow <- length(xs)
  do.call(grid.arrange,plotList)
  
  # Return the plot list without printing
  return(invisible(retVal))
}
