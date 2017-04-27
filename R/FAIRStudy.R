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
  if (all(lapply(Cohort$SubjectList,function(Subject,LastSample){
    return(Subject$Status==0 || (LastSample %in% Subject$SubjectSampleTime))
  },max(Cohort$SamplingDesign))==TRUE)) {
    #browser()
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
#' @description Impute covariate - simple median imputation per individual (strIDVariable) of data frame (without correlation) of covariates in covariates vector.
#'
#' @param df A \code{data.frame} in which to impute missing covariate values.
#' @param StudyObj A FAIRsimulator \code{study} object
#' @param strID A string with the name of the columns in \code{df} containing the subject identifiers.
#' 
#' @return A \code{data.frame} with all missing covariates imputed.
#' @export
#'
#' @examples
#' \dontrun{}
ImputeCovariates <- function(df,StudyObj,strIDVariable = "ID") {

  for (i in 1:length(StudyObj$StudyDesignSettings$Covariates)) {
    strcov<-StudyObj$StudyDesignSettings$Covariates[i]
    df[is.na(df[[strcov]]),strcov]<-median(df[!duplicated(strIDVariable),strcov],na.rm=TRUE)
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
#' @param Cohort A FAIRsimulator \code{cohort} object
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A \code{data.frame}
#' @export
#'
#' @examples
#' \dontrun{}
GetCohortData<-function(Cohort,StudyObj) {
  DebugPrint(paste0("Assembly cohort ",Cohort$Name," data available at time ",StudyObj$CurrentTime),2,StudyObj)
  df<-data.frame()
  if (!is.null(Cohort$SubjectList)) {
    dflist<-lapply(Cohort$SubjectList,function(Subject){
      dfr<-data.frame()
      if (!is.null(Subject$Data)) {
        for (i in 1:length(Subject$Data)) {
          dfr<-rbind(dfr,data.frame(ID=Subject$StudyID,DATA=Subject$Data[[i]],AGE=Subject$SampleAge[[i]],TRT=Subject$TreatmentIndex,TRTS=Subject$Treatment,Subject$Covariates))
        }
      }
      return(dfr)
    })
    df<-as.data.frame(rbind(data.table::rbindlist(dflist)))
  }
  return (df)  
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
  df<-GetCohortData(Cohort,StudyObj) #Get Cohorts data up to this point in time
  df[df==-99]<-NA #Set -99 to missing
  df<-ImputeCovariates(df,StudyObj) #Impute missgin covariates
  
  ##### Make some covariate factors
  df$TRT<-as.factor(df$TRT)
  df$SEXN<-as.factor(df$SEXN)
  df$SANITATN<-as.factor(df$SANITATN)
  df$AGE<-df$AGE/(12*30) #Rescale time to years
  
  library(lme4)
  #### Perform LME estimation based on some covariates and treatment effects for each cohort
  lmefit <- lmer(paste0("DATA~1 + AGE + AGE:TRT + (1+AGE|ID) +",paste0(StudyObj$StudyDesignSettings$Covariates,collapse = " + ")),data=df,REML=FALSE)
  
  ##### Calculate new probabilites based on another cohort LME results
  lmecoef<-summary(lmefit)$coefficients[,1] #Get coefficicents from LME
  lmese<-summary(lmefit)$coefficients[,2] #Get SE from LME
  lmecoef<-lmecoef[regexpr('AGE:TRT.*',names(lmecoef))==1]
  lmese<-lmese[regexpr('AGE:TRT.*',names(lmese))==1]
  
  #DebugPrint(paste0("Estimated treatment effect in cohort ",Cohort$Name," at time ",StudyObj$CurrentTime),1,StudyObj)
  #DebugPrint(lmecoef,1,StudyObj)
  
  probs<-GetNewRandomizationProbabilities(trtcoeff=lmecoef,trtse=lmese,StudyDesignSettings$iNumPosteriorSamples) #Calculate randomization probs based on posterior distribution
  print(probs)
  
  #Update probabilitites to keep pre-defined portions 
  probs[Cohort$MinAllocationProbabilities!=0]<-Cohort$MinAllocationProbabilities[Cohort$MinAllocationProbabilities!=0]
  probspresum<-sum(probs[Cohort$MinAllocationProbabilities!=0])
  probssum<-sum(probs[Cohort$MinAllocationProbabilities==0])
  probs[Cohort$MinAllocationProbabilities==0]<-(1-probspresum)*probs[Cohort$MinAllocationProbabilities==0]/probssum
  DebugPrint(paste0("Recalculated randomization probabilities in ",Cohort$Name," at time ",StudyObj$CurrentTime),1,StudyObj)
  print(probs)
  
  Cohort$UpdateProbabilities<-probs #The latest probability updates
  StudyObj$CohortList[[cohortindex]]<-Cohort #Save the updated cohort
  
  for (j in 1:length(StudyObj$CohortList)) {#Update all dependent cohorts
    if (!is.null(StudyObj$CohortList[[j]]$ProbabilityCohort) && cohortindex!=j && cohortindex==StudyObj$CohortList[[j]]$ProbabilityCohort) {#If cohort j should be updated based on prob in cohort i
      DebugPrint(paste0("Updating probabilities in cohort ",StudyObj$CohortList[[j]]$Name," based on analysis in cohort ",Cohort$Name," at time ",StudyObj$CurrentTime),1,StudyObj)
      RandProbs<-list() #Save previous randomization probabilities on cohort
      RandProbs$CohortTime<-StudyObj$CohortList[[j]]$CurrentTime #The time until the probability was valid
      RandProbs$StudyTime<-StudyObj$CurrentTime
      RandProbs$FromCohort<-cohortindex
      RandProbs$RandomizationProbabilities<-StudyObj$CohortList[[j]]$RandomizationProbabilities
      StudyObj$CohortList[[j]]$PreviousRandomizationProbabilities[[length(StudyObj$CohortList[[j]]$PreviousRandomizationProbabilities)+1]]<-RandProbs
      StudyObj$CohortList[[j]]$RandomizationProbabilities<-probs #Update the probabilities for the child cohort
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
        cohortInterimAnalysis<-InterimAnalyzesTime(Cohort,StudyObj)
        if (cohortCompleted || cohortInterimAnalysis) { #If we should update probabilities of child cohorts
          StudyObj<-UpdateProbabilities(Cohort,StudyObj,i)
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
  Crit<-StudyObj$CurrentTime>30*29
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
  Cohort$RecruitmentRateFunction<-RecruitmentRatefunction # The recruitment rate function
  if (!is.null(CohortNum)) {
    Cohort$MaxNumberOfSubjects<-StudyObj$StudyDesignSettings$MaxNumberofSubjects[[CohortNum]] #The maximum number of subject in this cohort
    Cohort$SamplingDesign<-StudyObj$StudyDesignSettings$SamplingDesigns[[CohortNum]] #The sampling design for this cohort
    Cohort$RandomizationProbabilities<-StudyObj$StudyDesignSettings$RandomizationProbabilities[[CohortNum]] #Get the randomization probabilities
    Cohort$UpdateProbabilities<-Cohort$RandomizationProbabilities #The Current updated probabilities
    Cohort$MinAllocationProbabilities<-StudyObj$StudyDesignSettings$MinAllocationProbabilities[[CohortNum]] #Get the minimum randomization probabilities
    Cohort$Treatments<-StudyObj$StudyDesignSettings$Treatments[[CohortNum]] #Get the treatments (treatment codes)
    Cohort$EffSizes<-StudyObj$StudyDesignSettings$EffSizes[[CohortNum]] #Get the effect sizes for each treatment
    Cohort$RandomizationAgeRange<-StudyObj$StudyDesignSettings$CohortAgeRange[[CohortNum]] #The age range at randomization for this cohort
    Cohort$DropoutRate<-StudyObj$StudyDesignSettings$CohortDropoutRate[CohortNum] #The dropout rate for this cohort
    Cohort$NewCohortLink<-StudyDesignSettings$NewCohortLink[[CohortNum]] #The link to another cohort to get the randomization probabilities
    Cohort$Name<-paste0("C-",CohortNum,"-1 [",Cohort$RandomizationAgeRange[1]/30,"-",Cohort$RandomizationAgeRange[2]/30,"m @ rand]")
    Cohort$StartNum<-CohortNum #The cohort start number
    Cohort$CycleNum<-1 #The cycle number
  }
  class(Cohort)<-"cohort"
  return(Cohort)  
}