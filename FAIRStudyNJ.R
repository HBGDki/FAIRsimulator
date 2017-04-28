library(FAIRsimulator)

# RecruitmentRatefunction<-function(StudyObj,Cohort) {
#    return(5000) #Instantaneous randomization
#  # return(20) #20 randomized subjects per time unit
# }



# StudyObj <- createStudy(nCohorts = 4,cohortStartTimes = c(2,1,0,6*30+20),samplingDesign = list(c(0,1,2,3,4,5,6)*30,c(0/30,3,6)*30,c(0/30,3,6)*30,c(0,1,2,3,4,5,6)*30),
#                         nSubjects = c(320,320,320,320),dropoutRates = c(0.2/(6*30),0.2/(6*30),0.2/(6*30),0.2/(6*30)), recruitmentAges = list(c(0,0.5)*30,c(6,6.5)*30,c(12,12.5)*30,c(0,0.5)*30))

createStudy <- function(nCohorts = 3, cohortStartTimes = c(2,1,0), samplingDesign = list(c(0,1,2,3,4,5,6)*30,c(0/30,3,6)*30,c(0/30,3,6)*30),nSubjects = c(320,320,320), 
                        dropoutRates = c(0.2/(6*30),0.2/(6*30),0.2/(6*30)),randomizationProbabilities = list(c(0.25,0.25,0.25,0.25),c(0.25,0.25,0.25,0.25),c(0.25,0.25,0.25,0.25)),
                        minAllocationProbabilities = list(c(0.25,0,0,0),c(0.25,0,0,0),c(0.25,0,0,0)),
                        treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3"),c("SoC-2","TRT-4","TRT-5","TRT-6"),c("SoC-3","TRT-7","TRT-8","TRT-9")),
                        effSizes = list(c(0,0.05,0.1,0.25),c(0,0,0.05,0.25),c(0,0.05,0.25,0.3)),
                        newCohortLink = list(2,3,NULL),
                        recruitmentAges = list(c(0,1)*30,c(6,7)*30,c(12,13)*30),RecruitmentFunction=RecruitmentRatefunction,
                        debugLevel=1,studyStopTime=30*28,latestTimeForNewBirthCohorts=0,
                        strCovariates=c("BIRTHWT","MAGE","MHTCM","SEXN","SANITATN"),
                        currentDate = Sys.Date()) {
  
  ## Create the study design setting list ##
  
  StudyDesignSettings                            <- list()
  StudyDesignSettings$CohortNumbers              <- nCohorts
  StudyDesignSettings$MaxNumberofSubjects        <- nSubjects
  
  StudyDesignSettings$CohortStartTimes           <- cohortStartTimes
  
  StudyDesignSettings$RandomizationProbabilities <- randomizationProbabilities #Initial randomization probabilities
  StudyDesignSettings$MinAllocationProbabilities <- minAllocationProbabilities #Minimum allocation probabilities for each treatment
                                                       

  StudyDesignSettings$iNumPosteriorSamples       <- 10000 #The number of samples to calculate prob of beeing best
  
  StudyDesignSettings$Treatments <- treatments #The treatment
  
  StudyDesignSettings$EffSizes <- effSizes #The actual effect sizes
  
  #The age ranges for each cohort
  StudyDesignSettings$CohortAgeRange <- recruitmentAges
  
  #The sampling design for each pre-defined cohort
  StudyDesignSettings$SamplingDesigns <- samplingDesign
  
  #When a cohort is evolving to a new cohort, the information about randomization should be based on the cohort in the list, NULL= no evolving
  StudyDesignSettings$NewCohortLink <- newCohortLink 
  
  #If the last sample will be the baseline sample (subject cohorttime = 0) in the new cohort 
  StudyDesignSettings$MoveLastSampleToNewCohort <- TRUE 
  
  #The dropout rates
  StudyDesignSettings$CohortDropoutRate <- dropoutRates
  
  #The stop time of the study
  StudyDesignSettings$StudyStopTime <- studyStopTime
  
  #The latest time to start new birth cohorts, set to 0 if no new birth cohorts should be added
  StudyDesignSettings$LatestTimeForNewBirthCohorts<-latestTimeForNewBirthCohorts
  
  #The covariates that should be stored on each subject, e.g. to be used in the lme analysis
  StudyDesignSettings$Covariates <- strCovariates  
  
  #The current ID number, i.e. global counter of ID number
  StudyDesignSettings$CurrentID  <- 1 
  
  ## Create the study object ##
  StudyObj <- list()
  
  StudyObj$CurrentTime <- 0 #Let the study be time driven
  StudyObj$CurrentDate <- currentDate #The date that the study will start
  
  ## Add recruitmentrate function
  StudyObj$Recruitmentfunction <- RecruitmentRatefunction
  
  ## Define events
  StudyObj$InitEvent           <- InitEvent
  StudyObj$StudyIncrementEvent <- StudyIncrementEvent
  StudyObj$StopEvent           <- StopEvent
  
  ## Set the debug level
  StudyObj$DebugLevel <- debugLevel # 4=Extreme output, #3=Print everything important, 2=Print events, 1=Sparse print, 0=Print nothing
  
  ## Create the list of Events
  EventList <- list()
  EventList[[length(EventList)+1]] <- AddCohortEvent #Add the AddCohort event
  EventList[[length(EventList)+1]] <- DropoutEvent #Dropout event
  EventList[[length(EventList)+1]] <- RecruitmentEvent #Add the Recruitment event
  EventList[[length(EventList)+1]] <- SimulateDataEvent #Add a Simulate data event
  EventList[[length(EventList)+1]] <- MoveCompletedSubjects #Move completed subjects event
  EventList[[length(EventList)+1]] <- AnalyzeDataEvent #Add a Analyze data event
  EventList[[length(EventList)+1]] <- UpdateProbabilitiesEvent #Update all probabilities, even for non directly connected cohorts
  EventList[[length(EventList)+1]] <- AddNewBirthCohortEvent #Add a new birth cohort
  
  StudyObj$EventList           <- EventList
  StudyObj$StudyDesignSettings <- StudyDesignSettings #Add the specific study design settings into the Study object
  
  class(StudyObj) <- "study"
  
  return(StudyObj)
}

StudyObj <- createStudy(RecruitmentFunction = RecruitmentRatefunction)
StudyObj<-AdaptiveStudy(StudyObj)

plotStudyCohorts(StudyObj)
plotActiveSubjects(StudyObj)
plotHAZ(StudyObj)

