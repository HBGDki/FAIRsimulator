
### Create the study object with all the design settign to send in to the Adaptive study

StudyDesignSettings<-list()
StudyDesignSettings$CohortNumbers<-3
StudyDesignSettings$MaxNumberofSubjects<-c(320,320,320)
#StudyDesignSettings$CohortStartTimes<-c(25,30,120)
StudyDesignSettings$CohortStartTimes<-c(4*30,2*30,0)

StudyDesignSettings$RandomizationProbabilities<-list(c(0.25,0.25,0.25,0.25), #Randomization probabilities for all possible treatments all cohorts
                                                     c(0.25,0.25,0.25,0.25),
                                                     c(0.25,0.25,0.25,0.25))
StudyDesignSettings$MinAllocationProbabilities<-list(c(0.25,0,0,0), #Minimum allocation probabilities for each treatment
                                                     c(0.25,0,0,0),
                                                     c(0.25,0,0,0))
StudyDesignSettings$iNumPosteriorSamples<-10000 #The number of samples to calculate prob of beeing best

StudyDesignSettings$Treatments<-list(c("SoC-1","TRT-1","TRT-2","TRT-3"),c("SoC-2","TRT-4","TRT-5","TRT-6"),c("SoC-3","TRT-7","TRT-8","TRT-9")) #Treatment codes
StudyDesignSettings$EffSizes<-list(c(0,0.05,0.1,0.25),c(0,0,0.05,0.25),c(0,0.05,0.25,0.3)) #EffectSizes for HAZ at 6 month of each treatment
StudyDesignSettings$CohortAgeRange<-list(c(0,1)*30,c(6,7)*30,c(12,13)*30) #The age ranges for each cohort

StudyDesignSettings$SamplingDesigns<-list(c(0,1,2,3,4,5,6)*30,c(0/30,3,6)*30,c(0/30,3,6,9,12)*30) #The sampling design for each pre-defined cohort

StudyDesignSettings$NewCohortLink<-list(2,3,NULL) #When a cohort is evolving to a new cohort, the information about randomization should be based on the cohort in the list, NULL= no evolving
StudyDesignSettings$MoveLastSampleToNewCohort<-TRUE #If the last sample will be the baseline sample (subject cohorttime = 0) in the new cohort 

StudyDesignSettings$CohortDropoutRate<-c(0.2/(6*30),0.2/(6*30),0.2/(6*30))
StudyDesignSettings$Covariates<-c("BIRTHWT","MAGE","MHTCM","SEXN","SANITATN")  #The covariates that should be stored on each subject
StudyDesignSettings$CurrentID<-1 #The current ID number, i.e. global counter of ID number



StudyObj<-list()
StudyObj$CurrentTime<-0 #Let the study be time driven
StudyObj$CurrentDate<-Sys.Date()

StudyObj$InitEvent<-InitEvent
StudyObj$StudyIncrementEvent<-StudyIncrementEvent
StudyObj$StopEvent<-StopEvent

StudyObj$DebugLevel<-1 #4=Extreme output, #3=Print everything important, 2=Print events, 1=Sparse print, 0=Print nothing

EventList<-list()
EventList[[length(EventList)+1]]<-AddCohortEvent #Add the AddCohort event
###Should subject move to another cohort?
  ### If yes, move subject

EventList[[length(EventList)+1]]<-DropoutEvent #Dropout event
EventList[[length(EventList)+1]]<-RecruitmentEvent #Add the Recruitment event
EventList[[length(EventList)+1]]<-SimulateDataEvent #Add a Simulate data event

EventList[[length(EventList)+1]]<-MoveCompletedSubjects #Move completed subjects event

EventList[[length(EventList)+1]]<-AnalyzeDataEvent #Add a Analyze data event

#EventList[[length(EventList)+1]]<-UpdateProbsEvent #UpdateAllProbabilities


StudyObj$EventList<-EventList

StudyObj$StudyDesignSettings<-StudyDesignSettings #Add the specific study design settings into the Study object

class(StudyObj) <- "study"


###### Call the Adaptive Study

StudyObj<-AdaptiveStudy(StudyObj)

print(paste0("The study stopped at time: ",StudyObj$CurrentTime, " i.e. ",StudyObj$CurrentDate))

## Do some plotting of the results
treatments <- sort(unique(unlist(StudyObj$StudyDesignSettings$Treatments)))


