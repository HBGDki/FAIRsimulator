
library(FAIRsimulator) ## Access functionality to run the simulations. See help(package="FAIRsimulator").

#source('./AdaptiveStudy.R') # Now part of the FAIRsimulator package

# source('./../../GAT-Growth-PMX-1/Supermodel/Analysis/SupportFunctions/calcHAZVector.R')  # Now part of the FAIRsimulator package



#The Study increment event




















#Get the number of recruited subjects at current time
GetRecruitedSubjects<-function(StudyObj,Cohort) { 
  pcount <- rpois(n=1, lambda=Cohort$RecruitmentRateFunction(StudyObj,Cohort)) #Recruit according to a Poisson distribution
  pcount<-min(pcount,Cohort$MaxNumberOfSubject-Cohort$NumberOfRecruitedSubjects) #Make sure we're not recruting too many subjects
  return(pcount)
}

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

#Get the treatment randomizations
GetTreatmentRandomizations<-function(Num,Cohort,StudyObj) { #Simple uniform random allocation between treatments based on randomization probabilities
  DebugPrint(paste0("Allocated ",Num," new subjects to treatments"),3,StudyObj)
  rnd<-runif(Num)
  treatmentIndex<-rep(NA,Num)
  for (i in 1:Num){
    lowprob<-0
    highprob<-0
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

#Returns the uniform age at randomization given an age range
GetAgesAtRandomization<-function(NumRecruited,Cohort,StudyObj) {
  DebugPrint(paste0("Simulate age at randomization for ",NumRecruited," subjects uniformly within [",Cohort$RandomizationAgeRange[1]/30,"-",Cohort$RandomizationAgeRange[2]/30,"] month"),3,StudyObj)
  return(runif(n=NumRecruited,min = Cohort$RandomizationAgeRange[1],max=Cohort$RandomizationAgeRange[2]))
}

#Constant hazard dropout function
DropoutTime<-function(lambda,randnum=NULL){
  if (is.null(randnum)) randnum<-runif(1)
  return(-log(randnum)/lambda)
}

GetSubject<-function(ID,TRTIndex,AgeAtRand,Cohort,StudyObj) {
  library(MASS)
  
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

#
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

RecruitmentEvent<-function(StudyObj) { #Event that recruites subjects to cohort
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

AddEffect<-function(y,StudyObj,Subject,time,effecttime=6*30,cumeffect=NULL) {
  DebugPrint(paste0("Adding HAZ effect for subject ",Subject$StudyID," at sample time = ",time,""),4,StudyObj)
  cumeff<-0
  #If we have a cumulative effect from previous treatments, i.e. at time=0
  if (!is.null(Subject$CumulativeEffect)) cumeff<-Subject$CumulativeEffect
  if (!is.null(cumeffect)) cumeff<-cumeffect #Use this as previous effect instead
return(y+Subject$TreatmentEff/effecttime*time + cumeff) #Assume linear effect (and additive)
}

SimulateHAZ<-function(StudyObj,Subject,time,age) {
  DebugPrint(paste0("Simulate HAZ for subject ",Subject$StudyID," at study time: ",StudyObj$CurrentTime," (sample time = ",time,")"),4,StudyObj)
  library(MASS)
  res_err<-mvrnorm(n = length(age), 0, StudyObj$sig , tol = 1e-6, empirical = FALSE, EISPACK = FALSE) #Simulate residual error on HAZ
  y = StudyObj$calcHAZVector(age = age/(30*12),basethetas = StudyObj$thbasevector,covthetas = Subject$FREMCoeff,etas = Subject$IndSamples) #Calculate Y without residual error
  yres<-y+res_err #Add residual error
  yres<-AddEffect(yres,StudyObj,Subject,time) #Add HAZ Effect
  DebugPrint(paste0("HAZ observation ",yres," simulated for subject with age ",age),4,StudyObj)
  return(yres)
}

SimulateDataEvent<-function(StudyObj) { #Simulate data for all subjects which should have a sample at this date
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

DropoutEvent<-function(StudyObj) { #Event that check if subject have dropped out
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

MoveSubjects<-function(FromCohort,ToCohort,StudyObj) { #Move subjects from FromCohort to ToCohort which are completed (Status==2), Re-randomize treatments based on ToCohort rand probabilities
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

#Return the cohortlist index of a certain start number
GetStartNumIndex <- function(StudyObj,StartNum) {
for (i in 1:length(StudyObj$CohortList)) {
  if (StartNum==StudyObj$CohortList[[i]]$StartNum) return(i)
}
return (NA)
}

#Event that check if completed subjects should be moved to new Cohort
#Only subjects that are in a cohort that is linked to another cohort is affected
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
          NewChildCohort$UpdateProbabilities<-NewChildCohort$RandomizationProbabilities
          NewChildCohort$MinAllocationProbabilities<-StudyObj$CohortList[[NewCohortLinkIndex]]$MinAllocationProbabilities #Get the minimum randomization probabilities
          NewChildCohort$Treatments<-StudyObj$CohortList[[NewCohortLinkIndex]]$Treatments #Get the treatments (treatment codes)
          NewChildCohort$EffSizes<-StudyObj$CohortList[[NewCohortLinkIndex]]$EffSizes #Get the effect sizes for each treatment
          NewChildCohort$RandomizationAgeRange<-Cohort$RandomizationAgeRange #The age range at randomization for this cohort
          NewChildCohort$DropoutRate<-StudyObj$CohortList[[NewCohortLinkIndex]]$DropoutRate #The dropout rate for this cohort
          NewChildCohort$NewCohortLink<-StudyObj$CohortList[[NewCohortLinkIndex]]$NewCohortLink #The link to another cohort to get the randomization probabilities

          NewChildCohort$Name<-paste0("C-",Cohort$StartNum,"-",Cohort$CycleNum+1," [",Cohort$RandomizationAgeRange[1]/30,"-",Cohort$RandomizationAgeRange[2]/30,"m @ rand]")
          NewChildCohort$StartNum<-Cohort$StartNum #The Cohort starting number
          NewChildCohort$CycleNum<-Cohort$CycleNum+1 #The Cohort cycle number
          StudyObj$CohortList[[i]]$ChildCohort<-length(StudyObj$CohortList)+length(NewChildCohortList)+1 #Add reference to the new cohort as a child cohort
          NewChildCohort$ParentCohort<-i #Add a reference to old Cohort as a parent cohort
          NewChildCohort$ProbabilityCohort<-NewCohortLinkIndex #The cohort where to update the probabilities from
          
          NewChildCohort<-MoveSubjects(Cohort,NewChildCohort,StudyObj)          
          
          NewChildCohortList[[length(NewChildCohortList)+1]]<-NewChildCohort
        } else {
          StudyObj$CohortList[[Cohort$ChildCohort]]<-MoveSubjects(Cohort,StudyObj$CohortList[[Cohort$ChildCohort]],StudyObj)          
        }       
      }
    }
    
    ### Add NewChildCohorts to Cohort list
    if (length(NewChildCohortList)!=0) {
      for (i in 1:length(NewChildCohortList)) {
        StudyObj$CohortList[[(length(StudyObj$CohortList)+1)]]<-NewChildCohortList[[i]]
      }
    }
  }
  return(StudyObj)
}


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


