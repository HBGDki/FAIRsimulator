library(FAIRsimulator)

RecruitmentRatefunction<-function(StudyObj,Cohort) {
   return(5000) #Instantaneous randomization
 # return(20) #20 randomized subjects per time unit
}

StudyObj <- createStudy(RecruitmentFunction = RecruitmentRatefunction)
StudyObj<-AdaptiveStudy(StudyObj)

plotStudyCohorts(StudyObj)
plotActiveSubjects(StudyObj)
plotHAZ(StudyObj)

