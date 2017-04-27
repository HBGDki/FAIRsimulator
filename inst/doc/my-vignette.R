## ----message=FALSE, warning=FALSE----------------------------------------
library(FAIRsimulator)

## ------------------------------------------------------------------------

## Set the recrutment rate to something really fast
RecruitmentRatefunction<-function(StudyObj,Cohort) {
    return(5000) #Instantaneous randomization
}
StudyObj <- createStudy(RecruitmentFunction = RecruitmentRatefunction)

## ------------------------------------------------------------------------
StudyObj <- AdaptiveStudy(StudyObj)

## ----fig.width=6---------------------------------------------------------

plotStudyCohorts(StudyObj)


## ----fig.width=6---------------------------------------------------------
plotActiveSubjects(StudyObj)

## ----fig.width=6---------------------------------------------------------
plotHAZ(StudyObj)

## ----fig.width=8,fig.height=8--------------------------------------------
plotHAZTreatmentEff(StudyObj)

## ------------------------------------------------------------------------
kable(StudyObj$CohortList %listmap% "RandomizationProbabilities",digits = 2)

## ------------------------------------------------------------------------
kable(StudyObj$CohortList %listmap% "UpdateProbabilities",digits=2)

