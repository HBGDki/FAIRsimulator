## ----message=FALSE, warning=FALSE----------------------------------------
library(FAIRsimulator)

## ------------------------------------------------------------------------
# Run the simulations in the exampel file
#source("../FAIRstudy.R",echo=FALSE)
load("../studyobj.Rdata")

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

