library(FAIRsimulator)


StudyObj <- createStudy(latestTimeForNewBirthCohorts=18*30,studyStopTime = 32*30,
                        randomizationProbabilities = list(rep(0.25,5),rep(0.25,5),rep(0.25,5)),
                        minAllocationProbabilities = list(c(0.25,rep(0,4)),c(0.25,rep(0,4)),c(0.25,rep(0,4))),
                        treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
                        effSizes = list(c(0,0.05,0.07,0.1,0.25),c(0,0,0.05,0.1,0.25),c(0,0,0.05,0.25,0.3)))

StudyObj<-AdaptiveStudy(StudyObj)

plotStudyCohorts(StudyObj)
plotHAZTreatmentEff(StudyObj)
plotActiveSubjects(StudyObj)
plotHAZ(StudyObj)

