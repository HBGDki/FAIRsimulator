## This example simulates a parallel design with 5 treatments and 3 cohorts all with the same sampling design (monthly samples).
## Immediate recruitment and dropout of 20% for each cohort.
## No residual error.
## A final analysis of each cohort is performed at end of study.

library(FAIRsimulator)

set.seed(324124)

StudyObjIni <- createStudy(latestTimeForNewBirthCohorts=0,studyStopTime = 6*30+1,
                        nSubjects = c(1000,1000,1000),cohortStartTimes = c(0,0,0),
                        samplingDesign = list(c(0,1,2,3,4,5,6)*30,c(0,1,2,3,4,5,6)*30,c(0,1,2,3,4,5,6)*30),
                        randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
                        minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
                        treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
                        effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
                        Recruitmentfunction=function(...) {return(5000)},
                        newCohortLink = list(NULL,NULL,NULL),
                        dropoutRates = rep(0.2/(6*30),3))

#Set residual error to 0
StudyObjIni$InitEvent <- function(StudyObj) {
  StudyObj<-FAIRsimulator::InitEvent(StudyObj)
  StudyObj$sig<-0 #Set residual error to 0
  return(StudyObj)
}


StudyObj<-AdaptiveStudy(StudyObjIni)

## Plot the design
plotStudyCohorts(StudyObj)

## plot the number of active subjects per treatment cycle and cohort
plotActiveSubjects(StudyObj)

## Plot the HAZ profiles versus age.
plotHAZ(StudyObj)

## Plot the HAZ data and the treatment effects
plotHAZTreatmentEff(StudyObj)

# Plot the randomization probabilities at end of study
plotProbs(StudyObj,strProb = "UnWeightedUpdateProbabilities")
print("LMER Coefficients for each cohort:")
StudyObj$CohortList[[1]]$UpdateCoefficients
StudyObj$CohortList[[2]]$UpdateCoefficients
StudyObj$CohortList[[3]]$UpdateCoefficients

getProbData(StudyObj,strProb = "UnWeightedUpdateProbabilities")

## Running the study multiple times to get statistics
iter   <- 20
ncores <- 7

system.time(myMultStud <- runMultiSim(StudyOnjIni,iter=iter,ncores=ncores))

probDf           <- myMultStud[[2]]  # The Randomization probabilities
probDfUnweighted <- getMultiProbList(myMultStud[[1]],ncores=ncores,strProb="UnWeightedUpdateProbabilities")  # The UnweightedRandomization probabilities

plotMultiProb(probDf)
plotMultiProb(probDfUnweighted,ylb="Unweighted randomization probabilities")

## Summary of the end of study probabilities
sumData <- 
  probDfUnweighted %>% 
  group_by(CohortAge) %>% 
  filter(RandStudyTime == max(RandStudyTime)) %>% 
  group_by(TreatmentName,CohortAge) %>% 
  summarise(Mean=mean(Prob),Low=quantile(Prob,p=0.025),High=quantile(Prob,p=0.975))

xs <- split(sumData,f = sumData$CohortAge)
plotList <- list()

plotList[[1]] <- 
  ggplot(xs[[1]],aes(TreatmentName,Mean,color=TreatmentName)) +
  geom_point() +
  geom_errorbar(aes(ymin=Low,ymax=High),width=0.2) +
  facet_wrap(~CohortAge,scales = "free_x")

for(i in 2:length(xs)) {
  plotList[[i]] <- plotList[[1]] %+% xs[[i]]
}

## print the plots  
plotList$nrow <- length(xs)
do.call(grid.arrange,plotList)


