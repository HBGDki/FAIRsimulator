## This example simulates the same design multiple times in parallel.

library(FAIRsimulator)

set.seed(32423)


probTemperation <- function(probs) {
  probs <- sqrt(probs)/sum(sqrt(probs))
  return(probs)
}

StudyObjIni <- createStudy(latestTimeForNewBirthCohorts=18*30,studyStopTime = 32*30,
                        nSubjects = c(320,320,320),
                        randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
                        minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
                        treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
                        effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
                        Recruitmentfunction=function(...) {return(5000)},
                        probTemperationFunction = probTemperation)


iter   <- 7
ncores <- 7

system.time(myMultStud <- runMultiSim(StudyObjIni,iter=iter,ncores=ncores))

probDf           <- myMultStud[[2]]  # The Randomization probabilities
probDfUnweighted <- getMultiProbList(myMultStud[[1]],ncores=ncores,strProb="UnWeightedRandomizationProbabilities")  # The UnweightedRandomization probabilities
probDfUnweightedUpdate <- getMultiProbList(myMultStud[[1]],ncores=ncores,strProb="UpdateProbabilities")  # The UnweightedRandomization probabilities

plotMultiProb(probDf)
plotMultiProb(probDfUnweighted,ylb="Unweighted randomization probabilities")
plotMultiProb(probDfUnweightedUpdate,ylb="Unweighted update probabilities")

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
  
  
