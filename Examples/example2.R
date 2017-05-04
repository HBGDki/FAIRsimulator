## This example simulates the same design multiple times in parallel.

library(FAIRsimulator)
library(foreach)
library(doParallel)
library(purrr)

set.seed(32423)


probTemperation <- function(probs) {
  probs <- sqrt(probs)/sum(sqrt(probs))
  return(probs)
}

StudyObjIni <- createStudy(latestTimeForNewBirthCohorts=18*30,studyStopTime = 32*30,
                        nSubjects = c(320,320,320),
                        randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
                        #minAllocationProbabilities = list(c(0,rep(0,4)),c(0,rep(0,4)),c(0,rep(0,4))),
                        minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
                        treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
                        effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
                        Recruitmentfunction=function(...) {return(5000)},
                        minSubjects = 10,
                        probTemperationFunction = probTemperation)


iter   <- 28
ncores <- 7

runMultiSim <- function(StudyOnjIni,extractProbs=TRUE,iter=1,ncores=1,...) {
  
  ## Start the parallell engine if ncores > 1
  if(ncores>1) {registerDoParallel(cores=ncores)}
  
  ## Run the simulations
  myRes <- foreach(i=1:iter) %dopar% AdaptiveStudy(StudyObjIni)
  
  ## Extract the probabilities if requested
  if(extractProbs) {
    probDf  <- getMultiProbList(myRes,ncores=ncores,...)
    retList <- list(studList = myRes,probDf = probDf)
  } else {
    retList <- list(studList = myRes)
  }
  
  return(retList)
  
}

getMultiProbList <- function(multiStudObj,ncores=1,...) {

  ## Start the parallell engine if ncores > 1
  if(ncores>1) {registerDoParallel(cores=ncores)}
  
  ## get a list of all probability data
  probList <- foreach(i=1:length(multiStudObj)) %dopar% getProbData(multiStudObj[[i]],...)
  
  ## Create one data frame of all the probability data frames
  probDf <- bind_rows(probList)
  
  ## Add an iteration variable (Iter)
  probDf$Iter <- rep(1:length(multiStudObj),each=nrow(probDf)/length(multiStudObj))

  return(probDf)  
}

myMultStud <- runMultiSim(StudyOnjIni,iter=iter,ncores=ncores)

probDf           <- myMultStud[[2]]  # The Randomization probabilities
probDfUnweighted <- getMultiProbList(myMultStud[[1]],ncores=ncores,strProb="UnWeightedRandomizationProbabilities")  # The UnweightedRandomization probabilities

plotMultiProb(probDf)

plotMultiProb(probDfUnweighted)

plotMultiProb <- function(probDf,ylb="Randomization probability",pup = 0.95,pdo=0.05) {

  ## Define functions to compute the upper and lower limits of the prediction intervals of the probabilities
  up95 <- function(x,PRup=pup) quantile(x,p=PRup)
  do05 <- function(x,PRdo=pdo) quantile(x,p=PRdo)

  ## Compute the plot data
  sumProbData <- probDf %>% 
    group_by(CohortAge,TreatmentName,RandStudyTime) %>% 
    summarise(Mean=mean(Prob),Up=up95(Prob),Down=do05(Prob)) 

  ## To get legends for each age cohort we need to create separate graphs for the. We'll split the data and create the graphs and then print them 
  ## with grid.arrange.
  
  xs <- split(sumProbData,f = sumProbData$CohortAge)

  plotList <- list()

  plotList[[1]] <- ggplot(xs[[1]],aes(RandStudyTime/30,Mean,group=TreatmentName,color=TreatmentName)) +
    geom_point() +
    geom_line() +
    geom_ribbon(aes(ymin=Down,ymax=Up,fill=TreatmentName),alpha=0.5) +
    theme(legend.position="top") +
    labs(color=NULL,fill=NULL) +
    xlab("Study time (months)") +
    ylab(ylb) +
    facet_grid(CohortAge~TreatmentName)
  
  for(i in 2:length(xs)) {
    plotList[[i]] <- plotList[[1]] %+% xs[[i]]
  }
  
  retVal <- plotList 
  
  ## print the plots  
  plotList$nrow <- length(xs)
  do.call(grid.arrange,plotList)
  
  # Return the plot list without printing
  return(invisible(retVal))
}
