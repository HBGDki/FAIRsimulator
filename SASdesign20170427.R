library(FAIRsimulator)
library(gridExtra)

StudyObj <- createStudy(latestTimeForNewBirthCohorts=18*30,studyStopTime = 32*30,
                        randomizationProbabilities = list(rep(0.20,5),rep(0.20,5),rep(0.20,5)),
                        #minAllocationProbabilities = list(c(0,rep(0,4)),c(0,rep(0,4)),c(0,rep(0,4))),
                        minAllocationProbabilities = list(c(0.2,rep(0,4)),c(0.2,rep(0,4)),c(0.2,rep(0,4))),
                        treatments =list(c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4"),c("SoC-2","TRT-5","TRT-6","TRT-7","TRT-8"),c("SoC-3","TRT-9","TRT-10","TRT-11","TRT-12")),
                        effSizes = list(c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25),c(0,0.05,0.1,0.15,0.25)),
                        Recruitmentfunction=function(...) {return(5000)})

StudyObj<-AdaptiveStudy(StudyObj)

plotStudyCohorts(StudyObj)
plotHAZTreatmentEff(StudyObj)
plotActiveSubjects(StudyObj)
plotHAZ(StudyObj)

plotProbs(StudyObj)
plotProbs(StudyObj,strProb="UnWeightedRandomizationProbabilities")

plotProbs <- function(StudyObj,strProb="RandomizationProbabilities") {
  
  allData <- getAllSubjectData(StudyObj)
  
  AgeRanges<-StudyObj$StudyDesignSettings$CohortAgeRange
  
  myCohorts <- cohorts(StudyObj)
  
  resDf <- data.frame()
  for (i in 1:length(myCohorts))  {
    Cohort <- myCohorts[[i]]
    for (j in 1:length(AgeRanges)) {
      if (AgeRanges[[j]][1] == Cohort$RandomizationAgeRange[1] && AgeRanges[[j]][2] == Cohort$RandomizationAgeRange[2]) {
        resDf[i,"CohortName"] <- Cohort$Name
        resDf[i,"AgeRange"] <- paste0(AgeRanges[[j+Cohort$CycleNum-1]][1]/30, "-",AgeRanges[[j+Cohort$CycleNum-1]][2]/30," month")
      }
    }
  }
  
  allData <- left_join(allData,resDf,by="CohortName") %>% mutate(CohortAge = ifelse(AgeRange=="0-1 month","0-6 months",ifelse(AgeRange=="6-7 month","6-12 months","12-18 months")))
  
  cohrtNams    <- dimnames(StudyObj$CohortList %listmap% strProb)[[2]]
  probs        <- data.frame(StudyObj$CohortList %listmap% strProb) 
  names(probs) <- cohrtNams
  
  treatnames        <- data.frame(StudyObj$CohortList %listmap% "Treatments")
  names(treatnames) <- cohrtNams
  
  
  probs      <- probs      %>%   gather(key="CohortName",value="Prob") 
  treatnames <- treatnames %>%   gather(key="CohortName",value="TreatmentName") 
  
  probs$TreatmentName <- treatnames$TreatmentName
  
  randProbs <- allData %>% distinct(CohortAge,CohortName,.keep_all=TRUE) %>% dplyr::select(CohortAge,CohortName,RandStudyTime) %>% left_join(probs,.,by="CohortName")
  randProbs$TreatmentName <- factor(randProbs$TreatmentName,levels=c("SoC-1","TRT-1","TRT-2","TRT-3","TRT-4","SoC-2","TRT-5","TRT-6","TRT-7","TRT-8","SoC-3","TRT-9","TRT-10","TRT-11","TRT-12"))

  xs <- split(randProbs,f = randProbs$CohortAge)
  
  p1 <- ggplot(xs$'0-6 months',aes(RandStudyTime,Prob,group=TreatmentName,color=TreatmentName)) +
    geom_point() +
    geom_line() +
    theme(legend.position="top") +
    labs(color=NULL) +
    facet_wrap(~CohortAge)
  
  p2 <- p1 %+% xs$'6-12 months'
  p3 <- p1 %+% xs$'12-18 months'
  
  grid.arrange(p1,p2,p3,ncol=3)
}