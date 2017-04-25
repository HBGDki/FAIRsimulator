#library(tidyverse)
library(FAIRsimulator)
library(zoo)



tmp <- getSubjectItems(StudyObj[["CohortList"]][[1]][["SubjectList"]][[1]])
tmp <- getItemsFromSubjects(StudyObj[["CohortList"]][[1]])
tmp <- getAllSubjectData(StudyObj,scalarItems="AgeAtRand",covariates=NULL)


extractInclusionData <- function(StudyObj) {
  cohortDetails <- getCohortDetails(StudyObj)
  allData       <- getAllSubjectData(StudyObj,scalarItems=c("RandCohortTime","DropoutCohortTime"),covariates=NULL) 
  
  allRes <- NULL
  for(i in 1:nrow(cohortDetails)) {
    
    cycleNum <- cohortDetails[i,"CycleNum"]
    cohortID <- cohortDetails[i,"CohortID"]
    cohortStartTime <- cohortDetails[i,"CohortStartTime"]
    
    cohort1 <- allData %>% 
      filter(CohortID==cohortID,CycleNum==cycleNum) 
    
    recruited <- 
      cohort1 %>% 
      group_by(RandCohortTime) %>% 
      distinct(StudyID) %>%  
      tally %>% 
      ungroup %>% 
      mutate(Recruited = cumsum(n)) %>% 
      dplyr::select(-n) %>% 
      rename(Day=RandCohortTime)
    
    dropout <- cohort1 %>% 
      filter(!is.na(DropoutCohortTime)) %>% 
      group_by(DropoutCohortTime) %>% 
      distinct(StudyID) %>%  
      tally %>% 
      ungroup %>% 
      mutate(Dropout = cumsum(n)) %>% 
      dplyr::select(-n) %>% 
      rename(Day=DropoutCohortTime)
    
    completed <- 
      cohort1 %>% 
      group_by(StudyID) %>% 
      filter(CohortSampleTime == max(CohortSampleTime), is.na(DropoutCohortTime)) %>% 
      ungroup %>% 
      group_by(CohortSampleTime) %>% 
      tally %>% 
      ungroup %>% 
      mutate(Completed = cumsum(n)) %>% 
      dplyr::select(-n) %>% 
      rename(Day = CohortSampleTime)
    
    
    resDf <- data.frame(Day=0:max(cohort1$CohortSampleTime))
    resDf <- left_join(resDf,recruited,by="Day") %>% na.locf %>% mutate(Recruited = ifelse(is.na(Recruited),0,Recruited))
    resDf <- left_join(resDf,dropout,by="Day")   %>% na.locf %>% mutate(Dropout = ifelse(is.na(Dropout),0,Dropout))
    resDf <- left_join(resDf,completed,by="Day") %>% na.locf %>% mutate(Completed = ifelse(is.na(Completed),0,Completed))
    
    resDf <- resDf %>% mutate(Remaining = Recruited-Dropout-Completed,StudyDay = Day + cohortStartTime,CohortID=cohortID,CycleNum=cycleNum,CohortStartTime=cohortStartTime)
    
    if(is.null(allRes)) {
      allRes <- resDf
    } else {
      allRes <- bind_rows(allRes,resDf)
    }
    
  }
  
  return(allRes)
  
}

extractInclusionData2 <- function(StudyObj) {
  cohortDetails <- getCohortDetails(StudyObj)
  allData       <- getAllSubjectData(StudyObj,scalarItems=c("RandStudyTime","DropoutStudyTime"),covariates=NULL) 
  maxStudyTime  <- max(allData$StudySampleTime)
  
  allRes <- NULL
  for(i in 1:nrow(cohortDetails)) {
    cycleNum <- cohortDetails[i,"CycleNum"]
    cohortID <- cohortDetails[i,"CohortID"]
    cohortStartTime <- cohortDetails[i,"CohortStartTime"]
    
    cohort1 <- allData %>% 
      filter(CohortID==cohortID,CycleNum==cycleNum) 
    
    recruited <- 
      cohort1 %>% 
      group_by(RandStudyTime) %>% 
      distinct(StudyID) %>%  
      tally %>% 
      ungroup %>% 
      mutate(Recruited = cumsum(n)) %>% 
      dplyr::select(-n) %>% 
      rename(Day=RandStudyTime)
    
    dropout <- cohort1 %>% 
      filter(!is.na(DropoutStudyTime)) %>% 
      group_by(DropoutStudyTime) %>% 
      distinct(StudyID) %>%  
      tally %>% 
      ungroup %>% 
      mutate(Dropout = cumsum(n)) %>% 
      dplyr::select(-n) %>% 
      rename(Day=DropoutStudyTime)
    
    completed <- 
      cohort1 %>% 
      group_by(StudyID) %>% 
      filter(StudySampleTime == max(StudySampleTime), is.na(DropoutStudyTime)) %>% 
      ungroup %>% 
      group_by(StudySampleTime) %>% 
      tally %>% 
      ungroup %>% 
      mutate(Completed = cumsum(n)) %>% 
      dplyr::select(-n) %>% 
      rename(Day = StudySampleTime)
    
    
    resDf <- data.frame(Day=0:maxStudyTime)
    resDf <- left_join(resDf,recruited,by="Day") %>% na.locf %>% mutate(Recruited = ifelse(is.na(Recruited),0,Recruited))
    resDf <- left_join(resDf,dropout,by="Day")   %>% na.locf %>% mutate(Dropout = ifelse(is.na(Dropout),0,Dropout))
    resDf <- left_join(resDf,completed,by="Day") %>% na.locf %>% mutate(Completed = ifelse(is.na(Completed),0,Completed))
    
    resDf <- resDf %>% 
      mutate(Remaining = Recruited-Dropout-Completed,
             CohortDay = Day - cohortStartTime,CohortID=cohortID,
             CycleNum=cycleNum,CohortStartTime=cohortStartTime) %>% 
      rename(StudyDay = Day)
    
    if(is.null(allRes)) {
      allRes <- resDf
    } else {
      allRes <- bind_rows(allRes,resDf)
    }
    
  }
  
  return(allRes)
  
}

indInfo <- extractInclusionData2(StudyObj)


indInfoLong <- indInfo %>% gather(key,value,-StudyDay,-CohortDay,-CohortID,-CycleNum,-CohortStartTime) %>% 
  arrange(StudyDay,-value) %>% 
  group_by(key) %>% 
  distinct(StudyDay,.keep_all=TRUE)

View(indInfoLong %>% filter(key=="Remaining",CohortID==1))
View(indInfoLong %>% filter(CohortID==1))

ggplot(indInfoLong %>% filter(key=="Remaining",CohortID==1),aes(StudyDay,value,group=key,color=key)) +
  geom_line() +
  facet_wrap(~CohortID)


ggplot(indInfoLong ,aes(StudyDay,value,group=key,color=key)) +
  geom_line() +
  facet_grid(CohortID~CycleNum)
