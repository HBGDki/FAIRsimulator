indInfo <- extractInclusionData(StudyObj) 

indInfo %>% filter(CohortDay>=-1) %>% filter(Cycle=="Cycle 2",Cohort=="Cohort 2") %>% View

p1 <- ggplot(indInfo %>% filter(CohortDay>=-1),aes(StudyDay/30,Remaining)) +
  geom_line() +
  xlab("Study time (months)") +
  ylab("Number of active subjects") +
  scale_x_continuous(breaks = seq(0,48,by=6),limits=c(0,max((indInfo %>% filter(CohortID==1,CohortDay>=0))$StudyDay/30))) +
  facet_grid(Cohort~Cycle)

cohortDetails <- getCohortDetails(StudyObj)
allData       <- getAllSubjectData(StudyObj) 
maxStudyTime  <- max(allData$StudySampleTime)

allRes <- NULL
for(i in 1:nrow(cohortDetails)) {
   i <- 4
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