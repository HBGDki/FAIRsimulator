#' extractInclusionData
#'
#' @description Extracts statistics on recruitment, dropout and completion of subjetcs from each cohort.
#' @param StudyObj A FAIRsimulator \code{study} object
#'
#' @return A \code{data.frame} with columns: StudyDay, Recruited, Dropout, Completed, Remaining, CohortDay, CohortID, CycleNum, CohortStartTime, Cycle and Cohort.
#' @export
#'
#' @examples
#' \dontrun{
#' indInfo <- extractInclusionData(StudyObj) 
#' }
extractInclusionData <- function(StudyObj) {
  cohortDetails <- getCohortDetails(StudyObj)
  allData       <- getAllSubjectData(StudyObj) 
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
      dplyr::rename(Day=RandStudyTime)
    
    dropout <- cohort1 %>% 
      filter(!is.na(DropoutStudyTime)) %>% 
      group_by(DropoutStudyTime) %>% 
      distinct(StudyID) %>%  
      tally %>% 
      ungroup %>% 
      mutate(Dropout = cumsum(n)) %>% 
      dplyr::select(-n) %>% 
      dplyr::rename(Day=DropoutStudyTime)
    
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
      dplyr::rename(Day = StudySampleTime)
    
    
    resDf <- data.frame(Day=0:maxStudyTime)
    resDf <- left_join(resDf,recruited,by="Day") %>% na.locf %>% mutate(Recruited = ifelse(is.na(Recruited),0,Recruited))
    resDf <- left_join(resDf,dropout,by="Day")   %>% na.locf %>% mutate(Dropout = ifelse(is.na(Dropout),0,Dropout))
    resDf <- left_join(resDf,completed,by="Day") %>% na.locf %>% mutate(Completed = ifelse(is.na(Completed),0,Completed))
    
    resDf <- resDf %>% 
      mutate(Remaining = Recruited-Dropout-Completed,
             CohortDay = Day - cohortStartTime,CohortID=cohortID,
             CycleNum=cycleNum,CohortStartTime=cohortStartTime) %>% 
      dplyr::rename(StudyDay = Day)
    
    if(is.null(allRes)) {
      allRes <- resDf
    } else {
      allRes <- bind_rows(allRes,resDf)
    }
    
  }
  
  allRes <- allRes %>% 
    mutate(Cycle  = factor(CycleNum,levels=sort(unique(CycleNum)),labels=paste("Cycle",sort(unique(CycleNum)))),
           Cohort = factor(CohortID,levels=sort(unique(CohortID)),labels=paste("Cohort",sort(unique(CohortID)))))
  
  return(allRes)
  
}
