#' getCohortDetails
#'
#' @description Extract scalar statistics for the cohorts in a FAIRsimulator \code{study} object.
#' @param StudyObj A FAIRsimulator \code{study} object.
#'
#' @return A data.frame with values
#' @export
#'
#' @examples
#' \dontrun{
#' getCohortDetails(StudyObj)
#' }
getCohortDetails <- function(StudyObj) {
  if(class(StudyObj) != "study") stop("studyObj needs to be a study object")
  
  myCohorts <- cohorts(StudyObj)
  df        <- data.frame(
    Name = myCohorts %listmap% "Name",
    NSubjects = myCohorts %listmap% "NumberOfRecruitedSubjects",
    MaxSubjects =  myCohorts %listmap% "MaxNumberOfSubjects",
    CohortStartTime = myCohorts %listmap% "CohortStartTime",
    CohortID =  myCohorts %listmap% "StartNum",
    CycleNum =  myCohorts %listmap% "CycleNum",
    RandomizationAge = (myCohorts %listmap% "RandomizationAgeRange")[1,],
    EndRandomizationAge = (myCohorts %listmap% "RandomizationAgeRange")[2,],
    CohortDuration = myCohorts %listmap% "CurrentTime",
    Level = myCohorts %listmap% "Level"
  )
  
  row.names(df) <- NULL
  
  df <- df %>%  mutate(Cycle  = factor(CycleNum,levels=sort(unique(CycleNum)),labels=paste("Cycle",sort(unique(CycleNum)))),
                       Cohort = factor(CohortID,levels=sort(unique(CohortID)),labels=paste("Cohort",sort(unique(CohortID))))
  )
  
  return(df)
}
