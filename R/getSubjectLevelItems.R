#' getSubjectLevelItems
#'
#' @description Extract scalar subject level items from all subjects in all cohorts in a FAIRsimulator \code{study} object
#'
#' @param studyObj A FAIRsimulator \code{study} object
#' @param items A vector of the scalar subject level items to extract for each subject in each cohort.
#'
#' @return A data.frame with the columns: cohort name, SubjID and each of the entries in \code{item}.
#' @export
#'
#' @examples
#' \dontrun{
#' tmp <- getSubjectLevelItem(StudyObj,item=c("AgeAtRand","RandCohortTime"))
#' }
getSubjectLevelItems <- function(studyObj,items="AgeAtRand") {
  if(class(studyObj) != "study") stop("obj needs to be a study object")
  
  myCohorts <- cohorts(studyObj)
  nams      <- names(myCohorts)
  df <- data.frame()
  
  for(i in 1:length(myCohorts)) {
    mySubjects <- subjects(myCohorts[[i]])
    
    dfSubjCohort <- data.frame(Cohort=nams[i],ID=mySubjects %listmap% "StudyID")
    
    for(item in c(items)) {
      dfSubjCohort[item] <- mySubjects %listmap% item
    }
    
    if(nrow(df)==0) {
      df <-  dfSubjCohort
    } else {
      df <- rbind(df,dfSubjCohort)
    }
  }
  
  return(df)
}