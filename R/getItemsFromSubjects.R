#' getItemsFromSubjects
#'
#' @description Extracts data from all subjects in a FAIRsimulator cohort object
#'
#' @param cohortObj A FAIRsimulator cohort object
#' @param ... Parameters passed to \code{\link{getSubjectItems}}
#'
#' @return A data frame
#' @seealso \code{\link{getSubjectItems}}
#' @export
#'
#' @examples
#' \dontrun{
#' tmp <- getItemsFromSubjects(StudyObj[["CohortList"]][[1]])
#' }
getItemsFromSubjects <- function(cohortObj,...) {
  if(class(cohortObj) != "cohort") stop("cohortObj needs to be a cohort object")
  
  subjList <- subjects(cohortObj)
  dataList<-list()
  i<-1
  for(subj in subjList) {
    dataList[[i]]<-getSubjectItems(subj,...)
    i<-i+1
  }
  resDf<-as.data.frame(data.table::rbindlist(dataList))

  if(nrow(resDf) == 0) return(NULL) # If there aren't any subjects in the cohort

  

  resDf$CohortName       <- cohortObj$Name
  resDf$CycleNum         <- cohortObj$CycleNum
  resDf$StartNum         <- cohortObj$StartNum
  resDf$RandomizationAge <- cohortObj$RandomizationAgeRange[1]
  resDf$Level            <- cohortObj$Level
  resDf$CohortStartTime  <- cohortObj$CohortStartTime

  return(resDf)
}