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
  # browser()
  resDf <- NULL
  for(subj in subjList) {
    if(is.null(resDf)) {
      resDf <- getSubjectItems(subj,...)
    } else {
      resDf <- rbind(resDf,getSubjectItems(subj,...))
    }
  }
  
  resDf$CohortName       <- cohortObj$Name
  resDf$CycleNum         <- cohortObj$CycleNum
  resDf$StartNum         <- cohortObj$StartNum
  resDf$RandomizationAge <- cohortObj$RandomizationAgeRange[1]
  
  return(resDf)
}