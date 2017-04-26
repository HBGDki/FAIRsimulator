##Accessor for SubjectsList
#' subjects
#'
#' @description Extract the list of subjects from a FAIRsimulator \code{cohort} object
#' @param obj A FAIRsimulator \code{cohort} object
#'
#' @return The list of cohorts in the \code{cohort} object.
#' @export
#'
#' @examples
#' \dontrun{
#' mySubjects <- subjects(StudyObj)
#' }
subjects <- function(obj) {
  if(class(obj) != "cohort") stop("obj needs to be a cohort object")
  return(obj[["SubjectList"]])
}
