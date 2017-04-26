#' getCohortNames
#'
#' @description Extracts the names of the cohorts in the CohortList is a FAIRsimulator \code{study} object.
#' @param StudyObj A FAIRsimulator \code{study} object.
#' @return A character vector
#' @details This is an internal function that extracts the values of the \code{Name} slot in the \code{CohortList} list of a FAIRsimulator \code{study} object. 
#' It is typically used to set the names of the \code{CohortList} components, so that they can be extracted using \code{names(CohortList)} and referred to by name, e.g. 
#' \code{CohortList[[cohortName]]}. 
#' @export
#'
#' @examples
#' \dontrun{
#' names(StudyObj$CohortList) <- getCohortNames(StudyObj)
#' }
getCohortNames <- function(StudyObj) {
  if(class(StudyObj) != "study") stop("studyObj needs to be a study object")
  return(cohorts(StudyObj) %listmap% "Name")
}
