##Accessor for CohortList
#' cohorts
#' @description Extract the list of cohorts from a FAIRsimulator \code{study} object
#' @param obj A FAIRsimulator \code{study} object
#'
#' @return The list of cohorts in the \code{study} object.
#' @export
#'
#' @examples
#' \dontrun{
#' myCohorts <- cohorts(StudyObj)
#' }
cohorts <- function(obj) {
  if(class(obj) != "study") stop("obj needs to be a study object")
  return(obj[["CohortList"]])
}
