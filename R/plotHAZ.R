#' plotHAZ
#' 
#' @description Plots the simulated data in a FAIRsimulator \code{study} object.
#' @param StudyObj A FAIRsimulator \code{study} object. 
#' @return A \code{ggplot} object
#' @export
#'
#' @examples
#' \dontrun{
#' plotHAZ(StudyObj)
#' }
plotHAZ <- function(StudyObj) {
  allData <- getAllSubjectData(StudyObj,scalarItems="AgeAtRand",covariates=NULL) 
  
  p1 <- ggplot(allData,aes(Age/30,HAZ,group=StudyID,color=Cycle)) +
    geom_point() +
    geom_line() +
    labs(color=NULL) + 
    theme(legend.position="top") +
    xlab("Age (months)") +
    facet_wrap(~Cohort)
  
  return(p1)
}