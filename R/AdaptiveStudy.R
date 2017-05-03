#' AdaptiveStudy
#'
#' @description Function that runs a FAIRsimulator simulation.
#' @param StudyObj An initialized FAIRsimulator \code{study} object.
#'
#' @return A FAIRsimulator \code{study} object complete with all aspects generated during the simulations.
#' @export
#'
#' @examples
#' \dontrun{
#' StudyObj<-AdaptiveStudy(StudyObj)
#' }
AdaptiveStudy<-function(StudyObj){
  
  StudyObj<-StudyObj$InitEvent(StudyObj) #More thing to intialize before the loop starts

  while (StudyObj$StopEvent(StudyObj)==FALSE) { #Start the study
    
    if (!is.null(StudyObj$EventList)) {
      for (i in 1:length(StudyObj$EventList)) {
        StudyObj<-StudyObj$EventList[[i]](StudyObj) #Call all events in the event list
      }
    }
    StudyObj<-StudyObj$StudyIncrementEvent(StudyObj) #Increment study, e.g. update study day to study day + 1
  }
  
  names(StudyObj$CohortList) <- getCohortNames(StudyObj) # Set the names of the cohorts in CohortList based on their Name slots
  
  return(StudyObj) #Study simulations are finished, return StudyObj
}