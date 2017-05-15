#' getAllSubjectData
#'
#' @description Extracts data framm all subjects in all cohorts in a FAIRsimulator study object
#' @param StudyObj A FAIRsimulator study object
#' @param ... Parameters passed to \code{\link{getSubjectItems}}
#' @seealso \code{\link{getSubjectItems}}, \code{\link{getItemsFromSubjects}}
#' @return A data frame
#' @export
#'
#' @examples
#' \dontrun{
#' tmp <- getAllSubjectData(StudyObj)
#' }
getAllSubjectData <- function(StudyObj,...) {
  
  if(class(StudyObj) != "study") stop("studyObj needs to be a study object")
  
  myCohorts <- cohorts(StudyObj)
  i<-1
  dataList<-list()
  for(cohrt in myCohorts) {
    dataList[[i]]<-getItemsFromSubjects(cohrt,...)
    i<-i+1
  }
  
  resDf<-as.data.frame(data.table::rbindlist(dataList))
  
  ## Rename some columns to be easier to understand and create factors
  resDf <- resDf %>% rename(Age=SampleAge,CohortID=StartNum,HAZ=Data) %>% 
    mutate(Cycle  = factor(CycleNum,levels=sort(unique(CycleNum)),labels=paste("Cycle",sort(unique(CycleNum)))),
           Cohort = factor(CohortID,levels=sort(unique(CohortID)),labels=paste("Cohort",sort(unique(CohortID))))
    )
  
  return(resDf)
}