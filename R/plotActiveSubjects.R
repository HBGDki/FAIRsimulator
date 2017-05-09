#' plotActiveSubjects
#'
#' @description Plot the number of active subjects in each cohort and cycle.
#' 
#' @param StudyObject A FAIRsimulator \code{study} object.
#'
#' @details The function returns a graph of the number of subjects in each cohort and cycle. The number is computed as the 
#' cumulated number of recruted subjects minus the cumulated number subjects that 
#' have dropped out minus the cumulated number of subjects that have completed the cohort.
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' plotActiveSubjects(StudyObject)
#' }
plotActiveSubjects <- function(StudyObj) {
  
  indInfo <- extractInclusionData(StudyObj) 
  
  p1 <- ggplot(indInfo %>% filter(CohortDay>=-1),aes(StudyDay/30,Remaining)) +
    geom_line() +
    xlab("Study time (months)") +
    ylab("Number of active subjects") +
    scale_x_continuous(breaks = seq(0,48,by=6),limits=c(0,max((indInfo %>% filter(CohortID==1,CohortDay>=0))$StudyDay/30))) +
    facet_grid(Cohort~Cycle)
  
  return(p1)
}