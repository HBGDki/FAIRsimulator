#' plotStudyCohorts
#'
#' @description Plots the cohort durations/cycles in a FAIRsimulator \code{study} object.
#' @param StudyObj A FAIRsimulator \code{study} object.
#' @param cohortlength The length (in days) of cohorts that are evolving 
#'
#' @return A \code{ggplot} object
#' @export
#'
#' @examples
#' \dontrun{
#' plotStudyCohorts(StudyObj)
#' }
plotStudyCohorts <- function(StudyObj,cohortlength=6*30,wrapVariable=NULL) {
  
  if(class(StudyObj) != "study") stop("studyObj needs to be a study object")
  
  cohortdetails<-getCohortDetails(StudyObj)
  
  cohortdetails$EndRandomizationAge[cohortdetails$CycleNum>1]<-cohortdetails$EndRandomizationAge[cohortdetails$CycleNum>1]+(cohortdetails$CycleNum[cohortdetails$CycleNum>1]-1)*cohortlength
  cohortdetails$RandomizationAge[cohortdetails$CycleNum>1]<-cohortdetails$RandomizationAge[cohortdetails$CycleNum>1]+(cohortdetails$CycleNum[cohortdetails$CycleNum>1]-1)*cohortlength
  
  p <- ggplot(data=cohortdetails)
  p <- p + geom_rect(aes(xmin=CohortStartTime/30,
                         xmax=(CohortStartTime+CohortDuration)/30,
                         ymin=RandomizationAge/30,
                         ymax=EndRandomizationAge/30,
                         fill=Cohort,linetype=Cohort),
                     color="black",alpha=0.7)
  p <- p + xlab("Study time (month)")
  p <- p + ylab("Age at randomization (months)")
  p <- p + labs(color=NULL,linetype=NULL,fill=NULL)
  p <- p + scale_y_continuous(breaks = seq(0,48,by=1))
  
  p <- p + scale_x_continuous(breaks = seq(0,48,by=6),limits= c(0, max((cohortdetails$CohortStartTime+cohortdetails$CohortDuration)/30)))
  if(!is.null(wrapVariable)) p <- p + facet_wrap(as.formula(paste("~", wrapVariable)))
  
  return(p)
}
