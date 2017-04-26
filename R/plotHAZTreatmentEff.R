#' plotHAZTreatmentEff
#' 
#' @description Plots the simulated data in a FAIRsimulator \code{study} object with treatment effects.
#' @param StudyObj A FAIRsimulator \code{study} object. 
#' @details Will produce a graph of individual change from cohort baseline HAZ. Treatment effects are visualised by linear regression lines
#' for each treatment group.
#' @return A \code{ggplot} object
#' @export
#'
#' @examples
#' \dontrun{
#' plotHAZTreatmentEff(StudyObj)
#' }
plotHAZTreatmentEff <- function(StudyObj) {

  allData <- getAllSubjectData(StudyObj,covariates=NULL) 
  allData$Treatment <- factor(allData$Treatment,levels=sort(unique(allData$Treatment)),labels=paste("Treatment",sort(unique(allData$Treatment))))
  allData <- allData %>% group_by(StudyID,Cycle) %>% mutate(HAZcfb=HAZ-HAZ[1])
  
  p1 <- ggplot(allData,aes(CohortSampleTime/30,HAZcfb,group=StudyID)) +
    geom_point(color="lightgrey") +
    geom_line(color="lightgrey") +
    geom_smooth(aes(group=Treatment,color=Treatment),se=F,method="lm")+
    labs(color=NULL) + 
    theme(legend.position="top") +
    xlab("Cohort time (months)") +
    ylab("HAZ - change from individual cohort baseline")+
    facet_grid(Cycle~Cohort)
  
  return(p1)
}