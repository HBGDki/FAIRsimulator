#' plotHAZTreatmentEff
#' 
#' @description Plots the simulated data in a FAIRsimulator \code{study} object with treatment effects.
#' @param StudyObj A FAIRsimulator \code{study} object. 
#' @param treatmentNames A vector of treatment names. Default is the treatment names assigned at study initialisation.
#' @param treatmentColors A vector of colors for the treatment effects. Has to be the same length as the \code{treatmentNames}.
#' @param SoCtrt A vector of indicies for the \code{treatmentNames} indicating Standard of Care treatment arms. The default is set by grepping \code{treatmentNames} for the string 'SoC'.
#' @param SoCcol A vector of color names to assign to the SoC treatment arms.
#' @details Will produce a graph of individual change from cohort baseline HAZ. Treatment effects are visualised by linear regression lines
#' for each treatment group.
#' @return A \code{ggplot} object
#' @export
#'
#' @examples
#' \dontrun{
#' plotHAZTreatmentEff(StudyObj)
#' }
plotHAZTreatmentEff <- function(StudyObj,
                                treatmentNames = unique(unlist(StudyObj$StudyDesignSettings$Treatments)),
                                treatmentColors = ggplotColours(length(treatmentNames)),
                                SoCtrt = grep(pattern = "SoC",unique(unlist(StudyObj$StudyDesignSettings$Treatments))),
                                SoCcol = "black") {


  if(length(treatmentNames) != length(treatmentColors)) stop("The number of treatment names myst be the same as the number of treatment colors.")
  names(treatmentColors) <- treatmentNames
  
  treatmentColors[SoCtrt] <- SoCcol
  
  allData <- getAllSubjectData(StudyObj,covariates=NULL) 
  allData <- allData %>% group_by(StudyID,Cycle) %>% mutate(HAZcfb=HAZ-HAZ[1])
  
  p1 <- ggplot(allData,aes(CohortSampleTime/30,HAZcfb,group=StudyID)) +
    geom_point(color="lightgrey") +
    geom_line(color="lightgrey") +
    geom_smooth(aes(group=Treatment,color=Treatment),se=F,method="lm")+
    labs(color=NULL) + 
    theme(legend.position="top") +
    xlab("Cohort time (months)") +
    ylab("HAZ - change from individual cohort baseline")+
    facet_grid(Cohort~Cycle) + 
    scale_color_manual(values=treatmentColors)
  
  return(p1)
}