#' plotStudyCohorts
#'
#' @description Plots the cohort durations/cycles in a FAIRsimulator \code{study} object.
#' @param StudyObj A FAIRsimulator \code{study} object.
#' @param cohortlength The length (in days) of cohorts that are evolving 
#' @param plotAnaTimes Logical. Should the analysis times be indicated in the plot.
#' @param wrapVariable The name (as a string) of the variable to create facets with, using \code{facet_wrap}.
#' @param shiftByLevel Shifting each cohort upwards the y-axis by level and this amount, default= NULL (no shift)
#'
#' @return A \code{ggplot} object
#' @export
#'
#' @examples
#' \dontrun{
#' plotStudyCohorts(StudyObj)
#' }
plotStudyCohorts <- function(StudyObj,cohortlength=6*30,plotAnaTimes=FALSE,wrapVariable=NULL,shiftWithinLevel=NULL) {
  
  if(class(StudyObj) != "study") stop("studyObj needs to be a study object")
  cohortdetails<-getCohortDetails(StudyObj)
  
  cohortdetails$EndRandomizationAge[cohortdetails$CycleNum>1]<-cohortdetails$EndRandomizationAge[cohortdetails$CycleNum>1]+(cohortdetails$CycleNum[cohortdetails$CycleNum>1]-1)*cohortlength
  cohortdetails$RandomizationAge[cohortdetails$CycleNum>1]<-cohortdetails$RandomizationAge[cohortdetails$CycleNum>1]+(cohortdetails$CycleNum[cohortdetails$CycleNum>1]-1)*cohortlength
  
  cohortdetails$ymin<-cohortdetails$RandomizationAge/30
  cohortdetails$ymax<-cohortdetails$EndRandomizationAge/30
  
  
  if (!is.null(shiftWithinLevel)) {
    for (i in 1:max(cohortdetails$Level)) {
      if (nrow(cohortdetails[cohortdetails$Level==i,])>0) {
        cohortdetails$ymin[cohortdetails$Level==i]<-cohortdetails$ymin[cohortdetails$Level==i]+(0:(nrow(cohortdetails[cohortdetails$Level==i,])-1))*shiftWithinLevel
        cohortdetails$ymax[cohortdetails$Level==i]<-cohortdetails$ymax[cohortdetails$Level==i]+(0:(nrow(cohortdetails[cohortdetails$Level==i,])-1))*shiftWithinLevel
      }
    }
  }
  
  p <- ggplot(data=cohortdetails)
  p <- p + geom_rect(aes(xmin=CohortStartTime/30,
                         xmax=(CohortStartTime+CohortDuration)/30,
                         ymin=ymin,
                         ymax=ymax,
                         fill=Cohort,linetype=Cohort),
                     color="black",alpha=0.7)
  p <- p + xlab("Study time (month)")
  p <- p + ylab("Age at randomization (months)")
  p <- p + labs(color=NULL,linetype=NULL,fill=NULL)
  p <- p + scale_y_continuous(breaks = seq(0,48,by=1))
  
  p <- p + scale_x_continuous(breaks = seq(0,48,by=6),limits= c(0, max((cohortdetails$CohortStartTime+cohortdetails$CohortDuration)/30)))
  
  if(plotAnaTimes) {
    
    anaTime        <- sapply(cohorts(StudyObj), `[[`, "AnalysisTime",simplify=FALSE)
    
    ## Need to handle when there are ab unequal number if IAs in different cohorts
    maxIA <- max(unlist(lapply(anaTime,function(x) length(x)))) # Count the max number of IAs
    
    myFun <- function(x,maxIA) {
      numElem <- length(x)
      
      if(is.null(x)) {
        x <- NA
      } else if(numElem < maxIA) {
        x[(numElem+1):maxIA] <- x[numElem]
      }
      return(x)
    }
    
    anaTime <- data.frame(lapply(anaTime,FUN=myFun,maxIA))
    
    #anaTime        <- data.frame(lapply(anaTime,function(x) {ifelse(is.null(x),return(NA),return(x))}))
    #anaTime        <- data.frame(sapply(cohorts(StudyObj), `[[`, "AnalysisTime",simplify=FALSE))
    names(anaTime) <- names(sapply(cohorts(StudyObj), `[[`, "AnalysisTime",simplify=FALSE))
    anaTime        <- gather(anaTime,"Name","AnalysisTimes") %>%  distinct(Name,AnalysisTimes)
    anaTime$Name   <- as.factor(anaTime$Name) 
    anaTime        <- left_join(anaTime,cohortdetails,by="Name") %>% mutate(midTime = (ymin + ymax)/2)
    
    p <- p + geom_point(data=anaTime,aes(AnalysisTimes/30,midTime),shape=19,size=4)
    p <- p + geom_point(data=anaTime,aes(AnalysisTimes/30,midTime,color=Cohort),shape=17,size=2)
  }
  
  if(!is.null(wrapVariable)) p <- p + facet_wrap(as.formula(paste("~", wrapVariable)))
  
  return(p)
}
