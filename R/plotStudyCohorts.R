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

plotStudyCohorts <- function(StudyObj,cohortlength=6*30,plotAnaTimes=FALSE,plotFutilityTimes=FALSE,plotFinalAnalysis=FALSE,wrapVariable=NULL,shiftWithinLevel=NULL,includeCohortStopTime=TRUE) {

  
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
    
    anaTime        <- data.frame(lapply(anaTime,FUN=myFun,maxIA))
    names(anaTime) <- names(sapply(cohorts(StudyObj), `[[`, "AnalysisTime",simplify=FALSE))
    anaTime        <- gather(anaTime,"Name","AnalysisTimes") %>%  distinct(Name,AnalysisTimes)
    
    ## Should we remove the analyses done at the end of cohorts
    if(!includeCohortStopTime) {
      cohortStopTimes      <- data.frame(StopTime =cohorts(StudyObj) %listmap% "CurrentTime" + cohorts(StudyObj) %listmap% "CohortStartTime")
      cohortStopTimes$Name <- row.names(cohortStopTimes)
      
      anaTime <- left_join(anaTime,cohortStopTimes,by="Name") %>% mutate(AnalysisTimes =ifelse(AnalysisTimes == StopTime,NA,AnalysisTimes)) %>% select(-StopTime) # filter(AnalysisTimes != StopTime) %>% select(-StopTime)
    }
    
    anaTime$Name   <- factor(anaTime$Name,levels=levels(cohortdetails$Name)) 
    
    yshift <- 0
    if(plotFutilityTimes) yshift = 0.25
    anaTime        <- left_join(anaTime,cohortdetails,by="Name") %>% mutate(midTime = (ymin + ymax)/2+yshift)
    
    p <- p + geom_point(data=anaTime,aes(AnalysisTimes/30,midTime),shape=19,size=4)
    p <- p + geom_point(data=anaTime,aes(AnalysisTimes/30,midTime,color=Cohort),shape=17,size=2)
  }
  

  ## Plot the futility time
  if(plotFutilityTimes) {
    futilityTime <- data.frame(futTime = StudyObj$FutilityList %listmap% "StudyTime",
                               Name = StudyObj$FutilityList %listmap% "CohortName",stringsAsFactors = FALSE)
    
    if(!includeCohortStopTime) {
      cohortStopTimes      <- data.frame(StopTime =cohorts(StudyObj) %listmap% "CurrentTime" + cohorts(StudyObj) %listmap% "CohortStartTime")
      cohortStopTimes$Name <- row.names(cohortStopTimes)
      
      futilityTime <- left_join(futilityTime,cohortStopTimes,by="Name") %>% mutate(futTime =ifelse(futTime == StopTime,NA,futTime)) %>% select(-StopTime) # filter(AnalysisTimes != StopTime) %>% select(-StopTime)
    }
    
    futilityTime$Name   <- factor(futilityTime$Name,levels=levels(cohortdetails$Name)) 
    
    futilityTime        <- left_join(futilityTime,cohortdetails,by="Name") %>% mutate(midTime = (ymin + ymax)/2-0.25)
    
    
    p <- p + geom_point(data=futilityTime,aes(futTime/30,midTime),shape=15,size=4)
    p <- p + geom_point(data=futilityTime,aes(futTime/30,midTime,color=Cohort),shape=20,size=2)
  }
  
  ## Plot final analysis
  if(plotFinalAnalysis) {
    p <- p + geom_vline(xintercept = StudyObj$CurrentTime/30,linetype="dashed")
  }
  
  if(!is.null(wrapVariable)) p <- p + facet_wrap(as.formula(paste("~", wrapVariable)))
  
  return(p)
}

