#' getRandomizationProbabilities
#' @description Extract randomization probabilities over time/treatment
#' @param StudyObj A FAIRsimulator study object
#'
#' @return A data frame
#' @export
#'
#' @examples
#' \dontrun{
#' }
getRandomizationProbabilities<-function(StudyObj) {
  df<-data.frame()
  
  for (i in 1:length(StudyObj$CohortList)) {
    cohort<-StudyObj$CohortList[[i]]
    if (!is.null(cohort$PreviousRandomizationProbabilities)) {
      for (j in 1:length(cohort$PreviousRandomizationProbabilities)) {
        PrevProb<-cohort$PreviousRandomizationProbabilities[[j]]
        df<-rbind(df,data.frame(Name=cohort$Name,CohortTime=PrevProb$CohortTime,StudyTime=PrevProb$StudyTime,RandomizationProbabilities=PrevProb$RandomizationProbabilities,Treatments=PrevProb$Treatments))
      }
    }
    df<-rbind(df,data.frame(Name=cohort$Name,CohortTime=cohort$CurrentTime,StudyTime=0,RandomizationProbabilities=cohort$RandomizationProbabilities,Treatments=cohort$Treatments))
  }
  return(df)
}