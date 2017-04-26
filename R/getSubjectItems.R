#' getSubjectItems
#' @description Extracts all subject level data as a data.frame from a FAIRsimulator individual object
#' @param subjectObj A FAIRsimulator individual object
#' @param scalarItems The scalar items to extract from the individual object. Set to \code{NULL} if no scalar items are to be extracted.
#' @param covariates Name of the item that holds the subject specific vector of covariates. Set to \code{NULL} if covariates should not be extracted.
#' @param longitudinalItems Name of the longitudinal items to be extracted. Set to \code{NULL} if no no longitudina items are to be extracted.
#'
#' @return A data frame with nrows=length of the longitudinal items. The scalar and covariate items are replicated to match the longitudinal items.
#' @export
#'
#' @examples
#' \dontrun{
#' getSubjectItems(StudyObj[["CohortList"]][[1]][["SubjectList"]][[1]])
#' #' }
getSubjectItems <- function(subjectObj,
                            scalarItems=c("StudyID","AgeAtRand","DateAtRand","RandStudyTime","RandCohortTime","CurrentAge",
                                          "CurrentCohortTime","TreatmentIndex","Treatment","TreatmentEff","RandNum","DropoutStudyTime","DropoutCohortTime"),
                            covariates = "Covariates",
                            longitudinalItems = c("SampleAge","CohortSampleTime","StudySampleTime","Data")) {
  
  if(class(subjectObj) != "individual") stop("subjectObj needs to be a individual object")
  
  myDf <- data.frame(StudyID=subjectObj[["StudyID"]])
  
  ## Extract the scalar items
  if(!is.null(scalarItems)) {
    for(i in scalarItems) {
      myDf[1,i] <- subjectObj[[i]]
    }
  }
  
  ## Extract the covariates
  if(!is.null(covariates)) {
    covs <- subjectObj[[covariates]]
    myDf[names(covs)] <- covs
  }
  
  
  ## Extract the longitudinal items
  longItems <- data.frame()
  if(!is.null(longitudinalItems)) {
    
    for(i in longitudinalItems) {
      if(is.null(subjectObj[[i]])) {
        longItems[1,i] <- NA
      } else {
        longItems[1:length(subjectObj[[i]]),i] <- subjectObj[[i]]
      }
    }
    
    ## Add StudyId
    longItems$StudyID <- subjectObj[["StudyID"]]
  }
  
  ## Merge the data
  myRetDf <- left_join(myDf,longItems,by="StudyID")
  return(myRetDf)
  
}
