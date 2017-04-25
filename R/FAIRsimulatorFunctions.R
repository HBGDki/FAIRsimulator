

##Accessor for CohortList
#' cohorts
#' @description Extract the list of cohorts from a FAIRsimulator \code{study} object
#' @param obj A FAIRsimulator \code{study} object
#'
#' @return The list of cohorts in the \code{study} object.
#' @export
#'
#' @examples
#' \dontrun{
#' myCohorts <- cohorts(StudyObj)
#' }
cohorts <- function(obj) {
  if(class(obj) != "study") stop("obj needs to be a study object")
  return(obj[["CohortList"]])
}

##Accessor for SubjectsList
#' subjects
#'
#' @description Extract the list of subjects from a FAIRsimulator \code{cohort} object
#' @param obj A FAIRsimulator \code{cohort} object
#'
#' @return The list of cohorts in the \code{cohort} object.
#' @export
#'
#' @examples
#' \dontrun{
#' mySubjects <- subjects(StudyObj)
#' }
subjects <- function(obj) {
  if(class(obj) != "cohort") stop("obj needs to be a cohort object")
  return(obj[["SubjectList"]])
}


#' getSubjectLevelItems
#'
#' @description Extract scalar subject level items from all subjects in all cohorts in a FAIRsimulator \code{study} object
#'
#' @param studyObj A FAIRsimulator \code{study} object
#' @param items A vector of the scalar subject level items to extract for each subject in each cohort.
#'
#' @return A data.frame with the columns: cohort name, SubjID and each of the entries in \code{item}.
#' @export
#'
#' @examples
#' \dontrun{
#' tmp <- getSubjectLevelItem(StudyObj,item=c("AgeAtRand","RandCohortTime"))
#' }
getSubjectLevelItems <- function(studyObj,items="AgeAtRand") {
  if(class(studyObj) != "study") stop("obj needs to be a study object")

  myCohorts <- cohorts(studyObj)
  nams      <- names(myCohorts)
  df <- data.frame()

  for(i in 1:length(myCohorts)) {
    mySubjects <- subjects(myCohorts[[i]])

    dfSubjCohort <- data.frame(Cohort=nams[i],ID=mySubjects %map% "StudyID")

    for(item in c(items)) {
      dfSubjCohort[item] <- mySubjects %map% item
    }

    if(nrow(df)==0) {
      df <-  dfSubjCohort
    } else {
      df <- rbind(df,dfSubjCohort)
    }
  }

  return(df)
}

#' map
#' @description An operator that extracts elements from a lisy
#' @param x The list to extract elements from
#' @param n The name of the element
#'
#' @return The element
#' @export
#'
#' @examples
#' \dontrun{
#' mySubjects %map% item
#' }
`%map%` <- function(x, n) {
sapply(x, `[[`, n)
}


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
                                          "CurrentCohortTime","TreatmentIndex","Treatment","TreatmentEff","RandNum"),
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


#' getItemsFromSubjects
#'
#' @description Extracts data from all subjects in a FAIRsimulator cohort object
#'
#' @param cohortObj A FAIRsimulator cohort object
#' @param ... Parameters passed to \code{\link{getSubjectItems}}
#'
#' @return A data frame
#' @seealso \code{\link{getSubjectItems}}
#' @export
#'
#' @examples
#' \dontrun{
#' tmp <- getItemsFromSubjects(StudyObj[["CohortList"]][[1]])
#' }
getItemsFromSubjects <- function(cohortObj,...) {
  if(class(cohortObj) != "cohort") stop("cohortObj needs to be a cohort object")
  
  subjList <- subjects(cohortObj)
  # browser()
  resDf <- NULL
  for(subj in subjList) {
    if(is.null(resDf)) {
      resDf <- getSubjectItems(subj,...)
    } else {
      resDf <- rbind(resDf,getSubjectItems(subj,...))
    }
  }
  
  resDf$CohortName       <- cohortObj$Name
  resDf$CycleNum         <- cohortObj$CycleNum
  resDf$StartNum         <- cohortObj$StartNum
  resDf$RandomizationAge <- cohortObj$RandomizationAgeRange[1]
  
  return(resDf)
}

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
  resDf     <- NULL
  
  for(cohrt in myCohorts) {
    if(is.null(resDf)) {
      resDf <- getItemsFromSubjects(cohrt,...)
    } else {
      resDf <- rbind(resDf,getItemsFromSubjects(cohrt,...))
    }
    #resDf$CohortStartTime <- cohrt$CohortStartTime
  }
  
  ## Rename some columns to be easier to understand and create factors
  resDf <- resDf %>% rename(Age=SampleAge,CohortID=StartNum,HAZ=Data) %>% 
    mutate(Cycle  = factor(CycleNum,levels=sort(unique(CycleNum)),labels=paste("Cycle",sort(unique(CycleNum)))),
           Cohort = factor(CohortID,levels=sort(unique(CohortID)),labels=paste("Cohort",sort(unique(CohortID))))
           )
  
  return(resDf)
}


#' getCohortNames
#'
#' @description Extracts the names of the cohorts in the CohortList is a FAIRsimulator \code{study} object.
#' @param StudyObj A FAIRsimulator \code{study} object.
#' @return A character vector
#' @details This is an internal function that extracts the values of the \code{Name} slot in the \code{CohortList} list of a FAIRsimulator \code{study} object. 
#' It is typically used to set the names of the \code{CohortList} components, so that they can be extracted using \code{names(CohortList)} and referred to by name, e.g. 
#' \code{CohortList[[cohortName]]}. 
#' @export
#'
#' @examples
#' \dontrun{
#' names(StudyObj$CohortList) <- getCohortNames(StudyObj)
#' }
getCohortNames <- function(StudyObj) {
  if(class(StudyObj) != "study") stop("studyObj needs to be a study object")
  return(cohorts(StudyObj) %map% "Name")
}


#' getCohortDetails
#'
#' @description Extract scalar statistics for the cohorts in a FAIRsimulator \code{study} object.
#' @param StudyObj A FAIRsimulator \code{study} object.
#'
#' @return A data.frame with values
#' @export
#'
#' @examples
#' \dontrun{
#' getCohortDetails(StudyObj)
#' }
getCohortDetails <- function(StudyObj) {
  if(class(StudyObj) != "study") stop("studyObj needs to be a study object")
  
  myCohorts <- cohorts(StudyObj)
  df        <- data.frame(
    Name = myCohorts %map% "Name",
    NSubjects = myCohorts %map% "NumberOfRecruitedSubjects",
    MaxSubjects =  myCohorts %map% "MaxNumberOfSubjects",
    CohortStartTime = myCohorts %map% "CohortStartTime",
    CohortID =  myCohorts %map% "StartNum",
    CycleNum =  myCohorts %map% "CycleNum",
    RandomizationAge = (myCohorts %map% "RandomizationAgeRange")[1,],
    EndRandomizationAge = (myCohorts %map% "RandomizationAgeRange")[2,]
  )
  
  row.names(df) <- NULL
  return(df)
}

#' plotHAZ
#' 
#' @description Plots the simulated data in a FAIRsimulator \code{study} object.
#' @param StudyObj A FAIRsimulator \code{study} object. 
#' @return A \code{ggplot} object
#' @export
#'
#' @examples
#' \dontrun{
#' plotHAZ(StudyObj)
#' }
plotHAZ <- function(StudyObj) {
  allData <- getAllSubjectData(StudyObj,scalarItems="AgeAtRand",covariates=NULL) 
  
  p1 <- ggplot(allData,aes(Age/30,HAZ,group=StudyID,color=Cycle)) +
    geom_point() +
    geom_line() +
    labs(color=NULL) + 
    theme(legend.position="top") +
    xlab("Age (months)") +
    facet_wrap(~Cohort)
  
  return(p1)
}

