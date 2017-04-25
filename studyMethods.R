#library(tidyverse)
library(FAIRsimulator)
library(zoo)



tmp <- getSubjectItems(StudyObj[["CohortList"]][[1]][["SubjectList"]][[1]])
tmp <- getItemsFromSubjects(StudyObj[["CohortList"]][[1]])
tmp <- getAllSubjectData(StudyObj,scalarItems="AgeAtRand",covariates=NULL)

plotActiveSubjects(StudyObject)


