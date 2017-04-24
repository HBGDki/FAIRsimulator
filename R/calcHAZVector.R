#' calHAZVector
#'
#' @description Computes a vector of individual HAZ values based on the FREM model.
#' @param age A vector of ages (in years) for which to compute the HAZ values.
#' @param basethetas A vector of base model thetas (fixed effects).
#' @param covthetas A vector with the FFEM, parameter specific covariate effetcs.
#' @param etas A vecor of individual specific random effects.
#' @param BASE The fixed effects parameter value for the BASE parameter in the FREM  model.
#' @param PLMAX The fixed effects parameter value for the PLMAX parameter in the FREM  model.
#' @param HLKON The fixed effects parameter value for the HLKON parameter in the FREM  model.
#' @param HLKOFF The fixed effects parameter value for the HLKOFF parameter in the FREM  model.
#' @param BASSL The fixed effects parameter value for the BASSL parameter in the FREM  model.
#' @param BP The fixed effects parameter value for the BP parameter in the FREM  model.
#' @param BASECOV The FFEM covariate effect for the BASE parameter in the FREM  model.
#' @param PLMAXCOV The FFEM covariate effect for the PLMAX parameter in the FREM  model.
#' @param HLKONCOV The FFEM covariate effect for the HLKON parameter in the FREM  model.
#' @param HLKOFFCOV The FFEM covariate effect for the HLKOFF parameter in the FREM  model.
#' @param BASSLCOV The FFEM covariate effect for the BASSL parameter in the FREM  model.
#' @param BPCOV The FFEM covariate effect for the BP parameter in the FREM  model.
#'
#' @details The basethetas+covthetas+etas and BASE+PLMAX+HLKON+HLKOFF+BASSL+BP+(the corrsponding covariate effects) are mutually exclusive. 
#' If basethetas is provided, then the former set of arguments will be used.
#' @return A vector of HAZ values as long as the vector of \code{age}.
#' @export
#' 
#' @examples
#' \dontrun{
#' hazValues <- calcHAZVector(age = c(0,30,60,90)/(30*12),basethetas = thbasevector,covthetas = FREMCoeff,etas = IndEtas)
#' }
calcHAZVector <- calcHAZ <- function(age,basethetas=NULL,covthetas=rep(0,length(basethetas)),etas=rep(0,length(basethetas)),
                          BASE,PLMAX,HLKON,HLKOFF,BASSL,BP,
                          BASECOV=0,PLMAXCOV=0,HLKONCOV=0,HLKOFFCOV=0,BASSLCOV=0,BPCOV=0) {

  if(!is.null(basethetas)) {
    BASE    <- basethetas[1]          + covthetas[1] + etas[1]
    PLMAX   <- basethetas[2]          + covthetas[2] + etas[2]
    HLKON   <- exp(log(basethetas[3]) + covthetas[3] + etas[3])
    HLKOFF  <- exp(log(basethetas[4]) + covthetas[4] + etas[4])
    BASSL   <- basethetas[5]          + covthetas[5] + etas[5]
    BP      <- exp(log(basethetas[6])      + covthetas[6] + etas[6])
  } else {
    BASE    <- BASE            + BASECOV
    PLMAX   <- PLMAX           + PLMAXCOV
    HLKON   <- exp(log(HLKON)  + HLKONCOV)
    HLKOFF  <- exp(log(HLKOFF) + HLKOFFCOV)
    BASSL   <- BASSL           + BASSLCOV
    BP      <- exp(log(BP)          + BPCOV)
  }
  
  KON  <- log(2)/HLKON
  KOFF <- log(2)/HLKOFF
  
  BASFAC1   <- BASSL*age
  BASFAC2   <- BASSL*BP
  
  PHI      <- age**2.25/(age**2.25 + BP**2.25)
  
  MYBASE   <- BASE + (1-PHI)*BASFAC1 + PHI*BASFAC2
  
  IPRED <-  MYBASE-PLMAX*KON*(exp(-KOFF*age)-exp(-KON*age))/(KON-KOFF)

  return(IPRED)
}
  

