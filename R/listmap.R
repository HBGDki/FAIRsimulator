#' listmap
#' @description An operator that extracts elements from a lisy
#' @param x The list to extract elements from
#' @param n The name of the element
#'
#' @return The element
#' @export
#'
#' @examples
#' \dontrun{
#' mySubjects %listmap% item
#' }
`%listmap%` <- function(x, n) {
sapply(x, `[[`, n)
}
















