#' @title correctGLW3
#' @description Read GLW3 file
#' @return Magpie objects with results on cellular level, weight, unit and description.
#' @param x magpie object provided by the read function
#' @author Marcos Alves
#' @seealso
#'   \code{\link{readGLW3}}
#' @examples
#'
#' \dontrun{
#'   readSource("GLW3", subtype = "DA", convert="onlycorrect")
#' }
#'
#' @importFrom madrat toolConditionalReplace

correctGLW3<- function(x){

  x <- toolConditionalReplace(x, conditions = c("is.na()","<0"), replaceby = 0)

  return(x)
}
