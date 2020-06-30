#' @title calcClusterBase
#' @description Reads a series of MAgPIE files and combines them to a matrix
#' which is then used for calculating a clustering.
#' @return A matrix containing the data
#' @author Jan Philipp Dietrich
#' @seealso \code{\link{calcCluster}}
#' @importFrom magclass wrap read.magpie
#' @importFrom madrat toolMappingFile
calcClusterBase <- function() {

  d <- list()
  # read in data which should be used to determine cluster
  d$yld    <- calcOutput("Yields",     years=1995, aggregate=FALSE)
  d$airrig <- calcOutput("Irrigation", years=1995, aggregate=FALSE)
  d$td     <- calcOutput("TransportDistance",      aggregate=FALSE)[,,rep(1,floor(ndata(d$yld)/2))]

  cdata <- do.call(cbind,lapply(d,wrap,list(1,c(2,3))))
  cdata <- cdata[,dimSums(cdata,dim=1)!=0]
  cdata <- scale(cdata)
  colnames(cdata) <- paste0("i",1:ncol(cdata))

  #make sure that spatial dimension contains country information
  iso <- toolMappingFile("cell", "CountryToCellMapping.csv", readcsv = TRUE)$iso
  dimnames(cdata)[[1]] <- paste(iso, 1:length(iso), sep = ".")

  return(list(
    x = as.magpie(cdata),
    weight = NULL,
    unit = "1",
    description = "Similarity matrix as basis for clustering",
    isocountries = FALSE))
}
