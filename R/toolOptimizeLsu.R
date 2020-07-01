#' Tool: StoreMapping
#'
#' Function which creates a mapping file
#'
#'
#' @param map mapping data to be stored.
#' @param name File name of the mapping file. Supported file types are currently csv (, or ; separated)
#' and rda (which needs to have the data stored with the object name "data"!). Use code{\link{toolConvertMapping}}
#' to convert between both formats
#' @param type Mapping type (e.g. "regional", "cell", or "sectoral"). Can be set to NULL if file
#' is not stored in a type specific subfolder
#' @param error.existing Boolean which decides whether an error is returned if
#' the mapping file does not exist or not.
#' @author Jan Philipp Dietrich
#' @seealso \code{\link{calcOutput}}, \code{\link{toolConvertMapping}}
#' @importFrom tools file_ext
#' @importFrom utils write.table
#' @importFrom magclass is.magpie
#' @importFrom madrat getConfig
#' @export


# toolOptimizeLsu <- function(test_data, model,min_Lsu, max_Lsu,col_means,col_stddevs, range){
toolOptimizeLsu <- function(envi_data, model_h5, lsu_bounds = c(-2,2), inputs_vec = NULL) {
  if (attr(envi_data, "scaled:center")) {


    # for (i in 1:length(inputs_vec)) {
    #   print(grep(inputs_vec[i], append(colnames(envi_data),"lsu"), ignore.case=TRUE))
    # }
    #
    # lsu_index <- grep(paste(append(colnames(envi_data),"lsu"),collapse="|"), inputs_vec , ignore.case=TRUE, value=TRUE)
    # x <- as.data.frame(cbind(envi_data, "lsu" = 0))
    # x <- x[,lsu_index]
    # colnames(x)

    envi_data <- cbind(envi_data, "lsu" = 0)

    min_Lsu <- lsu_bounds[1]
    max_Lsu <- lsu_bounds[2]


    optimization <- function(envi_data) {
      network_lsu <- function(lsu) {
        envi_data["lsu"] <- lsu
        y <- predict(model_h5, t(envi_data))
        return(y)
      }

      opt <- optimize(network_lsu, interval = lsu_bounds , maximum = T)
      return(opt)
    }


    # optmizing all cells using apply function
    max <- apply(envi_data, 1, optimization)
  } else {
    stop("The cellular dataset must be scaled and in matrix format")
  }
  return(max)
}
