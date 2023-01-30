#' @title readGrassSoilEmu
#' @description Read files related to the training and optimization of the LPJml emulators.
#' @param subtype Subtype of file to be opened. Subtypes available:
#' 'weights', 'inputs',  'stddevs' and 'means'.
#' @return Magpie objects with a diverse information
#' @author Marcos Alves
#' @examples
#' \dontrun{
#' readSource("GrassSoilEmu",
#'            subtype = "ISIMIP3bv2:IPSL_CM6A_LR:ssp126:1965_2100:5f5fa2:weights",
#'            convert = FALSE)
#' }
#'
#' @import madrat
#' @importFrom stringr str_c
readGrassSoilEmu <- function(subtype = "ISIMIP3bv2:IPSL_CM6A_LR:ssp126:1965_2100:5f5fa2:stddevs_lab") {
  subtypeSplit <- toolSplitSubtype(subtype, list(version = NULL, climatemodel = NULL, scenario = NULL,
                                                 years = NULL, model = NULL, variable = NULL))
  file <- subtypeSplit$variable
  dirs <- list.dirs()
  folder <- grep(subtypeSplit$model, dirs, value = TRUE)
  if (length(dir.exists(file.path(folder))) != 0) {
    filesList <- list.files(folder)
    files <- filesList[grep(file, filesList)]

    if (file %in% "weights") {
      files <- files[grep(".rds", files, ignore.case = TRUE)]
      x <- readRDS(file.path(folder, files))
      xDims <- lapply(x, dim)
      xNames <- NULL
      for (i in seq_along(xDims)) {
        xNames[i] <- str_c(xDims[[i]], collapse = "_")
      }
      names(x) <- paste(paste0("l", seq_along(x)), paste0(xNames, "."), sep = ".")
      x <- as.magpie(unlist(x))
    }

    if (file %in% c("mean_lab", "stddevs_lab", "mean_col", "stddevs_col")) {
      files <- files[grep(".rds", files, ignore.case = TRUE)]
      x <- readRDS(file.path(folder, files))
      x <- as.magpie(x, data = 1)
    }

    if (file %in% c("inputs")) {
      files <- files[grep(".rds", files)]
      tmp <- readRDS(file.path(folder, files))
      x <- new.magpie(seq_along(tmp))
      getCells(x) <- tmp
    }
    return(x)
  }
}
