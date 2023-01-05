#' @title readGrassSoilEmu
#' @description Read files related to the training and optimization of the LPJml emulators.
#' @param subtype Subtype of file to be opened. Subtypes available:
#' 'weights', 'inputs',  'stddevs' and 'means'.
#' @return Magpie objects with a diverse inforamtion
#' @author Marcos Alves
#' @examples
#' \dontrun{
#' readSource("GrassSoilEmu",
#'   subtype =
#'     "ISIMIP3bv2:IPSL_CM6A_LR:ssp126:1965_2100:5f5fa2:weights", convert = F
#' )
#' }
#'
#' @import madrat
#' @importFrom stringr str_c
#'

readGrassSoilEmu <-
  function(subtype = "ISIMIP3bv2:IPSL_CM6A_LR:ssp126:1965_2100:5f5fa2:stddevs_lab") {
    subtype_split <- toolSplitSubtype(subtype, list(version = NULL, climatemodel = NULL, scenario = NULL, years = NULL, model = NULL, variable = NULL))
    file <- subtype_split$variable
    dirs <- list.dirs()
    folder <- grep(subtype_split$model, dirs, value = TRUE)
    if (length(dir.exists(file.path(folder))) != 0) {
      files_list <- list.files(folder)
      files <- files_list[grep(file, files_list)]

      if (file %in% "weights") {
        files <- files[grep(".rds", files, ignore.case = TRUE)]
        x <- readRDS(file.path(folder, files))
        x_dims <- lapply(x, dim)
        x_names <- NULL
        for (i in 1:length(x_dims)) {
          x_names[i] <- str_c(x_dims[[i]], collapse = "_")
        }
        names(x) <- paste(paste0("l", 1:length(x)), paste0(x_names, "."), sep = ".")
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
        x <- new.magpie(1:length(tmp))
        getCells(x) <- tmp
      }
      return(x)
    }
  }
