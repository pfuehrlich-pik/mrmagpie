#' @title readGrassYldEmu
#' @description Read files related to the training and optimization of the LPJml emulators.
#' @param subtype Subtype of file to be opened. Subtypes available:
#' 'max_harvest', 'weights', 'inputs',  'stddevs' and 'means'.
#' @return Magpie objects with a diverse inforamtion
#' @author Marcos Alves
#' @examples
#'
#' \dontrun{
#' readSource("GrassYldEmu", subtype = "GrassYldEmu:20f33a2280.weights", convert="onlycorrect")
#' }
#'
#' @import madrat
#'
#'

readGrassYldEmu <-
  function(subtype = "109325f71e.inputs") {
    if (grepl("\\.", subtype)) {
      type     <- strsplit(subtype, split = "\\.")
      folder      <- unlist(type)[1]
      subtype     <- unlist(type)[2]

    } else {

      stop("Subtype did not define either: RCP, climate model or file. Please check the defintion of your subtype parameter.")

    }

    if (dir.exists(file.path(folder))) {
      files_list <- list.files(folder)
      file <- files_list[grep(subtype, files_list)]
     } else {
       stop(paste("Path", folder, "does not exist. Check the defition of your",
                  "subtype or the folder structure you are trying to access."))
     }

    if (subtype == "max_harvest") {
      max_grass <- unlist(readRDS(file.path(folder,file)))
      nyears <- length(max_grass)/59199
      matrix <- matrix(max_grass, ncol = nyears)
      magpie <- as.magpie(matrix, spatial = 1, temporal = 2)

      x <-  collapseNames(magpie)
      getNames(x) <- subtype
    }

    if(subtype %in% c("weights", "inputs", "mean", "stddevs")){
      x <- readRDS(file.path(folder,file))
      x <- as.magpie(as.matrix(x), spatial = 1)
      getNames(x) <- subtype
    }

    return(x)

  }

