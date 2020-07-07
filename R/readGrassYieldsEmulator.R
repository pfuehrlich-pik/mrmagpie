#' @title readGrassYldEmu
#' @description Read the results of the grass harvest maximization performed for the grass havest emulators
#' @param subtype subtype of file to be opened
#' @return Magpie object on a celullar level with maximum gass harvest in gC/m²
#' @author Marcos Alves
#' @examples
#'
#' \dontrun{
#' readSource("GrassYldEmu", subtype = "GrassYldEmu:20f33a2280.max_harvest", convert="onlycorrect")
#' }
#'
#' @import madrat
#'
#'

readGrassYldEmu <-
  function(subtype = "GrassYldEmu:20f33a2280.max_harvest") {
    if (grepl("\\.", subtype) & grepl("\\:", subtype)) {
      type     <- strsplit(gsub(":", "/" , subtype), split = "\\.")
      folder      <- unlist(type)[1]
      subtype     <- unlist(type)[2]

    } else {

      stop("Subtype did not define either: RCP, climate model or file. Please check the defintion of your subtype parameter.")

    }

    if (dir.exists(path(folder))) {
      files_list <- list.files(folder)
      file <- files_list[grep(subtype, files_list)]
     } else {
       stop(paste("Path", folder, "does not exist. Check the defition of your",
                  "subtype or the folder structure you are trying to access."))
     }

    if (subtype == "max_harvest") {
      max_grass <- unlist(readRDS(path(folder,file)))
      nyears <- length(max_grass)/59199
      matrix <- matrix(max_grass, ncol = nyears)
      magpie <- as.magpie(matrix, spatial = 1, temporal = 2)

      x <-  collapseNames(magpie)
      getNames(x) <- subtype
    }

    return(x)

  }

