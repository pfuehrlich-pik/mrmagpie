#' Tool: StoreMapping
#'
#' Function which creates a mapping file
#'
#'
#' @param map mapping data to be stored.
#' @param name File name of the mapping file. Supported file types are currently csv (, or ; separated)
#' and rds. Use code{\link{toolConvertMapping}} to convert between both formats
#' @param type [only relevant when written to mapping folder, otherwise ignored!] Mapping type
#' (e.g. "regional", "cell", or "sectoral"). Can be set to NULL if file is not stored in a type
#' specific subfolder
#' @param where vector of locations the file should be written to. Available options are "mappingfolder"
#' (madrat mapping folder), "output" (madrat output folder) and "local" (current directory)
#' @param error.existing Boolean which decides whether an error is returned if
#' the mapping file does not exist or not.
#' @author Jan Philipp Dietrich
#' @seealso \code{\link{calcOutput}}, \code{\link{toolConvertMapping}}
#' @importFrom tools file_ext
#' @importFrom utils write.table
#' @importFrom magclass is.magpie
#' @importFrom madrat getConfig
#' @export
#'
toolStoreMapping <- function(map, name, type=NULL, where="mappingfolder", error.existing=TRUE) {
  fnames <- NULL
  if("mappingfolder" %in% where) {
    mf <- getConfig("mappingfolder")
    if(is.null(mf)) stop('No mappingfolder specified in used cfg! Please load a config with the corresponding information!')
    fname <- paste0(mf,"/",type,"/",name)
    if(error.existing && file.exists(fname)) {
      stop('Mapping "',name,'" exists already in mapping folder!')
    }
    fnames <- fname
  }
  if("outputfolder" %in% where) {
    mf <- getConfig("outputfolder")
    if(is.null(mf)) stop('No output folder specified in used cfg! Please load a config with the corresponding information!')
    fname <- paste0(mf,"/",name)
    if(error.existing && file.exists(fname)) {
      stop('Mapping "',name,'" exists already in output folder!')
    }
    fnames <- c(fnames,fname)
  }
  if("local" %in% where) {
    fname <- name
    if(error.existing && file.exists(fname)) {
      stop('Mapping "',name,'" exists already in local folder!')
    }
    fnames <- c(fnames,fname)
  }

  fnames <- gsub("/+","/",fnames)

  if(is.magpie(map)) {
    pattern <- "^(.*)\\.(.*)\\.(.*)\\.(.*)$"
    map <- data.frame(cell    = sub(pattern, "\\1.\\3", getCells(map)),
                      cluster = sub(pattern, "\\2.\\4", getCells(map)),
                      region  = sub(pattern, "\\2", getCells(map)),
                      country = sub(pattern, "\\1", getCells(map)),
                      global  = "GLO")
  } else stop("Cannot handle this mapping format!")

  for(fname in fnames) {
    filetype <- tolower(file_ext(fname))
    if(filetype=="csv") {
      write.table(map, fname, sep=";", quote=FALSE)
    } else if(filetype=="rds") {
      saveRDS(map, fname, compress = "xz")
    } else {
      stop("Unsupported filetype \"", filetype,"\"")
    }
  }
}
