#' @title readDrainage
#' @description Read drainage file from LPJmL inputs
#'
#' @param subtype Drainage file to be read: drainagestn.bin (default) or drainage.bin
#'
#' @return Magpie object with river structure information
#' @author Jens Heinke, Felicitas Beier
#'
#' @examples
#'
#' \dontrun{
#'   readSource("Drainage")
#' }
#'

readDrainage <- function(subtype="drainagestn"){

  ### Read in drainage
  # Number of cells
  NCELL    <- 67420
  NCRUCELL <- 67420
  # Read binary file
  zz <- file(paste0("/",subtype,".bin"))
  zz <- file("C:/Users/beier/Documents/doktorarbeit/MAgPIE_Water/River_Routing_Postprocessing/drainagestn.bin","rb")
  seek(zz,where=43,origin="start")
  x <- readBin(zz, integer(), n=2*NCRUCELL, size=4)
  close(zz)

  ### River structure
  nextcell <- x[c(1:NCRUCELL)*2-1]
  dist <- x[c(1:NCRUCELL)*2]

  nextcell[which(nextcell<0)] <- -1
  nextcell[which(nextcell>=0)] <- nextcell[which(nextcell>=0)] + 1

  # Determine downstream cell list
  dummy <- array(data=-9999,dim=c(NCRUCELL))
  c <- 1
  i <- 1
  dummy[i] <- nextcell[c]
  while(dummy[i]>0){
    i <- i + 1
    dummy[i] <- nextcell[dummy[i-1]]
  }
  dsclist <- list(dummy[0:(i-1)])

  for(c in 2:NCRUCELL){
    dummy[] <- -9999
    i <- 1
    dummy[i] <- nextcell[c]
    while(dummy[i]>0){
      i <- i + 1
      dummy[i] <- nextcell[dummy[i-1]]
    }
    dsclist[[length(dsclist)+1]] <- dummy[0:(i-1)]
  }

  # Determine endcell and dictance to endcell
  cellstoend <- array(data=0,dim=c(NCRUCELL))
  endcell <- array(data=0,dim=c(NCRUCELL))
  for(i in 1:NCRUCELL){
    endcell[i] <- i
    while(nextcell[endcell[i]]>0){
      endcell[i] <- nextcell[endcell[i]]
      cellstoend[i] <- cellstoend[i] + 1
    }
  }

  # Determine calcorder
  basinids <- unique(endcell)
  calcorder <- array(data=0,dim=c(NCRUCELL))
  for(b in 1:length(basinids)){
    basincells <- which(endcell==basinids[b])
    calcorder[basincells] <- (cellstoend[basincells] - max(cellstoend[basincells]) - 1)*(-1)
  }

  # Determine upstream cell list
  usclist <- list()
  for(c in 1:NCRUCELL){
    if(c%%1000 == 0) print(c)
    basincells <- which(endcell==endcell[c])
    dummy <- numeric()
    for(cell in basincells){
      if(is.element(c,dsclist[[cell]])) dummy <- c(dummy,cell)
    }
    usclist[[c]] <- dummy
  }




  return(x)
}


out <- new.magpie()


  library(magclass)
  cells <- c("1.1","1.2","1.3","2.1","2.2","3.1","3.4","3.10","4.1","4.2")
  out <- new.magpie(cells,fill=1)
  getSets(out,fulldim=FALSE)[1] <- "cell.upstreamcell"
  as.integer(getItems(mselect(out,cell="1"),"upstreamcell"))






  data <- list(list(nextcell=nextcell),list(downstreamcells=dsclist),list(upstreamcells=usclist),list(endcell=endcell),list(calcorder=calcorder))

  save(data,file="C:/Users/beier/Documents/doktorarbeit/MAgPIE_Water/River_Routing_Postprocessing/River_structure_stn.rda")





  mf <- getConfig("mappingfolder")
  if(is.null(mf)) stop('No mappingfolder specified in used cfg! Please load a config with the corresponding information!')
  fname <- paste0(mf,"/",type,"/",name)
  if(error.existing && file.exists(fname)) {
    stop('Mapping "',name,'" exists already!')
  }
  fname <- gsub("/+","/",fname)

  if(is.magpie(map)) {
    pattern <- "^(.*)\\.(.*)\\.(.*)\\.(.*)$"
    map <- data.frame(cell    = sub(pattern, "\\1.\\3", getCells(map)),
                      cluster = sub(pattern, "\\2.\\4", getCells(map)),
                      region  = sub(pattern, "\\2", getCells(map)),
                      country = sub(pattern, "\\1", getCells(map)),
                      global  = "GLO")
  } else stop("Cannot handle this mapping format!")

  filetype <- tolower(file_ext(fname))
  if(filetype=="csv") {
    write.table(map, fname, sep=";", quote=FALSE)
  } else {
    stop("Unsupported filetype \"", filetype,"\"")
  }
}
