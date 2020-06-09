#' @title calcCluster
#' @description This function calculates the aggregation mapping for a given cluster methodology
#' @param ctype aggregation clustering type, which is a combination of a single letter, indicating the cluster methodology, and a number,
#' indicating the number of resulting clusters. Available methodologies are hierarchical clustering (h), normalized k-means clustering
#' (n) and combined hierarchical/normalized k-means clustering (c). In the latter hierarchical clustering is used to determine the
#' cluster distribution among regions whereas normalized k-means is used for the clustering within a region.
#' @param regionscode regionscode of the regional mapping to be used. Must agree with the regionscode of the mapping
#' mentioned in the madrat config! Can be retrieved via \code{regionscode()}.
#' @param seed Seed for Random Number Generation. If set to NULL it is chosen automatically, if set to an integer it will
#' always return the same pseudo-random numbers (useful to get identical clusters under identical inputs for n and c
#' clustering)
#' @param weight Should specific regions be resolved with more or less detail? Values > 1 mean higher share, < 1 lower share
#' e.g. cfg$cluster_weight <- c(LAM=2) means that a higher level of detail for region LAM if set to NULL all weights will be
#' assumed to be 1 (examples: c(LAM=1.5,SSA=1.5,OAS=1.5), c(LAM=2,SSA=2,OAS=2))
#' @return magpie object in cellular resolution
#' @author Jan Philipp Dietrich
#'
#' @examples
#' \dontrun{ calcOutput("Cluster", type="c200", aggregate = FALSE) }
#' @importFrom madrat calcOutput

calcCluster <- function(ctype, regionscode=madrat::regionscode(), seed=42, weight=NULL){

  mode <- substr(ctype,0,1)
  ncluster <- as.integer(substring(ctype,2))

  if(mode=="n") {
    mapping <- calcOutput("ClusterKMeans", regionscode=regionscode, ncluster=ncluster, weight=weight,
                          seed=seed, aggregate=FALSE)
  } else if(mode=="h" | mode=="w" | mode=="s") {
    mapping <- calcOutput("ClusterHierarchical", regionscode=regionscode, ncluster=ncluster,
                          mode=mode, weight=weight, aggregate=FALSE)
  } else if(mode=="c"){
    calcCPR <- function(x) {
      clusters <- table(sub("\\..*$","",unique(sub("\\..*\\.",".",x))))
      cells <- table(sub("\\..*$","",x))
      return(cbind(cells,clusters))
    }
    tmpmap  <- calcOutput("ClusterHierarchical", regionscode=regionscode, ncluster=ncluster,
                          mode="h", weight=weight, aggregate=FALSE)
    mapping <- calcOutput("ClusterKMeans", regionscode=regionscode, ncluster=ncluster,
                          weight=weight, cpr=calcCPR(sub("^[^\\.]*\\.","",getCells(tmpmap))),
                          seed=seed, aggregate=FALSE)
  } else {
    stop("Unkown clustering mode ",mode,"!")
  }
  #wkey <- ifelse(is.null(weight), "", gsub(".","",paste0("_",names(weight),weight,collapse=""),fixed=TRUE))

  # !!! HOW TO FORWARD NAME INFORMATION? !!! #
  #write.spam(spam,file.path(ofolder,paste(hr,"-to-",lr,wkey,"_sum.spam",sep="")))
  #saveRDS(spam2mapping(spam,rownames(cdata)), file.path(ofolder,paste(hr,"-to-",lr,wkey,"_mapping.rds",sep="")))

  return(list(
    x=mapping,
    weight=NULL,
    unit="1",
    description="Mapping between cells and cluster",
    isocountries=FALSE))
}
