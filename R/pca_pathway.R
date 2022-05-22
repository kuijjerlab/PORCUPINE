#' Run PCA analysis for a list of pathways on network edges
#'
#' This function runs PCA analysis for a list of pathways
#' @param pathways_list list of pathways
#' @param RegNet Numeric matrix with samples in columns
#' @param edges Table, containing information on "reg" and "tar"
#' @param ncores A number of cores to use
#' @return Dataframe with pca results for pathways in a pathway file
#' @export

pca_pathway <- function(pathways_list, RegNet, edges, ncores) {
    res <- parallel::mclapply(pathways_list, function(pathway) {
    idx <- which(edges$tar %in% pathway)
    subnet <- RegNet[idx, ]
    pca_result <- run_pca(subnet)
    }, mc.cores = ncores)
    res <- as.data.frame(do.call("rbind", res))
    res$pathway <- rownames(res)
    res$pathway_size <- lengths(pathways_list)
    res <- res[,c(3,1,2,4)]
    return(res)
  }