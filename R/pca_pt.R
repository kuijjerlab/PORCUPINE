#' Run PCA analysis for a list of pathways on network edges
#'
#' This function runs PCA analysis for a list of pathways
#' @param pt.filt A list of filtered pathways 
#' @param net A numeric matrix or a data frame, representing network
#' @param ncores A number of cores to use
#' @return Dataframe with pca results for pathways in pathway file
#' @export

pca_pt <- function(pt.filt, net, ncores) {
  res <- parallel::mclapply(pt.filt, function(ptw) {
    edges <- net[net$tar %in% ptw,] 
    edges.t <- t(edges[,-c(1:3)])
    res.pca.comp <- prcomp(edges.t, scale. =T, center=T) ### perfrom scaling
    ntfs <- length(unique(net$reg))
    ptw.size <- nrow(edges)/ntfs
    pca.var <- c("ptw_size"=ptw.size, "pc1"=summary(res.pca.comp)$importance[2,1]*100,"n_edges"=nrow(edges))
    pca.var
  }, mc.cores =ncores)
  res <- as.data.frame(do.call("rbind", res))
}