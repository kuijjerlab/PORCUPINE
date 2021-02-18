#' Creates random gene sets and runs PCA analysis on a set of random gene sets
#'
#' This function creates random gene sets and runs PCA analysis on a set of random gene sets
#' @param net A numeric matrix or a data frame, representing network
#' @param res.ptws Output result table of pca_gmt function 
#' @param pt.filt A list of filtered pathways 
#' @param n.perm Number of permutations to create random gene set
#' @param ncores A number of cores to use
#' @return Dataframe with pca results for random gene sets
#' @export

pca_random <- function(net, res.ptws, pt.filt, n.perm, ncores) {
  sizes.ptw <- unique(res.ptws$ptw_size)
  gene.set <- unique(unlist(pt.filt))
  res.random.all <-NULL
  res.random <- NULL
  for (m in 1:length(sizes.ptw)) {
    ptw.size <- sizes.ptw[m] 
    cat("Pathways with size"," ", ptw.size, "\n")
    random_gene_set <- lapply(1:n.perm, function(z) sample(gene.set, size=ptw.size, replace=FALSE))
    res.random <- parallel::mclapply(random_gene_set, function(ptw) {
      edges <- net[net$tar %in% ptw,] 
      edges.t <- t(edges[,-c(1:3)])
      res.pca.comp <- prcomp(edges.t, scale. =T, center=T) ### perfrom scaling
      ntfs <- length(unique(net$reg))
      ptw.size <- nrow(edges)/ntfs
      pca.var <- c("ptw_size"=ptw.size, "pc1"=summary(res.pca.comp)$importance[2,1]*100,"n_edges"=nrow(edges))
      pca.var
    },  mc.cores = ncores)
    res.random.all <- c(res.random, res.random.all)
  }
  res.random.all <- as.data.frame(do.call("rbind", res.random.all))
  res.random.all
}
