#' Compares PCA results for each pathway versus PCA results of a set set of random gene sets of the same size, calculates p-value and effect size for each pathway
#'
#' This function compares results of PCA analyses for each pathway versus set of permutated pathways of the same size, calculates p-value and effect size for each pathway
#' @param res.ptws Output result table of pca_gmt function 
#' @param res.random Output result table of pca_random function 
#' @export

PORCUPINE <- function(res.ptws, res.random) {
  sizes.ptw <- unique(res.ptws$ptw_size)
  
  pval.dat.all <- NULL
  for (k in 1:length(sizes.ptw)) {
    ptw.size <- sizes.ptw[k]
    ptw.block <- res.ptws[res.ptws$ptw_size %in% ptw.size,]
    ptw.perm <- res.random[res.random$ptw_size %in% ptw.size,]
    
    pval.dat <- NULL
    pval.dat.block  <- NULL
    
    for (l in 1:length(rownames(ptw.block))) {
      pev1 <- ptw.block$pc1[l]
      pc1.pval <- t.test(ptw.perm$pc1,mu=pev1, alternative="less")$p.value
      pc1.es <- lsr::cohensD(ptw.perm$pc1, mu = pev1)
      pval.dat <- data.table::data.table("pathway"=rownames(ptw.block)[l],"size"=ptw.block$sptw_ize[l], "pc1.pval"=pc1.pval, "pc1.es"=pc1.es)
      pval.dat.block <- rbind(pval.dat.block, pval.dat)
    }
    
    pval.dat.all <-  rbind( pval.dat.all, pval.dat.block)
  }
  pval.dat.all
}
