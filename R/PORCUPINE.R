#' Determine statistical significance of PCA score of a pathway
#'
#' Compares the observed PCA score for a pathway to a set of PCA scores of random gene sets of the same size as pathway. Calculates p-value and effect size.
#' @param res_pca_pathway Result table of pca_pathway function for one pathway
#' @param res_pca_random_subset Result table of pca_random function for one pathway 
#' @return Table of statistics (pvalue and effect size) for each pathway
#' @export

calculate_statistics <- function(res_pca_pathway, res_pca_rndm_set, test = t.test) {
    pc1_rndm <- res_pca_rndm_set$pc1
    pathway <- res_pca_pathway[["pathway"]]
    path_size <- res_pca_pathway[["pathway_size"]]
    pc1_pathway <- as.numeric(res_pca_pathway[["pc1"]])
    pvalue <- test(pc1_rndm, mu = pc1_pathway, alternative = "less")$p.value
    es <- lsr::cohensD(pc1_rndm, mu = pc1_pathway)
    stat_res <- data.frame("pathway" = pathway,
                "pathway_size" = path_size, "pval" = pvalue,
                "es" = es)
    return(stat_res)
}

#' Determine statistical significance of PCA scores of pathways
#'
#' This function compares results of PCA analyses for each pathway versus set of permutated gene sets. Calculates p-value and effect size for each pathway.
#' @param res_pca_pathways Output result table of pca_pathway function
#' @param res_pca_random Output result table of pca_random function
#' @return Table of statistics (pvalue and effect size) for each pathway
#' @export

PORCUPINE <- function(res_pca_pathways, res_pca_rndm, test=t.test) {
  pathways_size <- unique(res_pca_pathways$pathway_size)
  res_all <- list()
  for (k in 1:length(pathways_size)) {
    path_size <- pathways_size[k]
    print(path_size)
    # select pca results for the pathways of the size path_size
    res_pathway_set <- res_pca_pathways %>%
                          dplyr::filter(pathway_size == path_size)
    # select pca results for random sets of the size of a pathway
    res_rndm_set <- res_pca_rndm %>%
                         dplyr::filter(pathway_size == path_size)
    # calculate pvalue and effect size for each pathway versus random gene sets 
    res <- apply(res_pathway_set, 1, function(x) calculate_statistics(x, res_rndm_set, test=t.test))
    names(res) <- NULL
    res_all <- c(res, res_all)
}
    res_all <- do.call("rbind", res_all)
    return(res_all)
}
