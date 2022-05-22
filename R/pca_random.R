#' Creates gene set for a universe of genes
#'
#' This function creates random gene sets 
#' @param Size Number of genes in a gene set
#' @param Universe Gene universe to sample random gene set from
#' @param n_perm Number of permutations to create a random gene set (default: 1000)
#' @return Random gene set
#' @export

create_gene_set <- function(universe, psize, n_perm = 1000) {
    gene_set <- lapply(1:n_perm, function(x) sample(universe, size = psize, replace = FALSE))
    return(gene_set)
}

#' Creates random gene sets and runs PCA analysis on network for a set of random gene sets
#'
#' This function creates random gene sets and runs PCA analysis on a set of random gene sets
#' @param RegNet Table of network with samples in columns
#' @param edges Table, containing information on "reg" and "tar" of RegNet
#' @param results_pca_pathways Output result table of pca_pathway function 
#' @param pathways_list A list of pathways
#' @param n_perm Number of permutations to create a random gene set (default: 1000)
#' @param ncores A number of cores to use (ncores: 1)
#' @return Dataframe with pca results for random gene sets
#' @export


pca_random <- function(RegNet, edges, results_pca_pathways, pathways_list, n_perm = 1000, ncores = 1) {
    pathways_size <- unique(results_pca_pathways$pathway_size)
    universe <- unique(unlist(pathways_list))
    res_pca_random <- list()
    for (m in 1:length(pathways_size)) {
    psize <- pathways_size[m]
    cat("Pathways with size", " ", psize, "\n")
    random_genes <- create_gene_set(universe, psize, n_perm = n_perm)
    res_pca <- pca_pathway(random_genes, RegNet, edges, ncores)
    res_pca_random[[m]] <- res_pca
    }
    res_pca_random_all <- as.data.frame(do.call("rbind", res_pca_random))
    return(res_pca_random_all)
    }