#' Filter a pathway list
#'
#' This function filters a list of pathways.
#' @param pathways_list of pathways
#' @param edges Table, containing information on "reg" and "tar"
#' @param minSize Minimum size for number of genes in a pathway (default: 5)
#' @param maxSize Maximum size for a number of genes in a pathway (default: 200)
#' @return A list of filtered pathways
#' @export

filter_pathways <- function(pathways_list, edges, minSize = 5, maxSize = 200) {
  pathways_filt <- plyr::ldply(pathways_list, data.frame) %>%
                   dplyr::rename(pathway = ".id", gene = "X..i..") %>%
                   dplyr::filter(gene %in% edges$tar) %>%
                   with(., split(gene, pathway))
  pathways_filt <- purrr::keep(pathways_filt, function(x) length(x) >= minSize)
  pathways_filt <- purrr::keep(pathways_filt, function(x) length(x) <= maxSize)
  return(pathways_filt)
}