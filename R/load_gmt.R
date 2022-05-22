#' Load a .gmt file
#'
#' This function loads a .gmt file.
#' @param gmt_file Path to the input .gmt file. Gmt file is downloaded from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
#' @return A list of pathways from the .gmt file
#' @export

load_gmt <- function(gmt_file) {
  pathways_list <- fgsea::gmtPathways(gmt_file)
  return(pathways_list)
}