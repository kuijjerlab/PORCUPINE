#' Filter a pathway file 
#'
#' This function filters a list of pathways.
#' @param pt.list List of pathways
#' @param net A numeric matrix or a data frame, representing network
#' @param length.thres Threshold to filter pathways, based on their length
#' @return A list of filtered pathways
#' @export

filter_pt <- function(pt.list, net, length.thres) {
  len <- lengths(pt.list)
  len.filt <- len[len<=length.thres] 
  pt.filt <- pt.list[names(pt.list) %in% names(len.filt)]
  df <- plyr::ldply(pt.filt, data.frame)
  colnames(df) <- c("ptw", "gene")
  df <- df[df$gene %in% net$tar,]
  pt.filt.net <- split(as.character(df$gene),df$ptw)
  pt.filt.net
}