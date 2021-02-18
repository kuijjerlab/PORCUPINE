#' Load a gmt file
#'
#' This function loads a gmt file.
#' @param gmt.file Path to the input gmt file
#' @return A list of the gmt.file
#' @export

load_gmt <- function(gmt.file) {
  lines = readLines(gmt.file)
  lines = strsplit(lines,split="\t")
  pt = list()
  for (i in seq(lines)) {
    pt[[lines[[i]][1]]] = lines[[i]][c(3:length(lines[[i]]))]
  }
  pt
}