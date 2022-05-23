#' Run PCA analysis on network data
#'
#' This function runs PCA analysis
#' @param data Numeric matrix with samples in columns
#' @return Dataframe with pca results
#' @export

run_pca <- function(data) {
    data_t <- Matrix::t(data)
    # Perform scaling and centering of the data
    res_pca <- stats::prcomp(data_t, scale. = T, center = T)
    # Extract variance explained by the first PC
    pc1 <- summary(res_pca)$importance[2, 1] * 100
    pca_result <- data.frame("pc1" = pc1, "n_edges" = nrow(data))
    return(pca_result)
    }
