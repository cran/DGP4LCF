#' @title Visualizing cross-correlations among factors.
#'
#' @param k A numeric scalar. Number of latent factors.
#' @param q A numeric scalar. Number of time points in the covariance matrix of factors.
#' @param cov_input A matrix of dimension (kq, kq). The covariance matrix of the vector obtained from vectorizing the matrix of latent factor scores.
#' @param title A character. Title for the plot.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#' vignette("bsfadgp_irregular_data_example", package = "DGP4LCF")
#'
#' @return Visualization of cross-correlations among factors.
#' @export
mcem_cov_plot<- function(k, q, cov_input, title){
  col <- grDevices::colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  cov_matrix<- diag(1,k)
  for (row_index in 1:(k-1)){
    for (col_index in ((row_index + 1):k)){
      cov_matrix[row_index, col_index]<- cov2cor(cov_input)[(row_index-1)*q+1, (col_index-1)*q+1]
      cov_matrix[col_index, row_index]<- cov_matrix[row_index, col_index]
    }
  }

  corrplot::corrplot(cov_matrix,
                     title = title,
                     type="upper",
                     method="square",
                     addCoef.col = "black",
                     tl.col="black",
                     col = col (200),
                     mar=c(0,0,1,0),
                     diag=FALSE)
}
