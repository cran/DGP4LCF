#' @title Plotting figures for factor score trajectory.
#'
#' @description This function is used to visualize results of factor score trajectories.
#'
#' @param factor_score_matrix A matrix of dimension (q, k, n), used to store results for factor scores.
#' @param factor_index A numeric scalar. Index of the factor of interest.
#' @param person_index A numeric scalar. Index of the person of interest.
#' @param trajectory_title A character. Title for the factor trajectory plot.
#' @param cex_main A numeric scalar. Text size of the title.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#'
#' @return  Trajectory of the designated person-factor.
#' @export
factor_score_trajectory<- function(factor_score_matrix, factor_index, person_index, trajectory_title, cex_main = 1){
  plot(factor_score_matrix[,factor_index, person_index],
       type ="l",
       ylab = "Value",
       xlab = "Time Index",
       main = trajectory_title,
       cex.main = cex_main)
}
