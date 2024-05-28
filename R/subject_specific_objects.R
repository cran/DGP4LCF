#' Constructing subject-specific objects required for Gibbs sampler (for subjects with incomplete observations only).
#'
#' @details This function is used to extract subject-specific factor covariance matrix from the complete factor covariance matrix, through constructing subject-specific indicator matrix, which indicates time indexes when gene expression are available.
#'
#' @param k A numeric scalar. Number of latent factors.
#' @param q A numeric scalar. Number of time points in the complete factor covariance matrix.
#' @param a_full A q-dimensional numeric vector. Complete time sorted from early to late.
#' @param a_avail A vector of time when gene expressions are available, sorted from early to late.
#' @param cor_all A matrix of dimension (kq, kq). Correlation matrix of latent factor scores.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#' vignette("bsfadgp_irregular_data_example", package = "DGP4LCF")
#'
#' @return Subject-specific objects needed for Gibbs sampler.
#' @export
subject_specific_objects<- function(k, q, a_full, a_avail, cor_all){

    # a full identity matrix
    full_matrix<- diag(1, nrow = (k*q), ncol = (k*q))

    # find indexes of missing/observed time (in the original a_full)

    missing_time_index<- which(a_full %in% setdiff(a_full, a_avail))

    missing_time_num<- length(a_full) - length(a_avail)

    missing_row_index<- NULL

    for (index in 1:missing_time_num){

      index_temp <- missing_time_index[index]

      for (factor_index in 1:k){
        missing_row_index<- c(missing_row_index, (factor_index-1)*q + index_temp)
      }

    }

    # observed row index: note that setdiff find the difference in the 1st argument compared to the 2nd argument
    obs_row_index<- setdiff(1:(k*q), missing_row_index)

    # remove missing rows
    obs_matrix<- full_matrix[-missing_row_index,]
    missing_matrix<- full_matrix[-obs_row_index,]

    ######################## compute prod_covnew_covestinv and cov_cond_dist for each individual #########################

    # n_i is the number of time points in the training data and m_i is the number of time points in the test data
    cov_est<-   obs_matrix%*%cor_all%*%t(obs_matrix) # kn_i * kn_i
    cov_est_inv<- solve(cov_est) # kn_i * kn_i

    cov_new<-    missing_matrix%*%cor_all%*%t(missing_matrix) # km_i * km_i

    cov_new_est<- missing_matrix%*%cor_all%*%t(obs_matrix) # km_i * kn_i

    # mean of the conditional distribution
    prod_covnew_covestinv<- cov_new_est%*%cov_est_inv # km_i * kn_i

    # covariance of the conditional distribution
    cov_cond_dist<-  cov_new- (prod_covnew_covestinv%*%t(cov_new_est)) # km_i * km_i


    # return results
    result<- list(obs_matrix = obs_matrix,
                  missing_matrix = missing_matrix,
                  missing_time_num = missing_time_num,
                  missing_time_index = missing_time_index,
                  prod_covnew_covestinv = prod_covnew_covestinv,
                  cov_cond_dist = cov_cond_dist,
                  cov_est_inv =  cov_est_inv)

}
