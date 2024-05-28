#' Generating different initials for multiple chains.
#'
#' @param ind_x A logical value. ind_x = TRUE uses the model including the intercept term for subject-gene mean in after-MCEM-Gibbs sampler; otherwise uses the model without the intercept term.
#' @param tot_chain A numeric scalar. Number of parallel chains.
#' @param mcem_parameter_setup_result A list of objects returned from the function 'mcem_parameter_setup'.
#' @param mcem_algorithm_result A list of objects returned from the function 'mcem_algorithm'.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#' vignette("bsfadgp_irregular_data_example", package = "DGP4LCF")
#'
#' @return Different initials for multiple chains.
#' @export
gibbs_after_mcem_diff_initials<- function(ind_x = TRUE,
                                          tot_chain = 5,
                                          mcem_parameter_setup_result,
                                          mcem_algorithm_result){

  #######################################################################################################################################################
  ################################################################# different initials for multiple chains ##############################################
  #######################################################################################################################################################
  # assign results to objects
  p<- mcem_parameter_setup_result$p
  k<- mcem_parameter_setup_result$k
  n<- mcem_parameter_setup_result$n
  q<- mcem_parameter_setup_result$q
  a_init = mcem_parameter_setup_result$big_a[[1]]
  z_init = mcem_parameter_setup_result$big_z[[1]]

  if (ind_x){
    gene_person_mean_est = mcem_parameter_setup_result$gene_person_mean_est
  }

  index_used = mcem_algorithm_result$index_used
  sigmay_record = mcem_algorithm_result$sigmay_record
  prior_sparsity<- mcem_algorithm_result$prior_sparsity

  # remove mcem_parameter_setup_result after assignment finished
  rm(mcem_parameter_setup_result,
     mcem_algorithm_result)


  # y_init: different for different chains
  sigma_y_init<- sigmay_record[index_used,,]

  y_init_multiple_chains<- array(0,dim=c(tot_chain, q, k, n))
  for (init_index in 1:tot_chain){
    for (n_index in 1:n){
      temp_vec<- mvtnorm::rmvnorm(1,mean=rep(0,times=(q*k)), sigma= sigma_y_init)
      y_init_multiple_chains[init_index,,,n_index]<- temp_vec
    }
  }

  a_init_multiple_chains<- array(0,dim=c(tot_chain, p, k))

  for (row_index in 1:p){
    for (col_index in 1:k){

      a_init_multiple_chains[, row_index, col_index]<- rnorm(tot_chain,
                                             mean = a_init[row_index, col_index],
                                             sd = 0.5)

    }
  }

  # z_init

  z_init_multiple_chains<- array(0,dim=c(tot_chain, p, k))

  for (init_index in 1:tot_chain){

    z_init_multiple_chains[init_index,,]<- z_init

  }

  # pai_init: different for different chains

  e0<- prior_sparsity*p
  f0<- (1-prior_sparsity)*p
  pai_init_multiple_chains<- rbeta(tot_chain, e0, f0)

  # phi_init: different for different chains
  phi_init_multiple_chains<- matrix(rgamma(p*tot_chain, 2, 0.05),
                    nrow=p,
                    ncol=tot_chain)

  # beta_init: may keep the same for different chains
  beta_init_multiple_chains<- rep(1, times = tot_chain)

  if(ind_x){
    gene_person_mean_init_multiple_chains<- array(0, dim = c(tot_chain, p, n))

    for (gene_index in 1:p){
      for (person_index in 1:n){

        #if ((gene_index == 14 && person_index == 12) || (gene_index == 14 && person_index == 69)){
        #  next
        #}

        gene_person_mean_init_multiple_chains[,gene_index,person_index]<- rnorm(tot_chain,
                                                                mean = gene_person_mean_est[gene_index, person_index],
                                                                sd = 0.5) # p*1 vector
      }
    }
  }

  #######################################################################################################################################################
  ################################################################# return results ######################################################################
  #######################################################################################################################################################
  if (!ind_x){
    result_list<- list(y_init_multiple_chains = y_init_multiple_chains,
                       a_init_multiple_chains = a_init_multiple_chains,
                       z_init_multiple_chains = z_init_multiple_chains,
                       pai_init_multiple_chains = pai_init_multiple_chains,
                       phi_init_multiple_chains = phi_init_multiple_chains,
                       beta_init_multiple_chains = beta_init_multiple_chains,
                       ind_x = ind_x)
  } else {
    result_list<- list(y_init_multiple_chains = y_init_multiple_chains,
                       a_init_multiple_chains = a_init_multiple_chains,
                       z_init_multiple_chains = z_init_multiple_chains,
                       pai_init_multiple_chains = pai_init_multiple_chains,
                       phi_init_multiple_chains = phi_init_multiple_chains,
                       beta_init_multiple_chains = beta_init_multiple_chains,
                       gene_person_mean_init_multiple_chains = gene_person_mean_init_multiple_chains,
                       ind_x = ind_x)
  }

  return(result_list)
}
