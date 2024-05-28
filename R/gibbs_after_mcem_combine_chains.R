#' Combining from all chains the posterior samples for parameters in the model and predicted gene expressions.
#'
#' @param tot_chain A numeric scalar. Total number of chains.
#' @param gibbs_after_mcem_algorithm_result A list of objects storing model constants. Should be the same as that input to the 'function gibbs_after_mcem_load_chains'.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#'
#' @return All saved posterior samples for parameters in the model and predicted gene expressions.
#' @export
gibbs_after_mcem_combine_chains<- function(tot_chain,
                                           gibbs_after_mcem_algorithm_result){

  # assign results to objects
  p =  gibbs_after_mcem_algorithm_result$p
  k =  gibbs_after_mcem_algorithm_result$k
  n =  gibbs_after_mcem_algorithm_result$n
  q =   gibbs_after_mcem_algorithm_result$q
  ind_x =  gibbs_after_mcem_algorithm_result$ind_x
  num_time_test =  gibbs_after_mcem_algorithm_result$num_time_test
  mc_num =  gibbs_after_mcem_algorithm_result$mc_num
  thin_step =  gibbs_after_mcem_algorithm_result$thin_step
  pathname =  gibbs_after_mcem_algorithm_result$pathname
  burnin = gibbs_after_mcem_algorithm_result$burnin
  pred_indicator = gibbs_after_mcem_algorithm_result$pred_indicator

  # remove finish assignment
  rm(gibbs_after_mcem_algorithm_result)

  num_sample<- ((mc_num - burnin)/thin_step) - 1

  big_a_final_array<- array(0,dim=c(p,k,(num_sample),tot_chain))

  big_z_final_array<- array(0,dim=c(p,k,(num_sample),tot_chain))

  beta_final_array<- array(0,dim=c(k,(num_sample),tot_chain))

  phi_final_array<- array(0,dim=c(p,(num_sample),tot_chain))

  pai_final_array<- array(0,dim=c(k,(num_sample),tot_chain))

  # number of time points stored in latent_y

  num_time_all<- (q+num_time_test)

  latent_y_final_array<- array(0,dim=c(num_time_all, k, n, (num_sample), tot_chain))

  if (ind_x){

    individual_mean_final_array<- array(0,dim=c(p,n,(num_sample),tot_chain))

    variance_g_final_array<- array(0,dim=c(p,(num_sample),tot_chain))

  }

  if (pred_indicator){

    pred_y_final_array<- array(0,dim=c(num_time_test,k,n,(num_sample),tot_chain))

    pred_x_final_array<- array(0,dim=c(num_time_test,p,n,(num_sample),tot_chain))

  }

  gibbs_after_mcem_load_chains_result<- NULL

  for (chain_index in 1:tot_chain){

    filename = paste0(pathname,"/","chain_", chain_index,"_result.RData")

    load(file = filename)

    # gibbs_after_mcem_load_chains_result<- get(filename) # get the object based on its name

    latent_y_final_array[,,,,chain_index]<- gibbs_after_mcem_load_chains_result$latent_y

    big_a_final_array[,,,chain_index]<- gibbs_after_mcem_load_chains_result$big_a

    big_z_final_array[,,,chain_index]<- gibbs_after_mcem_load_chains_result$big_z

    beta_final_array[,,chain_index]<- gibbs_after_mcem_load_chains_result$beta

    phi_final_array[,,chain_index]<- gibbs_after_mcem_load_chains_result$phi

    pai_final_array[,,chain_index]<- gibbs_after_mcem_load_chains_result$pai

    if (ind_x){

      individual_mean_final_array[,,,chain_index]<- gibbs_after_mcem_load_chains_result$individual_mean

      variance_g_final_array[,,chain_index]<- gibbs_after_mcem_load_chains_result$variance_g
    }

    if (pred_indicator){

      pred_x_final_array[,,,,chain_index]<- gibbs_after_mcem_load_chains_result$pred_x

      pred_y_final_array[,,,,chain_index]<- gibbs_after_mcem_load_chains_result$pred_y

    }

  }

  if(ind_x & pred_indicator){

    result_list<- list(pred_x = pred_x_final_array,
                       pred_y = pred_y_final_array,
                       latent_y = latent_y_final_array,
                       big_z = big_z_final_array,
                       big_a = big_a_final_array,
                       pai = pai_final_array,
                       phi = phi_final_array,
                       beta = beta_final_array,
                       individual_mean = individual_mean_final_array,
                       variance_g = variance_g_final_array,
                       num_sample = num_sample,
                       p = p,
                       k = k,
                       n = n,
                       q = q,
                       num_time_test = num_time_test,
                       tot_chain = tot_chain,
                       ind_x = ind_x,
                       pred_indicator = pred_indicator)

  } else if (!ind_x & pred_indicator){

    result_list<- list(pred_x = pred_x_final_array,
                       pred_y = pred_y_final_array,
                       latent_y = latent_y_final_array,
                       big_z = big_z_final_array,
                       big_a = big_a_final_array,
                       pai = pai_final_array,
                       phi = phi_final_array,
                       beta = beta_final_array,
                       num_sample = num_sample,
                       p = p,
                       k = k,
                       n = n,
                       q = q,
                       num_time_test = num_time_test,
                       tot_chain = tot_chain,
                       ind_x = ind_x,
                       pred_indicator = pred_indicator)

  } else if (ind_x & !pred_indicator){

    result_list<- list(latent_y = latent_y_final_array,
                       big_z = big_z_final_array,
                       big_a = big_a_final_array,
                       pai = pai_final_array,
                       phi = phi_final_array,
                       beta = beta_final_array,
                       individual_mean = individual_mean_final_array,
                       variance_g = variance_g_final_array,
                       num_sample = num_sample,
                       p = p,
                       k = k,
                       n = n,
                       q = q,
                       num_time_test = num_time_test,
                       tot_chain = tot_chain,
                       ind_x = ind_x,
                       pred_indicator = pred_indicator)

  } else {

    result_list<- list(latent_y = latent_y_final_array,
                       big_z = big_z_final_array,
                       big_a = big_a_final_array,
                       pai = pai_final_array,
                       phi = phi_final_array,
                       beta = beta_final_array,
                       num_sample = num_sample,
                       p = p,
                       k = k,
                       n = n,
                       q = q,
                       num_time_test = num_time_test,
                       tot_chain = tot_chain,
                       ind_x = ind_x,
                       pred_indicator = pred_indicator)

  }

  return(result_list)
}
