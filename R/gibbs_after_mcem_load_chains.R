#' Loading the saved posterior samples for parameters in the model and predicted gene expressions.
#'
#' @param chain_index A numeric scalar. Index of the chain.
#' @param gibbs_after_mcem_algorithm_result A list of objects storing model constants.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example",  package = "DGP4LCF")
#'
#' @return All saved posterior samples for parameters in the model and predicted gene expressions.
#' @export
gibbs_after_mcem_load_chains<- function(chain_index,
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

  ###################################################################################################################################################
  ######################################## combine results from all chains ##########################################################################
  ###################################################################################################################################################

  ######################################## create empty containers ##################################################################################

  num_sample<- (mc_num - burnin)/thin_step - 1

  big_a_final_array<- array(0,dim=c(p,k,(num_sample)))

  big_z_final_array<- array(0,dim=c(p,k,(num_sample)))

  num_time_all<- (q+num_time_test)

  latent_y_final_array<- array(0,dim=c(num_time_all, k, n, (num_sample)))

  # need to combine for extra parameters if ind_x = 1
  if (ind_x){

    individual_mean_final_array<- array(0,dim=c(p,n,(num_sample)))

  }

  if (pred_indicator){

    pred_x_final_array<- array(0,dim=c(num_time_test,p,n,(num_sample)))

    pred_y_final_array<- array(0,dim=c(num_time_test,k,n,(num_sample)))
  }

  ######################################## load results from chains ##################################################################################

  directory_name = paste0(pathname, "/gibbs_after_mcem_chain_", chain_index)

  # fast read
  phi_unshaped_scan<- scan(file = paste0(directory_name, "/phi.csv"), sep = ",")
  # phi_unshaped<- matrix(phi_unshaped_scan[(burnin*(p+1) + 1):(num_sample*(p+1))], nrow=(num_sample-burnin), ncol= (p+1), byrow= T)
  phi_unshaped<- matrix(phi_unshaped_scan[(1):(num_sample*(p+1))], nrow=(num_sample), ncol= (p+1), byrow= T)
  phi_final_array<- t(phi_unshaped[,-(p+1)])

  # fast read
  pai_unshaped_scan<- scan(file =  paste0(directory_name,"/pai.csv"), sep = ",")
  # transform it into a matrix with desired format
  # pai_unshaped<- matrix(pai_unshaped_scan[(burnin*(k+1) + 1):(num_sample*(k+1))], nrow= (num_sample-burnin), ncol= (k+1), byrow=T) # num_sample*k
  pai_unshaped<- matrix(pai_unshaped_scan[(1):(num_sample*(k+1))], nrow= (num_sample), ncol= (k+1), byrow=T) # num_sample*k
  pai_final_array<- t(pai_unshaped[,-(k+1)])

  # fast read
  beta_unshaped_scan<- scan(file =  paste0(directory_name,"/beta.csv"), sep = ",")
  # beta_unshaped<- matrix(beta_unshaped_scan[(burnin*(k+1) + 1):(num_sample*(k+1))], nrow= (num_sample-burnin), ncol= (k+1), byrow=T) # num_sample*k
  beta_unshaped<- matrix(beta_unshaped_scan[(1):(num_sample*(k+1))], nrow= (num_sample), ncol= (k+1), byrow=T) # num_sample*k
  beta_final_array<- t(beta_unshaped[,-(k+1)])

  # fast_read
  big_z_unshaped_scan<- scan(file = paste0(directory_name,"/big_z.csv"), sep=",")
  # big_z_unshaped<- matrix(big_z_unshaped_scan[(burnin*(p*k+1) + 1):(num_sample*(p*k+1))], nrow = (num_sample-burnin), ncol = (p*k+1), byrow = T)
  big_z_unshaped<- matrix(big_z_unshaped_scan[(1):(num_sample*(p*k+1))], nrow = (num_sample), ncol = (p*k+1), byrow = T)

  big_a_unshaped_scan<- scan(file = paste0(directory_name,"/big_a.csv"), sep=",")
  # big_a_unshaped<- matrix(big_a_unshaped_scan[(burnin*(p*k+1) + 1):(num_sample*(p*k+1))], nrow = (num_sample-burnin), ncol = (p*k+1), byrow = T)
  big_a_unshaped<- matrix(big_a_unshaped_scan[(1):(num_sample*(p*k+1))], nrow = (num_sample), ncol = (p*k+1), byrow = T)

  latent_y_unshaped_scan <- scan(file = paste0(directory_name,"/latent_y.csv"), sep = ",")

  # latent_y_unshaped<- matrix(latent_y_unshaped_scan[(burnin*(num_time_all*k*n+1) + 1):(num_sample*(num_time_all*k*n+1))], nrow = (num_sample-burnin), ncol = (num_time_all*k*n+1), byrow = T)
  latent_y_unshaped<- matrix(latent_y_unshaped_scan[(1):(num_sample*(num_time_all*k*n+1))], nrow = (num_sample), ncol = (num_time_all*k*n+1), byrow = T)

  # fast read
  if (ind_x){
    variance_g_unshaped_scan<- scan(file = paste0(directory_name,"/variance_g.csv"), sep = ",")
    # variance_g_unshaped<- matrix(variance_g_unshaped_scan[(burnin*(p+1) + 1):(num_sample*(p+1))], nrow = (num_sample-burnin), ncol= (p+1), byrow= T)
    variance_g_unshaped<- matrix(variance_g_unshaped_scan[(1):(num_sample*(p+1))], nrow = (num_sample), ncol= (p+1), byrow= T)
    variance_g_final_array<- t(variance_g_unshaped[,-(p+1)])

    individual_mean_unshaped_scan<-  scan(file = paste0(directory_name,"/individual_mean.csv"), sep=",")
    # individual_mean_unshaped<- matrix(individual_mean_unshaped_scan[(burnin*(p*n+1) + 1):(num_sample*(p*n+1))], nrow = (num_sample-burnin), ncol = (p*n+1), byrow = T)
    individual_mean_unshaped<- matrix(individual_mean_unshaped_scan[(1):(num_sample*(p*n+1))], nrow = (num_sample), ncol = (p*n+1), byrow = T)

  }

  if (pred_indicator){

    pred_y_unshaped_scan <- scan(file =  paste0(directory_name,"/pred_y.csv"), sep = ",")
    # pred_y_unshaped<- matrix(pred_y_unshaped_scan[(burnin*(num_time_test*k*n+1) + 1):(num_sample*(num_time_test*k*n+1))], nrow = (num_sample-burnin), ncol = (num_time_test*k*n+1), byrow = T)
    pred_y_unshaped<- matrix(pred_y_unshaped_scan[(1):(num_sample*(num_time_test*k*n+1))], nrow = (num_sample), ncol = (num_time_test*k*n+1), byrow = T)

    pred_x_unshaped_scan <- scan(file =  paste0(directory_name,"/pred_x.csv"), sep = ",")
    # pred_x_unshaped<- matrix(pred_x_unshaped_scan[(burnin*(num_time_test*p*n+1) + 1):(num_sample*(num_time_test*p*n+1))], nrow = (num_sample-burnin), ncol = (num_time_test*p*n+1), byrow = T)
    pred_x_unshaped<- matrix(pred_x_unshaped_scan[(1):(num_sample*(num_time_test*p*n+1))], nrow = (num_sample), ncol = (num_time_test*p*n+1), byrow = T)
  }

  message(paste0("this is to save the result for which chain:",chain_index))

  for(i in 1:(num_sample)){

    big_z_final_array[,,i] <- matrix(as.numeric(big_z_unshaped[i, -(p*k+1)]), nrow=p, ncol=k,byrow=T)

    big_a_final_array[,,i]<- matrix(as.numeric(big_a_unshaped[i, -(p*k+1)]), nrow=p, ncol=k,byrow=T)

    latent_y_final_array[,,,i]<- array(as.numeric(latent_y_unshaped[i, -(num_time_all*k*n+1)]), dim = c(num_time_all,k,n)) # the flatting process has no problem

    if (ind_x){
      individual_mean_final_array[,,i]<- matrix(as.numeric(individual_mean_unshaped[i, -(p*n+1)]), nrow=p, ncol=n,byrow=T)
    }

    if (pred_indicator){

      pred_y_final_array[,,,i]<- array(as.numeric(pred_y_unshaped[i, -(num_time_test*k*n+1)]), dim = c(num_time_test,k,n))

      pred_x_final_array[,,,i]<- array(as.numeric(pred_x_unshaped[i, -(num_time_test*p*n+1)]), dim = c(num_time_test,p,n))

    }

  }

  ###################################################################################################################################################
  ##################################################### return results ##############################################################################
  ###################################################################################################################################################

  # combined samples from all chains
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
                       num_sample = num_sample)

  } else if (!ind_x & pred_indicator){
    result_list<- list(pred_x = pred_x_final_array,
                       pred_y = pred_y_final_array,
                       latent_y = latent_y_final_array,
                       big_z = big_z_final_array,
                       big_a = big_a_final_array,
                       pai = pai_final_array,
                       phi = phi_final_array,
                       beta = beta_final_array,
                       num_sample = num_sample)

  } else if (ind_x & !pred_indicator){

    result_list<- list(latent_y = latent_y_final_array,
                       big_z = big_z_final_array,
                       big_a = big_a_final_array,
                       pai = pai_final_array,
                       phi = phi_final_array,
                       beta = beta_final_array,
                       individual_mean = individual_mean_final_array,
                       variance_g = variance_g_final_array,
                       num_sample = num_sample)

  } else {

    result_list<- list(latent_y = latent_y_final_array,
                       big_z = big_z_final_array,
                       big_a = big_a_final_array,
                       pai = pai_final_array,
                       phi = phi_final_array,
                       beta = beta_final_array,
                       num_sample = num_sample)

  }

  return(result_list)

}
