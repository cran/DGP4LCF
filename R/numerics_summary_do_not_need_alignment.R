#' Numerical summary for important continuous variables that do not need alignment.
#'
#' @details This function corresponds to Algorithm 2: Steps 3 and 4 in the main manuscript; therefore reader can consult the paper for more explanations.
#'
#' @param burnin A numeric scalar. The saved samples are already after burnin; therefore the default value for this parameter here is 0. Can discard further samples if needed.
#' @param thin_step A numeric scalar. The saved samples are already after thinning; therefore the default value for this parameter here is 1. Can be further thinned if needed.
#' @param pred_x_truth_indicator A logical value. pred_x_truth_indicator = TRUE means that truth of predicted gene expressions are available. The default value is FALSE.
#' @param pred_x_truth Only needed if pred_x_truth_inidcator = TRUE. An array of dimension (n, p, num_time_test), storing true gene expressions in the testing data.
#' @param gibbs_after_mcem_combine_chains_result A list of objects returned from the function 'gibbs_after_mcem_combine_chains'.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example",  package = "DGP4LCF")
#'
#' @return Convergence assessment for important continuous variables that do not need alignment, and posterior summary for predicted gene expressions.
#' @export
numerics_summary_do_not_need_alignment<- function(burnin = 0,
                                                  thin_step = 1,
                                                  pred_x_truth_indicator = FALSE,
                                                  pred_x_truth = NULL,
                                                  gibbs_after_mcem_combine_chains_result){

  p<- gibbs_after_mcem_combine_chains_result$p
  k<- gibbs_after_mcem_combine_chains_result$k
  n<- gibbs_after_mcem_combine_chains_result$n
  q<- gibbs_after_mcem_combine_chains_result$q
  tot_chain<- gibbs_after_mcem_combine_chains_result$tot_chain
  num_time_test<- gibbs_after_mcem_combine_chains_result$num_time_test
  ind_x<- gibbs_after_mcem_combine_chains_result$ind_x
  num_sample = gibbs_after_mcem_combine_chains_result$num_sample
  phi_final_array = gibbs_after_mcem_combine_chains_result$phi
  pred_indicator = gibbs_after_mcem_combine_chains_result$pred_indicator

  if (ind_x){
    individual_mean_final_array = gibbs_after_mcem_combine_chains_result$individual_mean
    variance_g_final_array = gibbs_after_mcem_combine_chains_result$variance_g
  }

  if (pred_indicator){

    pred_x_final_array = gibbs_after_mcem_combine_chains_result$pred_x

  }

  rm(gibbs_after_mcem_combine_chains_result)

  # retained samples
  tot_sample<- (num_sample - burnin)

  if (tot_sample%%2 == 0){
    keep_set1<- seq(from = 1, to = floor(tot_sample/2), by = 1)
    keep_set2<- seq(from = (floor(tot_sample/2) + 1), to = tot_sample, by = 1)
  } else {
    keep_set1<- seq(from = 1, to = floor(tot_sample/2), by = 1)
    keep_set2<- seq(from = (floor(tot_sample/2) + 1), to = (tot_sample-1), by = 1)
  }

  # print(paste0("This is tot_sample to be used for each chain: ",tot_sample))

  keep_set_index<- c(keep_set1, keep_set2)

  tot_sample_after_thinning<- length(c(keep_set1,keep_set2))

  # print(paste0("This is tot_sample_after_thinning: ",  tot_sample_after_thinning))

  # create a matrix that is used to restore the result

  if(ind_x & pred_indicator){
    convergence_result_table<- matrix(0, nrow = 4, ncol = 2)
    rownames(convergence_result_table)<- c("rhat_pred_x", "rhat_phi", "rhat_pred_individual_mean", "rhat_variance_g")
  } else if (!ind_x & pred_indicator){
    convergence_result_table<- matrix(0, nrow = 2, ncol = 2)
    rownames(convergence_result_table)<- c("rhat_pred_x", "rhat_phi")
  } else if (ind_x & !pred_indicator){
    convergence_result_table<- matrix(0, nrow = 3, ncol = 2)
    rownames(convergence_result_table)<- c("rhat_phi", "rhat_pred_individual_mean", "rhat_variance_g")
  } else {
    convergence_result_table<- matrix(0, nrow = 1, ncol = 2)
    rownames(convergence_result_table)<- c("rhat_phi")
  }

  colnames(convergence_result_table)<- c("Max Rhat", "Proportion of Rhat larger than 1.2")

  convergence_result_table_index<- 0

  ######################################################################################################
  ########################################### convergence assessment ###################################
  ######################################################################################################

  if (pred_indicator){

    ########################################### rhat for predicted x ###########################################

    rhat_pred_x<- matrix(0, nrow = (n*p*num_time_test), ncol = 5)
    colnames(rhat_pred_x)<- c("person index", "biomarker index", "time index", "point estimate", "upper")

    row_index<- 0

    for (person_index in 1:n){
      #print(paste0("This is for person_index: ",  person_index))
      for (biomarker_index in 1:p){
        for (time_index in 1:num_time_test){

          row_index<- row_index + 1

          rhat_pred_x[row_index,1]<- person_index
          rhat_pred_x[row_index,2]<- biomarker_index
          rhat_pred_x[row_index,3]<- time_index

          mh.list1<-   list(as.mcmc(pred_x_final_array[time_index, biomarker_index, person_index, keep_set1,1]),
                            as.mcmc(pred_x_final_array[time_index, biomarker_index, person_index, keep_set2,1]))

          for (chain_index in 2:tot_chain){
            mh.list1<- c(mh.list1,
                         list(as.mcmc(pred_x_final_array[time_index, biomarker_index, person_index, keep_set1,chain_index]),
                              as.mcmc(pred_x_final_array[time_index, biomarker_index, person_index, keep_set2,chain_index])))
          }

          mh.list <- mcmc.list(mh.list1)

          rhat_pred_x[row_index,(4:5)]<- as.numeric(gelman.diag(mh.list)$psrf)
        }
      }
    }

    convergence_result_table_index<- convergence_result_table_index + 1
    convergence_result_table[convergence_result_table_index, 1]<- round(max(rhat_pred_x[,4]),2)
    convergence_result_table[convergence_result_table_index, 2]<- round(sum((rhat_pred_x[,4])>1.2)/row_index,2)

    #print(paste0("This is convergence summary for x: ",  convergence_result_table[convergence_result_table_index,]))


    ######################################################################################################
    ########################### posterior summary for predicted x ###################################
    ######################################################################################################

    tot_sample_all_chains<- (tot_sample_after_thinning*tot_chain)

    pred_xt<- array(0,dim=c(tot_sample_all_chains, n, p, num_time_test))

    pred_x_final_array_original<- pred_x_final_array[,,,keep_set_index,]

    for (initial_index in 1:tot_chain){
      for (person_index in 1:n){
        for (time_index in 1:num_time_test){
          pred_xt[((tot_sample_after_thinning*(initial_index-1)+1):(tot_sample_after_thinning*initial_index)),person_index,,time_index]<-
            t(pred_x_final_array_original[time_index,,person_index,,initial_index]) # r.h.s. tot_sample_after_thinning * p
        }
      }
    }

    pred_lower_quantile<- array(0,dim=c(n,p,num_time_test)) # used to store lower quantile of predicted intervals - need to be used when drawing the final result
    pred_upper_quantile<- array(0,dim=c(n,p,num_time_test)) # used to store upper quantile of predicted intervals - need to be used when drawing the final result
    pred_median_quantile<- array(0,dim=c(n,p,num_time_test))

    for (person_index in 1:n){
      for (biomarker_index in 1:p){
        for (time_index in 1:num_time_test){
          # 95% credible intervals
          pred_lower_quantile[person_index,biomarker_index,time_index]<- as.numeric(quantile(pred_xt[,person_index,biomarker_index,time_index],c(0.025,0.50,0.975)))[1]
          pred_upper_quantile[person_index,biomarker_index,time_index]<- as.numeric(quantile(pred_xt[,person_index,biomarker_index,time_index],c(0.025,0.50,0.975)))[3]

          # median and error
          pred_median_quantile[person_index,biomarker_index,time_index]<- as.numeric(quantile(pred_xt[,person_index,biomarker_index,time_index],c(0.025,0.50,0.975)))[2]


        }
      }
    }

    if (pred_x_truth_indicator){

      pred_median_err<- abs(pred_median_quantile - pred_x_truth)

      sum_truth_within_interval_counter<- sum((pred_x_truth > pred_lower_quantile) & (pred_x_truth < pred_upper_quantile))

      #print(paste0("This is mae for x: ",  mean(pred_median_err)))
      #print(paste0("This is mwi for x: ",  mean(abs(pred_upper_quantile-pred_lower_quantile))))
      #print(paste0("This is pwi for x: ",  sum_truth_within_interval_counter/(n*p*num_time_test)))

      pred_x_result_list<- list(mae_using_median_est = mean(pred_median_err),
                                proportion_of_within_interval_biomarkers = sum_truth_within_interval_counter/(n*p*num_time_test),
                                mean_width_interval = mean(abs(pred_upper_quantile-pred_lower_quantile)),
                                pred_lower_quantile =  pred_lower_quantile,
                                pred_upper_quantile =  pred_upper_quantile,
                                pred_median_quantile = pred_median_quantile)

    } else {

      pred_x_result_list<- list(pred_lower_quantile =  pred_lower_quantile,
                                pred_upper_quantile =  pred_upper_quantile,
                                pred_median_quantile = pred_median_quantile)

    }

  }

  ########################################### rhat for phi ###########################################
  rhat_phi<- matrix(0, nrow = p, ncol = 3)
  colnames(rhat_phi)<- c("biomarker index", "point estimate", "upper")

  row_index_phi<- 0


  for (biomarker_index in 1:p){

    row_index_phi<- row_index_phi+ 1

    rhat_phi[ row_index_phi, 1]<- biomarker_index

    mh.list1<-   list(as.mcmc(phi_final_array[biomarker_index, keep_set1,1]),
                      as.mcmc(phi_final_array[biomarker_index, keep_set2,1]))

    for (chain_index in 2:tot_chain){
      mh.list1<- c(mh.list1,
                   list(as.mcmc(phi_final_array[biomarker_index, keep_set1,chain_index]),
                        as.mcmc(phi_final_array[biomarker_index, keep_set2,chain_index])))
    }

    mh.list <- mcmc.list(mh.list1)

    rhat_phi[row_index_phi,(2:3)]<- as.numeric(gelman.diag(mh.list)$psrf)

  }

  convergence_result_table_index<- convergence_result_table_index + 1
  convergence_result_table[convergence_result_table_index, 1]<- round(max(rhat_phi[,2]),2)
  convergence_result_table[convergence_result_table_index, 2]<- round(sum((rhat_phi[,2])>1.2)/row_index_phi,2)

  # print(paste0("This is convergence summary for phi: ",  convergence_result_table[convergence_result_table_index,]))

  if (ind_x){

    ########################################### rhat for subject-gene mean ###########################################
    rhat_pred_individual_mean<- matrix(0, nrow = (p*n), ncol = 4)
    colnames(rhat_pred_individual_mean)<- c("person index", "biomarker index", "point estimate", "upper")

    row_index_individual_mean<- 0

    for (person_index in 1:n){
      for (biomarker_index in 1:p){


        row_index_individual_mean<- row_index_individual_mean + 1

        rhat_pred_individual_mean[row_index_individual_mean,1]<- person_index
        rhat_pred_individual_mean[row_index_individual_mean,2]<- biomarker_index

        mh.list1<-   list(as.mcmc(individual_mean_final_array[biomarker_index, person_index, keep_set1,1]),
                          as.mcmc(individual_mean_final_array[biomarker_index, person_index, keep_set2,1]))

        for (chain_index in 2:tot_chain){
          mh.list1<- c(mh.list1,
                       list(as.mcmc(individual_mean_final_array[biomarker_index, person_index, keep_set1,chain_index]),
                            as.mcmc(individual_mean_final_array[biomarker_index, person_index, keep_set2,chain_index])))
        }


        mh.list <- mcmc.list(mh.list1)

        rhat_pred_individual_mean[row_index_individual_mean,(3:4)]<- as.numeric(gelman.diag(mh.list)$psrf)

      }
    }

    convergence_result_table_index<- convergence_result_table_index + 1
    convergence_result_table[convergence_result_table_index, 1]<- round(max(rhat_pred_individual_mean[,3]),2)
    convergence_result_table[convergence_result_table_index, 2]<- round(sum((rhat_pred_individual_mean[,3])>1.2)/row_index_individual_mean,2)

    # print(paste0("This is convergence summary for individual_mean: ",  convergence_result_table[convergence_result_table_index,]))

    ########################################### rhat for variance_g ###########################################
    rhat_variance_g<- matrix(0, nrow = p, ncol = 3)
    colnames(rhat_variance_g)<- c("biomarker index", "point estimate", "upper")

    row_index_variance_g<- 0


    for (biomarker_index in 1:p){

      row_index_variance_g<- row_index_variance_g+ 1

      rhat_variance_g[ row_index_variance_g, 1]<- biomarker_index

      mh.list1<-   list(as.mcmc(variance_g_final_array[biomarker_index, keep_set1,1]),
                        as.mcmc(variance_g_final_array[biomarker_index, keep_set2,1]))

      for (chain_index in 2:tot_chain){
        mh.list1<- c(mh.list1,
                     list(as.mcmc(variance_g_final_array[biomarker_index, keep_set1,chain_index]),
                          as.mcmc(variance_g_final_array[biomarker_index, keep_set2,chain_index])))
      }

      mh.list <- mcmc.list(mh.list1)

      rhat_variance_g[row_index_variance_g,(2:3)]<- as.numeric(gelman.diag(mh.list)$psrf)

    }

    convergence_result_table_index<- convergence_result_table_index + 1
    convergence_result_table[convergence_result_table_index, 1]<- round(max(rhat_variance_g[,2]),2)
    convergence_result_table[convergence_result_table_index, 2]<- round(sum((rhat_variance_g[,2])>1.2)/row_index_variance_g,2)

    # print(paste0("This is convergence summary for variance_g: ",  convergence_result_table[convergence_result_table_index,]))
  }

  ##################################################### return results ##############################################################################
  # pred_x_result_list: point estimate and intervals
  if (pred_indicator){
    result_list<- list(convergence_summary = convergence_result_table,
                       pred_x_result = pred_x_result_list)
  } else {
    result_list<- list(convergence_summary = convergence_result_table)
  }

  return(result_list)
}
