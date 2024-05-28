#' @title Numerical summary for factor loadings and factor scores, which need alignment.
#'
#' @details This function corresponds to Algorithm 2: Steps 2, 3 and 4 in the main manuscript; therefore reader can consult the paper for more explanations.
#'
#' @param burnin A numeric scalar. The saved samples are already after burnin; therefore the default value for this parameter here is 0. Can discard further samples if needed.
#' @param thin_step A numeric scalar. The saved samples are already after thinning; therefore the default value for this parameter here is 1. Can be further thinned if needed.
#' @param gibbs_after_mcem_combine_chains_result A list of objects returned from the function 'gibbs_after_mcem_combine_chains'.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#' vignette("bsfadgp_irregular_data_example", package = "DGP4LCF")
#'
#' @return Reordered posterior samples, convergence assessment, and summarized posterior results for factor loadings and factor scores.
#' @export
numerics_summary_need_alignment<- function(burnin = 0,
                                           thin_step = 1,
                                           gibbs_after_mcem_combine_chains_result){

  p<- gibbs_after_mcem_combine_chains_result$p
  k<- gibbs_after_mcem_combine_chains_result$k
  n<- gibbs_after_mcem_combine_chains_result$n
  q<- gibbs_after_mcem_combine_chains_result$q
  tot_chain<- gibbs_after_mcem_combine_chains_result$tot_chain
  num_time_test<- gibbs_after_mcem_combine_chains_result$num_time_test
  num_time_all<- (q+num_time_test)

  num_sample = gibbs_after_mcem_combine_chains_result$num_sample
  big_a_final_array = gibbs_after_mcem_combine_chains_result$big_a
  big_z_final_array = gibbs_after_mcem_combine_chains_result$big_z
  latent_y_final_array = gibbs_after_mcem_combine_chains_result$latent_y

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

  keep_set_index<- c(keep_set1, keep_set2)

  tot_sample_after_thinning<- length(c(keep_set1,keep_set2))

  convergence_result_table<- matrix(0, nrow = 2, ncol = 2)
  rownames(convergence_result_table)<- c("rhat_l", "rhat_y")
  colnames(convergence_result_table)<- c("Max Rhat (Point Estimate)", "Proportion of Rhat larger than 1.2")

  ######################################################################################################################################################
  ##################################################### Reorder ##############################################################################
  ######################################################################################################################################################

  ######################################### obtain factor loading matrix l ###########################################
  big_l_final_array<- big_a_final_array*big_z_final_array

  ######################################### reorder factor loading matrix l ###########################################
  ### create original posterior chains following the software's requirement
  originalPosterior<- vector('list', length=tot_chain)

  col_name<- rep(0, times = (p*k))
  for (biomarker_index in 1:p){
    for (factor_index in 1:k){
      col_name[(((biomarker_index-1)*k)+factor_index)]<- paste0("LambdaV", biomarker_index,"_", factor_index)
    }
  }

  for (chain_index in 1:tot_chain){

    temp_matrix<- matrix(0, nrow = tot_sample_after_thinning, ncol = (p*k))
    colnames(temp_matrix)<- col_name

    for (sample_index in 1:tot_sample_after_thinning){
      temp_matrix[sample_index, ]<- c(t(big_l_final_array[,,keep_set_index[sample_index], chain_index]))
    }

    originalPosterior[[chain_index]]<- temp_matrix
  }

  ### reorder within each chain
  reorderedPosterior <- vector('list',length=tot_chain)
  for(chain_index in 1:tot_chain){
    print(paste0("Relabeling within Chain ", chain_index, " for Factor Loadings"))
    reorderedPosterior[[chain_index]] <- factor.switching::rsp_exact( lambda_mcmc = originalPosterior[[chain_index]],
                                                    maxIter = 100,
                                                    threshold = ((1e-6)),
                                                    verbose=TRUE,
                                                    rotate = FALSE)
  }

  ### reorder across chain
  print(paste0("Relabeling Across Chains for Factor Loadings"))
  makeThemSimilar <- factor.switching::compareMultipleChains(rspObjectList=reorderedPosterior)

  #################################### reorder for z, a, y, which is of interest ###########################################
  #################################### note: a,z's transformation is exactly the same as that of l ########################
  # within chain
  signed_permutation_matrix_within_chain<- array(0, dim = c(tot_chain, tot_sample_after_thinning, k, k))
  signed_permutation_matrix_within_chain_inv<- array(0, dim = c(tot_chain, tot_sample_after_thinning, k, k))

  big_z_final_array_reordered_within_chain<- array(0, dim = c(p,k,tot_sample_after_thinning,tot_chain))
  big_a_final_array_reordered_within_chain<- array(0, dim = c(p,k,tot_sample_after_thinning,tot_chain))
  latent_y_final_array_reordered_within_chain<- array(0, dim = c(num_time_all, k, n, tot_sample_after_thinning, tot_chain))

  # across chain
  signed_permutation_matrix_across_chain<- array(0, dim = c(tot_chain, k, k))
  signed_permutation_matrix_across_chain_inv<- array(0, dim = c(tot_chain, k, k))

  big_z_final_array_reordered_across_chain<- array(0, dim = c(p,k,tot_sample_after_thinning,tot_chain))
  big_a_final_array_reordered_across_chain<- array(0, dim = c(p,k,tot_sample_after_thinning,tot_chain))
  latent_y_final_array_reordered_across_chain<- array(0, dim = c(num_time_all, k, n, tot_sample_after_thinning, tot_chain))

  permutation_result<- combinat::permn(1:k)
  permutation_result_length<- length(permutation_result)

  # construct all possible combinations of signed variable
  l <- rep(list(c(-1,1)), k)

  sign_vector<- expand.grid(l)

  for (chain_index in 1:tot_chain){

    for (permu_index in 1:permutation_result_length){

      permutated_matrix_temp<- diag(length(permutation_result[[permu_index]]))[permutation_result[[permu_index]],]

      for (sign_index in 1:(2^k)){

        sign_matrix_temp<- diag(sign_vector[sign_index,],k)

        # consider a specfic q_temp
        q_temp<- permutated_matrix_temp%*%sign_matrix_temp
        q_temp_indicator <- 1 # flag

        for (sample_index in 1:tot_sample_after_thinning){

          transformed_matrix_temp<- matrix(as.numeric((reorderedPosterior[[chain_index]]$lambda_reordered_mcmc)[sample_index,]), nrow = p, ncol = k, byrow=T)[(1:k),(1:k)]%*%q_temp

          if (all.equal(transformed_matrix_temp, matrix(as.numeric(makeThemSimilar[[chain_index]][sample_index,]), nrow = p, ncol = k, byrow=T)[(1:k),(1:k)]) !=TRUE){
            #print(paste0("This is not the correct transformed matrix! stop sample test at sample ", sample_index))
            q_temp_indicator <- (-1)
            break # once a sample cannot be mapped, stop the following sample test
          }

        }

        if (q_temp_indicator == (-1)){
          next
        } else if (q_temp_indicator == 1){

          q_final<- q_temp

          #print("Congratulations! you have found the correct transformed matrix, which passed all sample test")
          #print("Current Permutation Matrix:")
          #print(permutated_matrix_temp)

          #print(("Current Signed Matrix:"))
          #print(sign_matrix_temp)

          #print(("Current Transformed Matrix:"))
          #print(q_final)

          break
        }
      }

      if (q_temp_indicator == (-1)){
        next
      } else if (q_temp_indicator == 1){
        break
      }
    }

    signed_permutation_matrix_across_chain[chain_index,,]<- q_final

    signed_permutation_matrix_across_chain_inv[chain_index,,]<- solve(signed_permutation_matrix_across_chain[chain_index,,])

    for (sample_index in 1:tot_sample_after_thinning){

      signed_permutation_matrix_within_chain[chain_index, sample_index,,]<- diag(1,k)
      signed_permutation_matrix_within_chain_inv[chain_index, sample_index,,]<- diag(1,k)

      big_z_final_array_reordered_within_chain[,,sample_index, chain_index]<- abs(big_z_final_array[,,keep_set_index[sample_index], chain_index]%*%signed_permutation_matrix_within_chain[chain_index, (sample_index),,])
      big_a_final_array_reordered_within_chain[,,sample_index, chain_index]<- big_a_final_array[,,keep_set_index[sample_index], chain_index]%*%signed_permutation_matrix_within_chain[chain_index, (sample_index),,]

      # align across chains
      big_z_final_array_reordered_across_chain[,,sample_index, chain_index]<- abs(big_z_final_array_reordered_within_chain[,,sample_index, chain_index]%*%signed_permutation_matrix_across_chain[chain_index,,])
      big_a_final_array_reordered_across_chain[,,sample_index, chain_index]<- big_a_final_array_reordered_within_chain[,,sample_index, chain_index]%*%signed_permutation_matrix_across_chain[chain_index,,]

      # confirm the transformation on a and z is correct by comparing them with final reordered lmakethemsimilar
      #print(all.equal((big_a_final_array_reordered_across_chain[,,sample_index,chain_index])*(big_z_final_array_reordered_across_chain[,,sample_index,chain_index]),
      #                matrix(as.numeric(makeThemSimilar[[chain_index]][sample_index,]), nrow = p, ncol = k, byrow=T)))

      for (person_index in 1:n){

        temp_y_kq<- signed_permutation_matrix_within_chain_inv[chain_index, sample_index,,]%*%t(latent_y_final_array[,,person_index, keep_set_index[sample_index], chain_index]) # k*q
        latent_y_final_array_reordered_within_chain[,, person_index, sample_index, chain_index]<- t(temp_y_kq) # q*k

        temp_y_kq<- signed_permutation_matrix_across_chain_inv[chain_index,,]%*%t(latent_y_final_array_reordered_within_chain[,,person_index, sample_index, chain_index]) # k*q
        latent_y_final_array_reordered_across_chain[,,person_index, sample_index, chain_index]<- t(temp_y_kq) # q*k

      }

      #print(paste0("Relabeling for Chain ", chain_index, " Sample ", sample_index))
    }
  }

  ######################################### convergence assessment and posterior summary for factor loadings l ###########################################
  ### obtain l in array format
  factor_loading_l_final_array<- array(0, dim = c(p, k, tot_sample_after_thinning, tot_chain))

  for (chain_index in 1:tot_chain){
    for (sample_index in 1:tot_sample_after_thinning){
      factor_loading_l_final_array[,,sample_index, chain_index]<- matrix(as.numeric(makeThemSimilar[[chain_index]][sample_index,]),
                                                                         nrow = p,
                                                                         ncol = k,
                                                                         byrow = T)
    }
  }

  ### convergence assessment table
  rhat_l<- matrix(0, nrow = (p*k), ncol = 4)
  colnames(rhat_l)<- c("factor index", "biomarker index", "point estimate", "upper")
  row_index_l<- 0

  # for those whose binary counterpart z_ga is 0 - we cannot calculate its rhat since it is meaningless
  factor_loading_l_final_array_summary_for_zero<- array(0, dim = c(p, k, tot_chain))
  factor_loading_final_array_summary_for_zero_binary<- matrix(1, nrow = p, ncol = k)
  row_index_l_excluded<- NULL

  # for those whose binary counterpart z_ga is not 0 - it is continuous variable so we can calculate its rhat and the retained indexes are:
  # if (tot_sample_after_thinning%%2 ==0){
  #  keep_set1_for_l<- seq(from = 1, to = floor(tot_sample_after_thinning/2), by = 1)
  #  keep_set2_for_l<- seq(from = (floor(tot_sample_after_thinning/2) + 1), to = tot_sample_after_thinning, by = 1)
  #} else {
  #  keep_set1_for_l<- seq(from = 1, to = floor(tot_sample_after_thinning/2), by = 1)
  #  keep_set2_for_l<- seq(from = (floor(tot_sample_after_thinning/2) + 1), to = (tot_sample_after_thinning-1), by = 1)
  #}

  #keep_set_index_for_l<- c(keep_set1_for_l, keep_set2_for_l)

  # for those whose binary counterpart z_ga is not 0: summarize the point estimate if it is significantly different from 0
  factor_loading_l_final_array_summary_for_nonzero<- matrix(0, nrow = p, ncol =k)
  factor_loading_l_final_array_summary_for_nonzero_upper<- matrix(0, nrow = p, ncol = k)
  factor_loading_l_final_array_summary_for_nonzero_lower<- matrix(0, nrow = p, ncol = k)

  for (biomarker_index in 1:p){
    for (factor_index in 1:k){

      row_index_l<-row_index_l + 1

      rhat_l[row_index_l,1]<- factor_index
      rhat_l[row_index_l,2]<- biomarker_index

      # if all chains say that this l_ga is 0, then there is no need to summarize for it as the result suggests that its binary counterpart is 0; otherwise cont to assess its rhat
      for (chain_index in 1:tot_chain){
        factor_loading_l_final_array_summary_for_zero[biomarker_index, factor_index, chain_index]<- (sum(factor_loading_l_final_array[biomarker_index, factor_index,,chain_index]==0)/tot_sample_after_thinning)
      }

      if (sum(factor_loading_l_final_array_summary_for_zero[biomarker_index, factor_index, ]> 0.5) == tot_chain){
        factor_loading_final_array_summary_for_zero_binary[biomarker_index, factor_index]<- 0
        row_index_l_excluded<- c(row_index_l_excluded,row_index_l)
        next
      } else {

        # summarize to see if it is significantly different from 0
        factor_loading_l_final_array_summary_for_nonzero_lower[biomarker_index, factor_index]<-
          as.numeric(quantile(factor_loading_l_final_array[biomarker_index, factor_index, ,],c(0.025,0.50,0.975)))[1]

        factor_loading_l_final_array_summary_for_nonzero_upper[biomarker_index, factor_index]<-
          as.numeric(quantile(factor_loading_l_final_array[biomarker_index, factor_index, ,],c(0.025,0.50,0.975)))[3]

        if (factor_loading_l_final_array_summary_for_nonzero_lower[biomarker_index, factor_index]<0 &  factor_loading_l_final_array_summary_for_nonzero_upper[biomarker_index, factor_index]>0){
          # not significant from 0: summary for this l_ga is set as default 0
          row_index_l_excluded<- c(row_index_l_excluded,row_index_l)
          next
        } else {

          # significant from 0: summarize for this l_ga using median estimate
          factor_loading_l_final_array_summary_for_nonzero[biomarker_index, factor_index]<-
            as.numeric(quantile(factor_loading_l_final_array[biomarker_index, factor_index, ,],c(0.025,0.50,0.975)))[2]

          ## revised: when assessing convergence for l_ga|Z_ga, should only consider the iterations when Z_ga = 1

          length_index_nonzero<- sum(factor_loading_l_final_array[biomarker_index, factor_index, ,1]!=0)

          for (chain_index in 2:tot_chain){

            length_index_nonzero<- min(length_index_nonzero, sum(factor_loading_l_final_array[biomarker_index, factor_index, ,chain_index]!=0))

          }

          if ( length_index_nonzero %% 2 == 0){
            length_index_nonzero<-  length_index_nonzero
          } else {
            length_index_nonzero<-  length_index_nonzero - 1
          }

          keep_set1_index_nonzero<- (1:(length_index_nonzero/2))
          keep_set2_index_nonzero<- ((length_index_nonzero/2)+1):length_index_nonzero

          ###################################################################################################################

          index_nonzero<- which(factor_loading_l_final_array[biomarker_index, factor_index, ,1]!=0)

          mh.list1<-   list(as.mcmc(factor_loading_l_final_array[biomarker_index, factor_index,   index_nonzero[keep_set1_index_nonzero], 1]),
                            as.mcmc(factor_loading_l_final_array[biomarker_index, factor_index,   index_nonzero[keep_set2_index_nonzero], 1]))

          for (chain_index in 2:tot_chain){

            index_nonzero<- which(factor_loading_l_final_array[biomarker_index, factor_index, ,chain_index]!=0)

            mh.list1<- c(mh.list1,
                         list(as.mcmc(factor_loading_l_final_array[biomarker_index, factor_index, index_nonzero[keep_set1_index_nonzero], chain_index]),
                              as.mcmc(factor_loading_l_final_array[biomarker_index, factor_index, index_nonzero[keep_set2_index_nonzero], chain_index])))
          }

          mh.list <- mcmc.list(mh.list1)

          rhat_l[row_index_l,(3:4)]<- as.numeric(gelman.diag(mh.list)$psrf)

        }
      }
    }
  }

  if (sum(is.na(rhat_l[,3])!=0)){
    print("There are NA values of rhat_l")
    rhat_l[is.na(rhat_l[,3]),3]<- 5 # assign NA values a very large value to suggest non-convergence
  } else {
    print("There are no NA values of rhat_l")
  }

  convergence_result_table_index<- 1
  convergence_result_table[convergence_result_table_index,1]<- round(max(rhat_l[-row_index_l_excluded,3]),2)
  convergence_result_table[convergence_result_table_index,2]<- round(sum((rhat_l[-row_index_l_excluded,3])>1.2)/(row_index_l-length(row_index_l_excluded)),2)

  print(paste0("This is convergence summary for rhat_l: ",  convergence_result_table[convergence_result_table_index,]))

  ########################################### convergence assessment and posterior summary for factor scores y ###########################################
  rhat_y<- matrix(0, nrow = (n*k*num_time_all), ncol = 5)
  colnames(rhat_y)<- c("person index", "factor index", "time index" , "point estimate", "upper")

  row_index_y<- 0

  latent_y_final_array_reordered_across_chain_summary<- array(0, dim = c(num_time_all, k, n))

  for (person_index in 1:n){
    for (factor_index in 1:k){
      for (time_index in 1:num_time_all){

        row_index_y<- row_index_y + 1

        rhat_y[row_index_y,1]<- person_index
        rhat_y[row_index_y,2]<- factor_index
        rhat_y[row_index_y,3]<- time_index

        mh.list1<-   list(as.mcmc(latent_y_final_array_reordered_across_chain[time_index, factor_index, person_index, keep_set1, 1]),
                          as.mcmc(latent_y_final_array_reordered_across_chain[time_index, factor_index, person_index, keep_set2, 1]))

        for (chain_index in 2:tot_chain){
          mh.list1<- c(mh.list1,
                       list(as.mcmc(latent_y_final_array_reordered_across_chain[time_index, factor_index, person_index, keep_set1,chain_index]),
                            as.mcmc(latent_y_final_array_reordered_across_chain[time_index, factor_index, person_index, keep_set2,chain_index])))
        }

        mh.list <- mcmc.list(mh.list1)

        rhat_y[row_index_y,(4:5)]<- as.numeric(gelman.diag(mh.list)$psrf)

        latent_y_final_array_reordered_across_chain_summary[time_index, factor_index, person_index]<-
          median(latent_y_final_array_reordered_across_chain[time_index, factor_index, person_index,,])

      }
    }
  }

  if (sum(is.na(rhat_y[,4])!=0)){
    print("There are NA values of rhat_y.")
    rhat_y[is.na(rhat_y[,4]),4]<- 5 # assign NA values a very large value to suggest non-convergence
  } else {
    print("There are no NA values of rhat_y.")
  }

  convergence_result_table_index<- convergence_result_table_index + 1
  convergence_result_table[convergence_result_table_index,1]<- round(max(rhat_y[,4]),2)
  convergence_result_table[convergence_result_table_index,2]<- round(sum((rhat_y[,4])>1.2)/row_index_y,2)

  print(paste0("This is convergence summary for latent_y: ",  convergence_result_table[convergence_result_table_index,]))

  ###################################################################################################################################################
  ##################################################### return results ##############################################################################
  ###################################################################################################################################################
  # posterior samples
  reordered_list<- list(latent_y = latent_y_final_array_reordered_across_chain,
                        big_z = big_z_final_array_reordered_across_chain,
                        big_a = big_a_final_array_reordered_across_chain,
                        big_l = factor_loading_l_final_array,
                        big_l_summary_for_zero_binary = factor_loading_final_array_summary_for_zero_binary,
                        big_l_summary_for_zero_cont = factor_loading_l_final_array_summary_for_zero)

  # posterior median estimate
  reordered_summary_list<- list(latent_y = latent_y_final_array_reordered_across_chain_summary,
                                big_l = factor_loading_l_final_array_summary_for_nonzero)

  # combine all results
  result_list<- list(reordered_samples = reordered_list,
                     reordered_summary =  reordered_summary_list,
                     convergence_summary = convergence_result_table)

  return(result_list)

}
