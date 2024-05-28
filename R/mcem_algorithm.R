#' Monte Carlo Expectation Maximization (MCEM) algorithm to return the Maximum Likelihood Estimate (MLE) of DGP Parameters.
#'
#' @description This function is used to return the MLE of DGP parameters.
#'
#' @param ind_x A logical value. ind_x = TRUE uses the model including the intercept term for subject-gene mean in within-MCEM-Gibbs sampler; otherwise uses the model without the intercept term.
#' @param ig_parameter A numeric scalar. Hyper-parameters for the prior Inverse-Gamma distribution.
#' @param increasing_rate A numeric scalar. Rate of increasing the sample size.
#' @param prob_conf_interval A numeric scalar. The probability that the true change in the Q-function is larger than the lower bound.
#' @param iter_count_num A numeric scalar. Maximum number of increasing the sample size; a larger number than this would end the algorithm.
#' @param x  A list of n elements. Each element is a matrix of dimension (p, q_i), storing the gene expression observed at q_i time points for the ith subject.
#' @param mcem_parameter_setup_result A list of objects returned from the function 'mcem_parameter_setup'.
#' @param ipt_x A logical value. ipx_x = TRUE denotes the need to impute for NAs of gene expression. The default value is ipt_x = FALSE.
#' @param missing_num A vector of n elements. Each element corresponds to a single person's number of NAs that needs imputation.
#' @param missing_list A list of n elements. Each element is a matrix of dimension (missing_num, 2): each row corresponds to the position of one NA that needs imputation; first and second columns denote the row and column indexes, respectively, of the NA in the corresponding person's matrix of gene expression.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#' vignette("bsfadgp_irregular_data_example", package = "DGP4LCF")
#'
#' @return The MLE of DGP parameters.
#' @export
mcem_algorithm<- function(ind_x,
                          ig_parameter = 10^-2,
                          increasing_rate = 0.5,
                          prob_conf_interval = 0.90,
                          iter_count_num = 5,
                          x,
                          mcem_parameter_setup_result,
                          ipt_x = FALSE,
                          missing_list = NULL,
                          missing_num = NULL){

  ################################################################################################################################################
  ############################################ preparation for running the main body of the algorithm ############################################
  ################################################################################################################################################
  # assign results to objects
  p<- mcem_parameter_setup_result$p
  k<- mcem_parameter_setup_result$k
  n<- mcem_parameter_setup_result$n
  q<- mcem_parameter_setup_result$q
  prior_sparsity<- mcem_parameter_setup_result$prior_sparsity
  model_dgp<- mcem_parameter_setup_result$model_dgp

  ind_num<- mcem_parameter_setup_result$ind_num
  burn_in_prop <- mcem_parameter_setup_result$burn_in_prop
  thin_step <- mcem_parameter_setup_result$thin_step
  mc_num<- mcem_parameter_setup_result$mc_num
  em_num<- mcem_parameter_setup_result$em_num

  a_full<- mcem_parameter_setup_result$a_full
  train_index<- mcem_parameter_setup_result$train_index
  a_person = mcem_parameter_setup_result$a_person
  obs_time_index = mcem_parameter_setup_result$obs_time_index
  obs_time_num = mcem_parameter_setup_result$obs_time_num

  latent_y <- mcem_parameter_setup_result$latent_y
  big_a <- mcem_parameter_setup_result$big_a
  big_z <- mcem_parameter_setup_result$big_z
  phi <- mcem_parameter_setup_result$phi
  beta <- mcem_parameter_setup_result$beta
  pai <- mcem_parameter_setup_result$pai
  hyper_record <- mcem_parameter_setup_result$hyper_record
  sigmay_record <- mcem_parameter_setup_result$sigmay_record
  sigmay_inv_record <- mcem_parameter_setup_result$sigmay_inv_record

  if (ind_x){

    individual_mean_in <- mcem_parameter_setup_result$individual_mean_in
    variance_g <- mcem_parameter_setup_result$variance_g
    mu_g <- mcem_parameter_setup_result$mu_g

  } else {

    individual_mean_in <- matrix(0, nrow = p, ncol = n)
    variance_g<- rep(0, times = p)
    mu_g<- rep(0, times = p)

  }

  # remove mcem_parameter_setup_result after assignment finished
  rm(mcem_parameter_setup_result)

  # ensure the input type of variabels 'missing_num' and 'missing_list' always matches those required by the c++ code
  if (ipt_x == TRUE) {

    if (is.null(missing_list)) {
      stop("mcem_algorithm(): parameter 'missing_list' needs to be provided if ipt_x is true")
    }

    if (is.null(missing_num)) {
      stop("mcem_algorithm(): parameter 'missing_num' needs to be provided if ipt_x is true")
    }

  } else {
    missing_list <- list()
    missing_num <- vector()
  }

  # constants

  e0<- prior_sparsity*p
  f0<- (1-prior_sparsity)*p

  c0<- ig_parameter
  d0<- ig_parameter
  c1<- ig_parameter
  d1<- ig_parameter
  c2<- ig_parameter
  d2<- ig_parameter

  big_z_table<- as.matrix(table_generator(k))

  a_train<- a_full[train_index]

  # for gpfda fitting
  if (model_dgp){

    h3n2_data<- list()

    list_temp <- vector("list", k)

    for (list_index in 1:k){
      list_temp[[list_index]]<- a_train
    }

    h3n2_data$input<- list_temp

  } else {

   h3n2_data_igp<- list()

   h3n2_data_igp$input<- list(a_train)
  }

  # if (ipt_x){
  #
  #   ########################################### obtain NA positions and count the number ###########################################
  #   missing_list<- vector("list", n)
  #   missing_num<- rep(0, times = n)
  #
  #   for (person_index in 1:n){
  #
  #     missing_list[[person_index]]<- which(is.na(x[[person_index]]),
  #                                          arr.ind = TRUE)
  #
  #     missing_num[person_index]<- nrow(which(is.na(x[[person_index]]),
  #                                            arr.ind = TRUE))
  #   }
  #
  #   ########################################### impute for x using subject-gene mean (initialization) ###########################################
  #   x_impute<- x
  #
  #   for (person_index in 1:n){
  #     if (missing_num[person_index] > 0){
  #       for (missing_index in 1:missing_num[person_index]){
  #
  #         # obtain the position of the NA value
  #         row_index<- missing_list[[person_index]][missing_index,1]
  #         column_index<- missing_list[[person_index]][missing_index,2]
  #
  #         # calculate subject-gene mean
  #         temp_mean<- mean(x[[person_index]][row_index,],na.rm = TRUE)
  #
  #         # impute using the subject-gene mean
  #         x_impute[[person_index]][row_index, column_index]<- temp_mean
  #
  #       }
  #     } else {
  #       next
  #     }
  #   }
  #
  #   x<- x_impute
  #
  # }

  # if ind_x = TRUE, just keep the original scale

  # if ind_x = FALSE, meaning that the model does not include the intercept term, then need to use center_within_individual_x

  if (!ind_x){

    x_center_within_individual<- x

    for (person_index in 1:n){
       for (gene_index in 1:p){
         x_center_within_individual[[person_index]][gene_index, ]<-
           x[[person_index]][gene_index, ] - mean(x[[person_index]][gene_index, ], na.rm = TRUE)
       }
    }

    x<- x_center_within_individual

  }

  # column names for applying rsa_exact algorithm
  col_name<- rep(0, times = (p*k))
  for (biomarker_index in 1:p){
    for (factor_index in 1:k){
      col_name[(((biomarker_index-1)*k)+factor_index)]<- paste0("LambdaV", biomarker_index,"_", factor_index)
    }
  }

  #######################################################################################################################################################
  ########################################### main body of the MCEM algorithm ###########################################################################
  #######################################################################################################################################################
  iter_count_tot<- 0

  # calculate the speed of each main part
  gibbs_comp<- rep(0, times = em_num)
  gpfda_comp<- rep(0, times = em_num)
  relabel_comp<- rep(0, times = em_num)
  lb_comp<- rep(0, times = em_num)

  for (i_em in 2:em_num){

    mc_num_old<- 1 # starting index of the gibbs sampler
    lower_bound<- (-1)  # to force the while loop happen if em algorithm needs to continue
    # Algorithm 1: 2.1.
    iter_count<- 0 # count # of increasing the sample size within one EM iteration
    index_used<- 0 # for each EM iteration: the default value of this is 0; it is kept as 0 until the algorithm reaches the point to end.

    # initial gibbs sample under new iterated values should always use the most recent gibbs sample value
    if (i_em>2){

      latent_y[[1]]<- latent_y[[mc_num]]
      big_a[[1]]<- big_a[[mc_num]]
      big_z[[1]]<-  big_z[[mc_num]]

    }

    while(lower_bound<0){

      iter_count<- iter_count+1

      # if it is more than once within this iteration: need to increase the sample size

      if (iter_count>1) {

        mc_num_old<- mc_num

        ind_num<- floor((1 + increasing_rate)*ind_num)

        mc_num<- floor((thin_step/(1-burn_in_prop))*ind_num)

        iter_count_tot<- iter_count_tot + 1 # count the number of increasing 'ind_num'

        # end the whole algorithm if  iter_count_tot exceeds a pre-specified number
        if (iter_count_tot> iter_count_num){
          message(paste0("Attempt number in total: ", iter_count_tot))
          message("The total number of increasing the sample size has exceeded the pre-specified number; algorithm ends.")
          index_used<- (i_em-1)
          message(paste0("The final result is stored in this iteration: ", index_used))
          break
        }
        # end the whole algorithm if iter_count exceeds a pre-specified number

        # recreate empty containers

        big_at<- list()
        big_zt<- list()
        latent_yt<- list()

        big_at[[mc_num]]<- 0
        big_zt[[mc_num]]<- 0
        latent_yt[[mc_num]]<- 0

        # assign multiple values of a list to another list
        for (list_index in 1:mc_num_old){
          latent_yt[[list_index]]<- latent_y[[list_index]] # q*k*n
          big_at[[list_index]]<- big_a[[list_index]]
          big_zt[[list_index]]<- big_z[[list_index]]

        }

        # for using the algorithm
        latent_y <- latent_yt
        big_a <- big_at
        big_z <- big_zt

      }


      ################ Algorithm 1: 2.2. gibbs sampling to obtain posterior samples - for the purpose to update gp parameters ##########################
      ### construct individual-specific objects first
      # create empty objects

      full_person_index<- which(obs_time_num==q) # index of person who have complete observations
      full_person_num<- length(full_person_index)

      missing_person_index<- which(obs_time_num!=q)
      missing_person_num<- length(missing_person_index)

      missing_time_num<- rep(0, times =  missing_person_num)
      missing_time_index <- vector("list",  missing_person_num)
      prod_covnew_covestinv<- vector("list",  missing_person_num)
      cov_cond_dist<- vector("list",  missing_person_num)
      sigmay_inv<- vector("list",  missing_person_num)

      if (missing_person_num > 0){

        for (list_index in 1:missing_person_num){

          person_index<- missing_person_index[list_index]

          result<- subject_specific_objects(k, q, a_train, a_person[[person_index]], sigmay_record[(i_em-1),,])

          # results returned from the algorithm
          missing_time_num[list_index]<- result$missing_time_num
          missing_time_index[[list_index]]<- result$missing_time_index
          prod_covnew_covestinv[[list_index]]<- result$prod_covnew_covestinv
          cov_cond_dist[[list_index]]<- result$cov_cond_dist
          sigmay_inv[[list_index]]<- result$cov_est_inv

        }
      }

      # x: a list of n elements, each element is a p*q_i matrix
      # missing_list: tells the positions of imputed values that needs to be updated in each Gibbs iteration

      start_time<- Sys.time()

        out <-
          gibbs_within_mcem_irregular_time(latent_y,
                                           x,
                                           missing_list, missing_num, ipt_x,
                                                big_a, big_z, phi, pai, beta, k, n, p, q, c0, c1, d0, d1, e0, f0, mc_num_old, mc_num,
                      big_z_table, ind_x, individual_mean_in, mu_g, variance_g, c2, d2,
                      obs_time_num,
                      obs_time_index,
                      sigmay_inv_record[(i_em-1),,],
                      missing_time_num,
                      missing_time_index,
                      prod_covnew_covestinv,
                      cov_cond_dist,
                      sigmay_inv,
                      missing_person_num,
                      missing_person_index,
                      full_person_num,
                      full_person_index)

        gibbs_comp[i_em]<-  Sys.time() - start_time

        # print("Gibbs Sampler Finished.")


      ############################### Algorithm 1: 2.3. reorder gibbs samples before feed it to obtain covariance matrix ####################
      burn_in<- floor(burn_in_prop*mc_num)
      sample_index<- seq(from=(burn_in+1), to=mc_num, by=thin_step)
      sample_num<- length(sample_index)

      reorder_y_array<- array(0, dim = c(sample_num, k, q, n))
      temp_matrix<- matrix(0, nrow = sample_num, ncol = (p*k))
      colnames(temp_matrix)<- col_name

      # need all mcmc draws to obtain reorderedPosterior
      for (i in 1:sample_num) {
        sample_index_temp<- sample_index[i]
        temp_matrix[i, ]<- c(t(big_a[[sample_index_temp]]*big_z[[sample_index_temp]]))
      }

      # print("Relabling for factor loadings and scores...")

      start_time<- Sys.time()
      reorderedPosterior <-  factor.switching::rsp_exact(lambda_mcmc = temp_matrix, rotate = FALSE)
      relabel_comp[i_em]<- Sys.time() - start_time

      response <- list()

      if (!model_dgp){
        cov_matrix<- matrix(0, nrow=(k*q), ncol=(k*q))
        hp_record<- matrix(0, nrow=k, ncol=5) # each process has 5 parameters
      }

      for (j in 1:k){

        latent_y_list<- list() # each element is a q*n matrix
        #sample_count<- 0
        #for (i in sample_index){
        for (i in 1:sample_num){
          #sample_count<- sample_count+1
          sample_index_temp<- sample_index[i]

          #rot_matrix<- varimax((big_a[[sample_index_temp]]*big_z[[sample_index_temp]]), normalize = FALSE)$rotmat
          permutation_matrix<- diag(length(reorderedPosterior$permute_vectors[i,]))[reorderedPosterior$permute_vectors[i,],]
          sign_matrix<- diag(reorderedPosterior$sign_vectors[i,])

          # transformation matrix
          # q_matrix<- rot_matrix %*% permutation_matrix %*% sign_matrix
          q_matrix<- permutation_matrix %*% sign_matrix
          q_matrix_inv<- solve(q_matrix) # apply this to factor score can obtain corresponding transformed factor scores

          # latent_y_list[[sample_count]]<- latent_y[[i]][,j,] # list - each element is a q*k*n matrix
          # latent_y_list[[i]]<- q_matrix_inv%*%latent_y[[sample_index_temp]][,j,]

          for (person_index in 1:n){
            reorder_y_array[i,,,person_index]<- q_matrix_inv%*%t(latent_y[[sample_index_temp]][,,person_index]) # k*q matrix
          }

          latent_y_list[[i]]<- reorder_y_array[i,j,,]

        }

        response[[j]]<- do.call(cbind, latent_y_list)

        # IGP model only: for the jth response
        if (!model_dgp){

          response_temp<- list()
          response_temp[[1]]<- response[[j]]
          h3n2_data_igp$response<- response_temp

          res_temp<- mgpr(Data=h3n2_data_igp,meanModel = 0) # q*q-dimensional matrix

          cov_matrix[((j-1)*q+1):(j*q),((j-1)*q+1):(j*q)]<- res_temp$Cov
          hp_record[j,]<- res_temp$hyper
        }

      }

      if (model_dgp){

        h3n2_data$response<- response

        start_time<- Sys.time()

        h3n2_res <- GPFDA::mgpr(Data=h3n2_data, meanModel = 0)

        gpfda_comp[i_em]<- Sys.time() - start_time

        hyper_record[i_em,]<- h3n2_res$hyper

        sigmay_record[i_em,,]<- cov2cor(h3n2_res$Cov) # apply the constraint that the diagonal elements should be 1

      } else {

        hyper_record[i_em,,]<- hp_record

        sigmay_record[i_em,,]<- cov2cor(cov_matrix)

      }

      sigmay_inv_record[i_em,,]<-  solve(sigmay_record[i_em,,])

      ################# Algorithm 1: 2.4. calculate asymptotic lower bound to decide whether or not to accept this updated \theta #################
      start_time<- Sys.time()

      batch_size<- 1
      batch_num<- sample_num

      delta_batch_sum<- rep(0, times=batch_num)


        for (i in 1:batch_num){
          for (j in 1:n){
            delta_batch_sum[i]<- delta_batch_sum[i]+
              mvtnorm::dmvnorm(c(t(reorder_y_array[i,,,j])), rep(0, times=(q*k)), sigmay_record[i_em,,], log=TRUE)-
              mvtnorm::dmvnorm(c(t(reorder_y_array[i,,,j])), rep(0, times=(q*k)), sigmay_record[(i_em-1),,], log=TRUE)
          }
        }

      delta_batch_mean<- delta_batch_sum/batch_size
      delta_q_mean<- sum(delta_batch_sum)/(batch_size*batch_num)

      sigma2_est<- var(delta_batch_mean) # directly calculating the sample variance leads to the same result as sigma2_est when batch size is set to be 1

      ase<- sqrt(sigma2_est/(batch_size*batch_num))
      lower_bound<- delta_q_mean-(qnorm(prob_conf_interval)*ase)

      lb_comp[i_em]<- Sys.time() - start_time

      # check_list<- list(sigmay_record,
      #                   q,
      #                   k,
      #                   reorder_y_array,
      #                   latent_y,
      #                   big_a,
      #                   big_z,
      #                   delta_batch_mean,
      #                   delta_q_mean,
      #                   batch_size,
      #                   batch_num,
      #                   h3n2_data,
      #                   h3n2_res)
      #
      # return(check_list)

      ############################################### print some results ##############################################################################
      ### algorithm progress
      message(paste("EM iteration number:",(i_em-1)))
      message(paste("Lower bound:", lower_bound/(k*q*n)))
      message(paste("Attempt number within this EM iteration:", iter_count))
      message(paste("Attempt number in total:", iter_count_tot))

      ### computation complexity
      # print(paste("This is the sample size:",n))
      # print(paste("This is the independent sample size:",ind_num))
      # print(paste("This is the monte carlo sample size:",mc_num))

      # print(paste("This is the time cost for gibbs:",gibbs_comp[i_em]))
      # print(paste("This is the time cost for relabeling:",relabel_comp[i_em]))
      # print(paste("This is the time cost for gpfda:",gpfda_comp[i_em]))
      # print(paste("This is the time cost for lb:",lb_comp[i_em]))

      ################################ if not factor first but time first, then need to reorder the covariance matrix for next Gibbs #########################

    }

    # if one iterations successfully updates dgp parameters within the allowed number: index_used=0, hence algo will not exit and will continue to the next EM iteration
    # if one iterations fails to update dgp parameters within the allowed number: index_used!=0, hence algo will jump out the 'for' loop
    if (index_used!=0){
      break
    }
  }

  #######################################################################################################################################################
  ################################################################# return results ######################################################################
  #######################################################################################################################################################
  result_list<- list(index_used = index_used,
                     hyper_record = hyper_record,
                     sigmay_record = sigmay_record,
                     sigmay_inv_record = sigmay_inv_record,
                     prior_sparsity = prior_sparsity,
                     ig_parameter = ig_parameter,
                     missing_list = missing_list,
                     missing_num = missing_num,
                     ipt_x = ipt_x,
                     gibbs_comp = gibbs_comp,
                     relabel_comp = relabel_comp,
                     gpfda_comp = gpfda_comp,
                     lb_comp = lb_comp)

  return(result_list)
}

