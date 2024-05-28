#' Parameters' setup and initial value assignment for the Monte Carlo Expectation Maximization (MCEM) algorithm.
#'
#' @description This function is used to create R objects storing parameters in the desired format, and assign initial values so that they are ready to use in the MCEM algorithm.
#'
#' @details
#' The following parameters are worth particular attention, and users should tune these parameters according to the specific data.
#'
#'  'burn_in_prop' and 'thin_step' co-control the number of Gibbs samples needed in order to generate approximately 'ind_num' independent samples.
#'  The ultimate purpose of tuning these two parameters is to generate high-quality posterior samples for latent factor scores.
#'  Therefore: if initials of the Gibbs sampler are not good, readers may need to increase 'burn_in_prop' to discard more burn-in samples; if high-correlation is a potential concern, 'thin_step' may need to be larger.
#'
#' @param p A numeric scalar. Number of genes.
#' @param k A numeric scalar. Number of latent factors.
#' @param n A numeric scalar. Number of subjects.
#' @param q A numeric scalar. Complete number of time points in the training data.
#' @param ind_num A numeric scalar. Starting size of approximately independent samples for MCEM.
#' @param burn_in_prop A numeric scalar. Proportion of burnin, which be used to calculate size of Monte Carlo samples needed in the Gibbs sampler. Must be the same as that in the function 'mcem_algorithm_irregular_time'.
#' @param thin_step A numeric scalar. Thinning step, which be used to calculate size of Monte Carlo samples needed in the Gibbs sampler. Must be the same as that in the function 'mcem_algorithm_irregular_time'.
#' @param prior_sparsity A numeric scalar. Prior expected proportion of genes involved within each pathway.
#' @param em_num A numeric scalar. Maximum iterations of the expectation maximization (EM) algorithm allowed.
#' @param obs_time_num A n-dimensional vector. One element represents one person's observed number of time points in the training data.
#' @param obs_time_index A list of n elements. One element is a vector of observed time indexes for one person in the training data, sorted from early to late.
#' @param a_person A list of n elements. One element is a vector of observed time for one subject in the training data, sorted from early to late.
#' @param col_person_index A list of n elements. One element is a vector of column indexes for one subject in y_init.
#' @param y_init A matrix of dimension (k, sum(obs_time_num)). Initial values of the latent factor score. Can be obtained using BFRM software.
#' @param a_init A matrix of dimension (p, k). Initial values of the regression coefficients of factor loadings. Can be obtained using BFRM software.
#' @param z_init A matrix of dimension (p, k). Initials values of the binary variables of factor loadings. Can be obtained using BFRM software.
#' @param phi_init A p-dimensional column vector. Initials values of the variance for residuals when modeling gene expressions, corresponding to \eqn{\frac{1}{\phi^2}}{1/(\phi^2)} in the manuscript. Can be obtained using BFRM software.
#' @param a_full A numeric vector. Complete time observed, sorted from early to late.
#' @param train_index A q-dimensional column vector. Index of time points used in the training data.
#' @param x A list of n elements. Each element is a matrix of dimension (p, q_i), storing the gene expressions for the ith subject.
#' @param model_dgp A logical value. model_dgp = TRUE (default setting) uses the Dependent Gaussian Process to model latent factor trajectories, otherwise the Independent Gaussian Process is used.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#' vignette("bsfadgp_irregular_data_example", package = "DGP4LCF")
#'
#' @return A list of R objects required in the MCEM algorithm.
#' @export
mcem_parameter_setup<- function(p, k, n, q, ind_num = 10, burn_in_prop = 0.20, thin_step = 5, prior_sparsity = 0.1, em_num = 50,
                                               obs_time_num,
                                               obs_time_index,
                                               a_person,
                                               col_person_index,
                                               y_init,
                                               a_init,
                                               z_init,
                                               phi_init,
                                               a_full,
                                               train_index,
                                               x,
                                               model_dgp = TRUE){

  # set.seed(seed)
  #################################################################################################################################################################
  ########################################### create empty objects for storing parameters ############################################################################
  #################################################################################################################################################################
  # decide the starting sample size for the gibbs sampler
  mc_num<- floor((thin_step/(1-burn_in_prop))*ind_num)

  # parameter arrays
  big_a<- list()
  big_z<- list()
  latent_y<- list()

  # need to set the length of lists so they can be used as C++ output
  big_a[[mc_num]]<- 0
  big_z[[mc_num]]<- 0
  latent_y[[mc_num]]<- 0

  a_train<- a_full[train_index] # training data

  #################################################################################################################################################################
  ###################################################### initial value assignment #################################################################################
  #################################################################################################################################################################

  ######################################################## load results from two-step approach  ###################################################################

  ################################################################# assign initial values to DGP parameters ######################################################

  h3n2_response_gpfda <- list()  # create desired data for gpfda fitting: q common training points among all n people, so need to transform each latent factor score into a matrix of dimension q*n

  full_num<- sum(obs_time_num==q) # how many people have complete observations

  if (full_num > 0){

    full_index<- which(obs_time_num==q) # index of person who have complete observations

    a_train_init<- a_train

  } else {

    full_index<- which.max(obs_time_num) # this will list only one person has the maximum number of observations
    full_num<- 1
    a_train_init<- a_person[[full_index]]
  }

  col_person_index_combine<- NULL

  for (person_index in full_index){
     col_person_index_combine<- c(col_person_index_combine, col_person_index[[person_index]])
  }

  for (l in 1:k){

    # transform a vector into a matrix (q*n)

    h3n2_response_gpfda[[l]]<- matrix(as.numeric(y_init[l, col_person_index_combine]), max(obs_time_num), ncol = full_num)

  }

  # assign time points following the chosen model

  if (model_dgp){

    h3n2_data_dgp<- list()

    list_temp <- vector("list", k)

    for (list_index in 1:k){
      list_temp[[list_index]]<- a_train_init
    }

    h3n2_data_dgp$input<- list_temp

  } else {

    h3n2_data_igp<- list()

    h3n2_data_igp$input<- list(a_train_init)
  }

  ### assign response according to the y

  if (model_dgp){

    h3n2_data_dgp$response<- h3n2_response_gpfda

    h3n2_res <- GPFDA::mgpr(Data=h3n2_data_dgp, meanModel = 0)

    gp_num<- length(h3n2_res$hyper)

    hyper_record<- matrix(0, nrow=em_num, ncol=gp_num)

    hyper_record[1,]<- h3n2_res$hyper


    if (max(obs_time_num) == q){

      psi_initial <- (cov2cor(h3n2_res$Cov)) # apply the constraint: variance of factors = 1

    } else {

      # need to construct psi_initial using full time

      h3n2_data_dgp_full<- list()

      list_temp <- vector("list", k)

      for (list_index in 1:k){

        list_temp[[list_index]]<- a_train

      }

      h3n2_data_dgp_full$input<- list_temp

      psi_initial<- mgpCovMat(Data = h3n2_data_dgp_full, hp = h3n2_res$hyper)

      psi_initial <- cov2cor(psi_initial)

    }

  } else {

    hyper_record<- array(0, dim=c(em_num,k,5))

    cov_matrix<- matrix(0, nrow=(k*q), ncol=(k*q))

    for (j in 1:k){

      response_temp<- list()

      response_temp[[1]]<- h3n2_response_gpfda[[j]]

      h3n2_data_igp$response<- response_temp

      res_temp<- mgpr(Data= h3n2_data_igp, meanModel = 0) # q*q-dimensional matrix

      hyper_record[1,j,]<- res_temp$hyper

      if (max(obs_time_num) == q){

        cov_matrix[((j-1)*q+1):(j*q),((j-1)*q+1):(j*q)]<- res_temp$Cov

      } else {

        h3n2_data_igp_full<- list()

        h3n2_data_igp_full$input<- list(a_train)

        cov_matrix[((j-1)*q+1):(j*q),((j-1)*q+1):(j*q)]<- mgpCovMat(Data = h3n2_data_igp_full, hp = res_temp$hyper)

      }
    }

    psi_initial <- cov2cor(cov_matrix)
  }

  ### assign results to the record matrix

  sigmay_record<- array(0,dim=c(em_num, (q*k), (q*k)))

  sigmay_inv_record<- array(0,dim=c(em_num, (q*k), (q*k)))

  sigmay_record[1,,]<- psi_initial

  sigmay_inv_record[1,,]<- solve(sigmay_record[1,,])

  # ################################################################# assign initial values to DGP parameters ######################################################

  ################################################### assign initial values to other parameters in the model #################################################

  big_a[[1]]<- a_init
  big_z[[1]]<- z_init

  # to obtain initials for y: only for people who have estimates from bfrm

  initial_y<- array(0,dim=c(q, k, n))

  for (person_index in 1:n){
    initial_y[obs_time_index[[person_index]],,person_index]<- t(data.matrix(y_init[,col_person_index[[person_index]]])) # q_i*k
  }

  latent_y[[1]]<- initial_y # q*k*n array

  phi<- phi_init # p*1 vector
  beta<- rep(1, times = k) # self-assign random initials as this is not directly available from BFRM

  pai<- rep(0, times = k)
  e0<- prior_sparsity*p
  f0<- (1-prior_sparsity)*p

  for (i in 1:k){
    sum_i<- sum(big_z[[1]][,i])
    e<- e0+ sum_i
    f<- f0+ (p-sum_i)
    pai[i]<- rbeta(1,e,f)
  }

  # need to include extra parameters if want to use the model with the intercept term

    gene_person_mean_est<- matrix(0, nrow = p, ncol =n)
    for (gene_index in 1:p){
      for (person_index in 1:n){
        gene_person_mean_est[gene_index, person_index]<- mean(x[[person_index]][gene_index,], na.rm = TRUE)
      }
    }

    individual_mean_in<- gene_person_mean_est

    mu_g<- rowMeans(gene_person_mean_est)

    variance_g<- apply(gene_person_mean_est, 1, var)

  #######################################################################################################################################################
  ################################################################# return results ######################################################################
  #######################################################################################################################################################

    result_list<- list(latent_y=latent_y,
                       big_a=big_a,
                       big_z=big_z,
                       phi=phi,
                       beta=beta,
                       pai=pai,
                       hyper_record=hyper_record,
                       sigmay_record=sigmay_record,
                       sigmay_inv_record=sigmay_inv_record,
                       mc_num = mc_num,
                       individual_mean_in=individual_mean_in,
                       variance_g=variance_g,
                       mu_g = mu_g,
                       gene_person_mean_est = gene_person_mean_est,
                       ind_num = ind_num,
                       burn_in_prop = burn_in_prop,
                       thin_step = thin_step,
                       em_num = em_num,
                       p = p,
                       k = k,
                       n = n,
                       q = q,
                       prior_sparsity = prior_sparsity,
                       a_full = a_full,
                       train_index = train_index,
                       a_person = a_person,
                       obs_time_index =  obs_time_index,
                       obs_time_num = obs_time_num,
                       model_dgp = model_dgp)

  return(result_list)
}
