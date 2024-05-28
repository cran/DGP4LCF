## ----setup--------------------------------------------------------------------
library(DGP4LCF)

## ----eval = FALSE-------------------------------------------------------------
#  
#  # for reproducibility purpose
#  set.seed(456)
#  
#  mcem_parameter_setup_irregular_time_result<-
#    mcem_parameter_setup(p = 100, k = 4, n = 17, q = 8,
#                                        obs_time_num = sim_fcs_truth$obs_time_num,
#                                        obs_time_index = sim_fcs_truth$obs_time_index,
#                                        a_person = sim_fcs_truth$a_person,
#                                        col_person_index = sim_fcs_truth$col_person_index,
#                                        y_init = sim_fcs_init$y_init_irregular,
#                                        a_init = sim_fcs_init$a_init_2,
#                                        z_init = sim_fcs_init$z_init_2,
#                                        phi_init = sim_fcs_init$phi_init_irregular,
#                                        a_full = sim_fcs_truth$a_full,
#                                        train_index = (1:8),
#                                        x = sim_fcs_truth$observed_x_train_irregular)
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  mcem_algorithm_irregular_time_result<-
#    mcem_algorithm(ind_x = TRUE,
#                   x = sim_fcs_truth$observed_x_train_irregular,
#                   mcem_parameter_setup_result = mcem_parameter_setup_irregular_time_result)
#  

## ----fig1, fig.height = 4, fig.width=7----------------------------------------

old_par <- par(no.readonly = TRUE)

par(mfrow = c(2,4))
for (em_index in 1:sim_fcs_results_irregular_6_8$mcem_algorithm_irregular_time_result$index_used){
  mcem_cov_plot(sim_fcs_results_irregular_6_8$mcem_algorithm_irregular_time_result$sigmay_record[em_index,,], k = 4, q = 8, title = paste0("MCEM Iteration ", em_index))
}

mcem_cov_plot(sim_fcs_truth$gp_sigmay_truth, k = 4, q = 10, title = "Truth: Correlated Factors")

par(old_par)


## ----eval = FALSE-------------------------------------------------------------
#  
#  gibbs_after_mcem_diff_initials_irregular_time_result<-
#    gibbs_after_mcem_diff_initials(ind_x = TRUE,
#                                   tot_chain = 5,
#                                   mcem_parameter_setup_result = mcem_parameter_setup_irregular_time_result,
#                                   mcem_algorithm_result = mcem_algorithm_irregular_time_result)
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  tot_chain<- 5
#  
#  for (chain_index in 1:tot_chain){
#  
#    gibbs_after_mcem_algorithm(chain_index = chain_index,
#                                              mc_num = 10000,
#                                              burnin = 3000,
#                                              thin_step = 10 ,
#                                              pathname = "path",
#                                              pred_indicator = TRUE,
#                                              pred_time_index = (9:10),
#                                              x = sim_fcs_truth$observed_x_train_irregular,
#                                              gibbs_after_mcem_diff_initials_result = gibbs_after_mcem_diff_initials_irregular_time_result,
#                                              mcem_algorithm_result = mcem_algorithm_irregular_time_result,
#                                              mcem_parameter_setup_result =  mcem_parameter_setup_irregular_time_result)
#  }
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  constant_list<- list(num_time_test = 2,
#                       mc_num = 10000,
#                       thin_step = 10,
#                       burnin = 3000,
#                       pathname = "path",
#                       p = 100,
#                       k = 4,
#                       n = 17,
#                       q = 8,
#                       ind_x = TRUE,
#                       pred_indicator = TRUE)
#  
#  for (chain_index in 1:tot_chain){
#  
#    gibbs_after_mcem_load_chains_result<- gibbs_after_mcem_load_chains(chain_index = chain_index,
#                                                                       gibbs_after_mcem_algorithm_result = constant_list)
#  
#    save(gibbs_after_mcem_load_chains_result,
#         file = paste0("path/chain_", chain_index,"_result.RData"))
#  
#  }
#  
#  gibbs_after_mcem_combine_chains_irregular_time_result<- gibbs_after_mcem_combine_chains(tot_chain = 5,
#                                                                                          gibbs_after_mcem_algorithm_result = constant_list)
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#   numerics_summary_do_not_need_alignment_irregular_time_result<-
#    numerics_summary_do_not_need_alignment(pred_x_truth =  sim_fcs_truth$observed_x_pred_reformat,
#                                               pred_x_truth_indicator = TRUE,
#                                           gibbs_after_mcem_combine_chains_result =  gibbs_after_mcem_combine_chains_irregular_time_result)

## -----------------------------------------------------------------------------

pred_result_overview<- matrix(c(sim_fcs_results_irregular_6_8$numerics_summary_do_not_need_alignment_irregular_time_result$pred_x_result$mae_using_median_est,
                                sim_fcs_results_irregular_6_8$numerics_summary_do_not_need_alignment_irregular_time_result$pred_x_result$mean_width_interval,        
                                sim_fcs_results_irregular_6_8$numerics_summary_do_not_need_alignment_irregular_time_result$pred_x_result$proportion_of_within_interval_biomarkers), 
                                nrow = 3, ncol = 1)

rownames(pred_result_overview)<- c("Mean Absolute Error", "Mean Width of Intervals", "Proportion of Coverage")
colnames(pred_result_overview)<- "Value"

pred_result_overview


## ----eval = FALSE-------------------------------------------------------------
#  
#  numerics_summary_need_alignment_irregular_time_result<-
#    numerics_summary_need_alignment(gibbs_after_mcem_combine_chains_result =  gibbs_after_mcem_combine_chains_irregular_time_result)
#  

## -----------------------------------------------------------------------------

q<- 8
k<- 4
n<- 17

a_train<-  sim_fcs_truth$a_full[(1:q)]

h3n2_data<- list()

list_temp <- vector("list", k)
for (list_index in 1:k){
  list_temp[[list_index]]<- a_train
}
h3n2_data$input<- list_temp

fcs_sigma_y_init_truth_for_train_data<- GPFDA::mgpCovMat(Data=h3n2_data, hp=sim_fcs_truth$gp_hp_truth)

# rescale truth
d_matrix<- diag(sqrt(diag(fcs_sigma_y_init_truth_for_train_data)))
d_matrix_inv<- solve(d_matrix)

fcs_real_y_rescaled<- array(0, dim = c(q,k,n))

for (person_index in 1:n){
  fcs_real_y_rescaled[,,person_index]<- matrix(d_matrix_inv%*%as.numeric(sim_fcs_truth$real_y[1:q,,person_index]),
                                               nrow = q,
                                               ncol = k)
}

# compare
mae_irregular_6_8<-  mean(abs(sim_fcs_results_irregular_6_8$numerics_summary_need_alignment_irregular_time_result$reordered_summary$latent_y[(1:q),,] - fcs_real_y_rescaled)) 

mae_irregular_6_8

