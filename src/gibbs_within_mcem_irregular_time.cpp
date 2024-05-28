// Use RcppArmadillo (a linear algebra library) because it provides a
// 3D cube/tensor object. It is not strictly necessary but many (most/)
// R installs that have Rcpp will also have RcppArmadillo

//#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iomanip>

using namespace Rcpp;
using namespace arma;
using std::endl;

// [[Rcpp::export]]
arma::mat gibbs_within_mcem_irregular_time(List& latent_y_in,
                   List& h3n2_response,
                   List& missing_list,
                   arma::vec& missing_num,
                   const bool ipt_x,
                   List& big_a_in,
                   List& big_z_in,
                   arma::vec& phi,
                   arma::vec& pai,
                   arma::vec& beta,
                   const int k, const int n, const int p, const int q,
                   const double c0, const double c1, const double d0,
                   const double d1, const double e0, const double f0,
                   const int start_iter,
                   const int iters,
                   const arma::mat& big_z_table,
                   const bool ind_x,
                   arma::mat individual_mean,
                   const arma::vec& mu_g,
                   arma::vec& variance_g,
                   const double c2,
                   const double d2,
                   const arma::vec& obs_time_num,
                   const List& obs_time_index,
                   const arma::mat& sigmay_inv_full,
                   const arma::vec& missing_time_num,
                   const List& missing_time_index,
                   const List& prod_covnew_covestinv,
                   const List& cov_cond_dist,
                   const List& sigmay_inv,
                   const int missing_person_num,
                   const arma::vec& missing_person_index,
                   const int full_person_num,
                   const arma::vec& full_person_index) {

  arma::mat retval(1, 1);

  arma::mat big_a(as<arma::mat>(big_a_in[start_iter - 1]));
  arma::mat big_z(as<arma::mat>(big_z_in[start_iter - 1]));

  NumericVector latent_y_tmp = latent_y_in[0];
  arma::cube latent_y(latent_y_tmp.begin(), q, k, n);

  /////////////////////////// note that everytime we use this, need to ensure time index is labelled!! ///////////////////////////////////

  // initialisation of 'h3n2_response_cube'
  // if ipt_x = FALSE: it will be constant
  // if ipt_x = TRUE : it will change across iterations

  arma::cube h3n2_response_cube(p, q, n);

  for (int person_index = 0; person_index < n; person_index ++){

    h3n2_response_cube.slice(person_index).cols(0, obs_time_num(person_index) - 1) = (as<arma::mat>(h3n2_response[person_index]));

  }

  for (int iter = start_iter; iter < iters; iter++) {

    Rcout << "this is which iteration of Gibbs sampling:" << iter <<"\n";

    // Calculate multivariate normal params and take random draws
    // ----------------------------------------------------------

    // Element-wise multiplication
    arma::mat l_temp = big_a % big_z; // p*k

    arma::mat fast_comp = (l_temp.each_col() % phi).t(); // k*p matrix

    // k x k, rather than kq x kq, consisting of all elements needed in sigma_pos
    arma::mat mini_sigma = fast_comp * l_temp;

    ///////////////////////////////////// for people with incomplete observations   /////////////////////////////////////

    if (missing_person_num > 0){

      for (int j = 0; j< missing_person_num; j++){

        // this person's index in the person list (under cpp system)

        int person_index = (missing_person_index[j]-1);

        arma::uvec obs_time_vector_temp(as<arma::uvec>(obs_time_index[person_index]));

        obs_time_vector_temp = obs_time_vector_temp - 1;

        arma::mat sigmay_inv_temp(as<arma::mat>(sigmay_inv[j]));

        arma::mat sigma_pos(sigmay_inv_temp); // initialize the sigma_pos matrix as the matrix sigmay_inv, which is dependent on the common sigmay_inv and the indicator matrix

        arma::mat sigma_posinv(size(sigmay_inv_temp));

        // Expand mini_sigma along the diagonal and add to sigmay_inv
        for (int x = 0; x < k; x++) {
          for (int y = 0; y < k; y++) {
            // Starting point in the output matrix
            int startx = x * obs_time_num(person_index), starty = y * obs_time_num(person_index);
            for (int z = 0; z < obs_time_num(person_index); z++)
              sigma_pos(startx + z, starty + z) += mini_sigma(x, y);

          }
        }

        // Handle ill-conditioned covariance matrices
        int tries = 0;
        // Equiv to sigma_pos = inv(sigma_pos)
        while (inv(sigma_posinv, sigma_pos) == false) {
          // inv() returns false if matrix appears to be singular
          sigma_pos.diag() += 1e-8;

          // Place a limit on the number of retries
          tries++;
          if (tries > 3) {
            Rcout << "Error: Iter " << iter << ": Inverting covariance fails after 3 tries\n";
            return retval;
          }
        }

        // Handle non-symmetry
        if(! sigma_posinv.is_symmetric(0.01)){
          Rcout << "Warning: sigma_posinv is not symmetric";
        }

        sigma_posinv= 0.5*(sigma_posinv + sigma_posinv.t());

        arma::vec mu_pos(k*obs_time_num(person_index));

        int pos = 0;

        // By iterating over q, we can multiply fast_comp by h3n2_response
        // without expanding or flattening

        for (int ik = 0; ik < k; ik++) {
          for (int iq = 0; iq < obs_time_num(person_index); iq++) {

            arma::vec h3n2vec = h3n2_response_cube.slice(person_index).col(iq) - individual_mean.col(person_index);

            mu_pos[pos++] = as_scalar(fast_comp.row(ik) * h3n2vec);

          }
        }

        mu_pos = sigma_posinv * mu_pos;

        arma::mat draw_obs = arma::mvnrnd(mu_pos, sigma_posinv);

        latent_y.slice(person_index).rows(obs_time_vector_temp) = reshape(draw_obs.col(0), obs_time_num(person_index), k);

        // Rcout << "this is which person:" << person_index <<"\n";

        ///////////////////////////////////////////// for factors without observed gene expressions ////////////////////////////////////////////////////////////////////////
        arma::mat prod_covnew_covestinv_temp(as<arma::mat>(prod_covnew_covestinv[j]));

        // Rcout << "prod_covnew_covestinv_temp:" << prod_covnew_covestinv_temp <<"\n";

        arma::mat cov_cond_dist_temp(as<arma::mat>(cov_cond_dist[j]));

        // Rcout << "cov_cond_dist_temp:" << cov_cond_dist_temp <<"\n";

        // calculate mean vector for the conditional distribution

        arma::vec mean_cond_dist(prod_covnew_covestinv_temp*draw_obs.col(0));

        // Rcout << "mean_cond_dist:" << mean_cond_dist <<"\n";

        arma::mat draw_miss = arma::mvnrnd(mean_cond_dist, cov_cond_dist_temp);

        // Rcout << "draw_miss:" << draw_miss <<"\n";

        arma::uvec missing_time_vector_temp(as<arma::uvec>(missing_time_index[j])); // uner r system

        // Rcout << "missing_time_vector_temp:" << missing_time_vector_temp <<"\n";

        missing_time_vector_temp = missing_time_vector_temp - 1; // under c++ system

        // Rcout << "this is missing_time_vector_temp:" << missing_time_vector_temp <<"\n";

        latent_y.slice(person_index).rows(missing_time_vector_temp) = reshape(draw_miss.col(0), missing_time_num(j), k);

        // Rcout << "this is latent_y.slice(person_index).rows(missing_time_vector_temp):" << latent_y.slice(person_index).rows(missing_time_vector_temp) <<"\n";

      }

    }


    if (full_person_num > 0 ){

      ///////////////////////////////////// for people with complete observations   /////////////////////////////////////
      arma::mat sigma_pos(sigmay_inv_full); // initialize the sigma_pos matrix as the matrix sigmay_inv

      arma::mat sigma_posinv(size(sigmay_inv_full));

      // Expand mini_sigma along the diagonal and add to sigmay_inv
      for (int x = 0; x < k; x++) {
        for (int y = 0; y < k; y++) {
          // Starting point in the output matrix
          int startx = x * q, starty = y * q;
          for (int z = 0; z < q; z++)
            sigma_pos(startx + z, starty + z) += mini_sigma(x, y);
        }
      }

      // Handle ill-conditioned covariance matrices
      int tries = 0;
      // Equiv to sigma_pos = inv(sigma_pos)
      while (inv(sigma_posinv, sigma_pos) == false) {
        // inv() returns false if matrix appears to be singular
        sigma_pos.diag() += 1e-8;

        // Place a limit on the number of retries
        tries++;
        if (tries > 3) {
          Rcout << "Error: Iter " << iter << ": Inverting covariance fails after 3 tries\n";
          return retval;
        }
      }

      // Handle non-symmetry
      if(! sigma_posinv.is_symmetric(0.01)){
        Rcout << "Warning: sigma_posinv is not symmetric";
      }

      sigma_posinv = 0.5*(sigma_posinv + sigma_posinv.t());

      for (int j = 0; j< full_person_num; j++) {

        arma::vec mu_pos(k*q);

        int pos = 0;

        int person_index = (full_person_index[j]-1);

        for (int ik = 0; ik < k; ik++) {

          for (int iq = 0; iq < q; iq++) {

              arma::vec h3n2vec = h3n2_response_cube.slice(person_index).col(iq) - individual_mean.col(person_index);

              mu_pos[pos++] = as_scalar(fast_comp.row(ik) * h3n2vec);

          }
        }

        mu_pos = sigma_posinv * mu_pos;

        arma::mat draw_obs = arma::mvnrnd(mu_pos, sigma_posinv);

        latent_y.slice(person_index) = reshape(draw_obs.col(0), q, k);

        // Rcout << "this is which person:" << person_index <<"\n";

      }
    }

    NumericVector latent_tmp2(Dimension(q, k, n));
    std::copy(latent_y.begin(), latent_y.end(), latent_tmp2.begin());
    latent_y_in[iter] = latent_tmp2;

    // update matrix z: updated version - block update for one row

    int num_z = pow(2,k);

    int final_index = num_z-1;

    IntegerVector all_index = seq(0,final_index); // result: 0,1,2,...2^k-1 = 31 when k=5

    arma::cube latent_y_pre_cal(q, k, n);

    for (int person_index = 0; person_index < n; person_index ++){

       arma::uvec obs_time_vector_temp(as<arma::uvec>(obs_time_index[person_index])); // this person's y's index in the full matrix (under R counting system)

       obs_time_vector_temp = obs_time_vector_temp - 1; // transform it to under c++ counting system

       latent_y_pre_cal.slice(person_index).rows(0, obs_time_num(person_index) - 1) = latent_y.slice(person_index).rows(obs_time_vector_temp); // n_i*k

    }

    for (int g = 0; g < p; g++){

      NumericVector ratio_z(num_z);

      NumericVector prob_z(num_z);

      arma::vec one_vector(k);

      for (int table_index = 1; table_index < num_z; table_index ++){

        arma::rowvec big_l_vector = big_z_table.row(table_index) % big_a.row(g);

        double first_part = 0;
        for (int iter_index = 0; iter_index < k; iter_index ++){
          first_part += big_z_table(table_index, iter_index) * (log(pai(iter_index) / (1 - pai(iter_index))));
        }

        double second_part = 0;


        for (int person_index = 0; person_index < n; person_index ++){

          for (int time_index = 0; time_index < obs_time_num(person_index); time_index ++){

            double sum_ly = sum(big_l_vector % (latent_y_pre_cal.slice(person_index)).row(time_index));

            second_part += ((2*((h3n2_response_cube.slice(person_index))(g,time_index) - individual_mean(g,person_index))*sum_ly) - pow(sum_ly, 2))*(phi(g)/2.0);

          }

        }

        ratio_z[table_index] = first_part + second_part;

      }

      int index_max_ratio_z = which_max(ratio_z);

      ratio_z = ratio_z - ratio_z[index_max_ratio_z];

      prob_z[index_max_ratio_z] = 1.0/sum(exp(ratio_z));

      prob_z =  prob_z[index_max_ratio_z] * exp(ratio_z);

      int sample_index = sample(all_index, 1, TRUE, prob_z)[0];

      big_z.row(g) = big_z_table.row(sample_index);

    }

    // Rcout << "this is big_z:" <<"\n";

    // Update matrix a
    // ---------------

    // This sum can be pre-calculated

    arma::mat sum_1(k, k, fill::zeros);

    for (int l = 0; l < n; l++) {

      // arma::uvec obs_time_vector_temp(as<arma::uvec>(obs_time_index[l])); // this person's y's index in the full matrix (under R counting system)

      // obs_time_vector_temp = obs_time_vector_temp - 1; // transform it to under c++ counting system

      // arma::mat tmp(latent_y.slice(l).rows(obs_time_vector_temp)); // n_i*k

      arma::mat tmp(latent_y_pre_cal.slice(l).rows(0, obs_time_num(l) - 1));

      sum_1 += tmp.t() * tmp; // k*k matrix
    }

    for (int g = 0; g < p; g++) {

      arma::rowvec sum_2(k, fill::zeros); // 1*k

      for (int l = 0; l < n; l++) {

        arma::mat tmp(latent_y_pre_cal.slice(l).rows(0, obs_time_num(l) - 1)); // n_i * k

        arma::rowvec tmp_individual_mean(obs_time_num(l), fill::value(individual_mean(g,l))); // 1*n_i

        sum_2 += (h3n2_response_cube.slice(l)(g, span(0, obs_time_num(l) - 1)) - tmp_individual_mean) * tmp;

      }

      arma::mat sum_1_g = sum_1 * phi(g);

      // sigma_ag<- solve(diag(big_z[i,g,])%*%sum_1%*%diag(big_z[i,g,])+diag(beta[(i-1),]))
      // 5-vector
      arma::mat zdiag(k, k, fill::zeros), betadiag(k, k, fill::zeros);
      zdiag.diag() = big_z.row(g);
      betadiag.diag() = beta;

      arma::mat sigma_ag = zdiag * sum_1_g * zdiag + betadiag;
      arma::mat sigma_aginv(size(sigma_ag));

      // Handle ill-conditioned covariance matrices
      int tries = 0;
      // Equiv to sigma_ag = inv(sigma_ag)
      while (inv(sigma_aginv, sigma_ag) == false) {
        // inv() returns false if matrix appears to be signular
        sigma_ag.diag() += 1e-8;

        // Place a limit on the number of retries
        tries++;
        if (tries >= 3) {
          Rcout << "Error: Iter " << iter << ": Inverting covariance fails after 3 tries\n";
          return retval;
        }
      }

      // Handle non-symmetry

      if(! sigma_aginv.is_symmetric(0.01)){
        Rcout << "Warning: the matrix is not symmetric";
      }

      sigma_aginv= 0.5*(sigma_aginv + sigma_aginv.t());

      sum_2 = sum_2 * phi(g);

      arma::vec mu_ag = sigma_aginv * zdiag * sum_2.t();

      arma::mat draw = arma::mvnrnd(mu_ag, sigma_aginv);

      big_a.row(g) = trans(draw.col(0));


    }

    // Rcout << "this is  big_a: " <<  big_a <<"\n"; // kq*1

    // Remaining params
    // ----------------

    // Update pai
    for (int j = 0; j < k; j++) {

      double sum_j = sum(big_z.col(j));

      double e = e0 + sum_j;

      double f = f0 + (p - sum_j);

      pai(j) = R::rbeta(e, f);
    }

    // Rcout << "this is pai:" <<"\n";

    // Update beta
    double beta_c = c1 + (0.5*p);

    for (int j = 0; j < k; j++) {

      double beta_d = d1 + (sum(big_a.col(j) % big_a.col(j))/2.0);

      beta(j) = R::rgamma(beta_c, 1.0/beta_d);

    }

    // Rcout << "this is beta:" <<"\n";

    // Update phi

    double phi_c = c0 + (0.5 * sum(obs_time_num));

    for (int g = 0; g < p; g++) {

      double sum_g = 0;

      for (int l = 0; l < n; l++) {

        arma::mat latent_slice(latent_y_pre_cal.slice(l).rows(0, obs_time_num(l) - 1));

        arma::rowvec tosum = (big_a.row(g) % big_z.row(g)) * latent_slice.t(); // 1*n_i

        arma::rowvec tmp_individual_mean(obs_time_num(l), fill::value(individual_mean(g,l))); // 1*n_i

        tosum = h3n2_response_cube.slice(l)(g, span(0, obs_time_num(l) - 1)) - tmp_individual_mean - tosum;

        sum_g += sum(tosum % tosum);
      }

      double phi_d = d0 + sum_g / 2.0;

      phi(g) = R::rgamma(phi_c, 1.0/phi_d);
    } // phi is (inverse) of (variance)

    // Rcout << "this is phi:" <<"\n";

    // update v_{g,i}, which is a p*n matrix, therefore can be treated like big_a and big_z

    if (ind_x){

      for (int g = 0; g<p; g++){

        for (int l = 0; l < n; l++) {

          double variance_g_pos_temp = 1.0/((1.0/variance_g(g)) + (obs_time_num(l)*phi(g)));

          arma::mat latent_slice(latent_y_pre_cal.slice(l).rows(0, obs_time_num(l) - 1));

          arma::rowvec tosum = (big_a.row(g) % big_z.row(g)) * latent_slice.t(); // 1*n_i

          tosum = h3n2_response_cube.slice(l)(g, span(0, obs_time_num(l) - 1)) - tosum;

          double mu_g_pos_temp = ((mu_g(g)/variance_g(g)) + (sum(tosum)*phi(g)))*variance_g_pos_temp;

          individual_mean(g,l) = as<NumericVector>(rnorm(1, mu_g_pos_temp, sqrt(variance_g_pos_temp)))[0]; // individual_mean always saves the updated values


        }

        double sigma2_c = c2 + (0.5 * n);


        double sum_g = 0;

        for (int l = 0; l < n; l++) {
          sum_g += pow((individual_mean(g,l) - mu_g(g)),2);
        }

        double sigma2_d = d2 + (0.5*sum_g);

        variance_g(g) = 1.0/R::rgamma(sigma2_c, 1.0/sigma2_d);


      }
    }

    // certain positions (i.e., imputed values) in 'h3n2_response_cube' needs to be updated
    if (ipt_x){

      arma::mat l_temp = big_a % big_z;

      for (int person_index = 0; person_index < n; person_index ++){

         // evaluate if need to be imputed
         if (missing_num(person_index)>0){

           for (int missing_index = 0; missing_index < missing_num(person_index); missing_index ++){

             arma::mat missing_list_temp(as<arma::mat>(missing_list[person_index]));

             // obtain the position of the NA value

             int row_index = (missing_list_temp(missing_index, 0) - 1); // gene_index
             int column_index = (missing_list_temp(missing_index, 1) - 1); // time_index

             // mean = l %*% y
             // latent_y: q,k,n
             double temp_mean = sum(l_temp.row(row_index) % latent_y.slice(person_index).row(column_index)) +  individual_mean(row_index, person_index); // (1*k) % (1*k)

             // impute using the subject-gene mean
             h3n2_response_cube(row_index, column_index, person_index) = as<NumericVector>(rnorm(1, temp_mean, sqrt(1/phi(row_index))))[0];

           }

         } else {
           continue;
         }

         // check if the ipt_x is doing its job - pass the test
         //if (person_index == 65){
         //   Rcout << "this is h3n2_response_cube for this person:" << h3n2_response_cube.slice(person_index).col(7) <<"\n";
         //}
      }

    }

    // Copy matrices to output lists
    big_a_in[iter] = wrap(big_a); //converting armadillo objects to rcpp objects
    big_z_in[iter] = wrap(big_z);

  }

  // return retval: a dummy matrix to return;
  return retval;
}

