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

// self-define a function that can: save vec results into a .csv file - for phi and pai
void print_irregular(std::ofstream& file, arma::vec& vec) {
  // n_elem
  //Rcpp::Rcout << "printing vector\n";
  for (int i = 0; i < vec.n_elem; i++)
    file << std::setprecision(9) << std::fixed << vec(i) << ","; // std::fixed makes sure that c++ won't ignore ending 0s and std::setprecision(9) makes sure the saved result is accurate to the 9th digit
  file << std::endl;
}

// self-define a function that can: save matrix results into a .csv file - for big_a and big_z
void print_irregular(std::ofstream& file, arma::mat& mat) {
  //Rcout << "printing matrix\n";
  for (int i = 0; i < mat.n_rows; i++) // write in rows rather than write in columns
    for (int j = 0; j < mat.n_cols; j++){
      // Rcout << row << "," << col;
      file << std::setprecision(9) << std::fixed << mat(i, j) << ",";
    }
  file << std::endl;
}

// self-define a function that can: save cube results into a .csv file - for latent_y
// cube: q*k*n
void print_irregular(std::ofstream& file, arma::cube& cube) {
  // n_rows, n_cols, n_slices
  //Rcout << "printing cubes\n";
  for (int slice = 0; slice < cube.n_slices; slice++) // first fix person
    for (int col = 0; col < cube.n_cols; col++) // then column
      for (int row = 0; row < cube.n_rows; row++) // then row
        file << std::setprecision(9) << std::fixed << cube(row, col, slice) << ",";
  file << std::endl;
}

// [[Rcpp::export]]
arma::mat gibbs_after_mcem_irregular_time(arma::cube latent_y,
                   List& h3n2_response,
                   List& missing_list,
                   arma::vec& missing_num,
                   const bool ipt_x,
                   arma::mat& big_a,
                   arma::mat& big_z,
                   arma::vec& phi,
                   arma::vec& pai,
                   arma::vec& beta,
                   const int k, const int n, const int p, const int q, const int num_time_test,
                   const double c0, const double c1, const double d0,
                   const double d1, const double e0, const double f0,
                   const int start_iter,
                   const int iters,
                   const bool pred_indicator,
                   const bool noise,
                   arma::cube pred_y,
                   arma::cube pred_x,
                   const arma::mat& big_z_table,
                   const bool ind_x,
                   arma::mat individual_mean,
                   const arma::vec& mu_g,
                   arma::vec& variance_g,
                   const double c2,
                   const double d2,
                   const arma::vec& obs_time_num,
                   const List& obs_time_index,
                   const arma::vec& missing_time_num,
                   const List& missing_time_index,
                   const List& prod_covnew_covestinv,
                   const List& cov_cond_dist,
                   const List& sigmay_inv,
                   const int thin_step,
                   const int burnin,
                   const arma::mat& sigmay_inv_full,
                   const int missing_person_num,
                   const arma::vec& missing_person_index,
                   const int full_person_num,
                   const arma::vec& full_person_index,
                   String directory_name) {

  // declare this at the beginning - but only open the relevant .csv file under 'prediction'
  std::ofstream f_phi;
  std::ofstream f_pai;
  std::ofstream f_beta;
  std::ofstream f_big_a;
  std::ofstream f_big_z;
  std::ofstream f_latent_y;
  std::ofstream f_pred_y;
  std::ofstream f_pred_x;
  std::ofstream f_individual_mean;
  std::ofstream f_variance_g;

  // set printing precision at the beginning of the file
  Rcout.precision(10);
  
  f_phi.open("phi.csv");
  f_pai.open("pai.csv");
  f_beta.open("beta.csv");
  f_big_a.open("big_a.csv");
  f_big_z.open("big_z.csv");
  
  f_latent_y.open("latent_y.csv");
  
  if (ind_x){
    f_individual_mean.open("individual_mean.csv");
    f_variance_g.open("variance_g.csv");
  }
  
  if (pred_indicator){
    f_pred_y.open("pred_y.csv");
    f_pred_x.open("pred_x.csv");
  }

    // String dir_phi = directory_name;
    // f_phi.open(dir_phi+="phi.csv");
    // 
    // String dir_pai = directory_name;
    // f_pai.open(dir_pai+="pai.csv");
    // 
    // String dir_beta = directory_name;
    // f_beta.open(dir_beta+="beta.csv");
    // 
    // String dir_a = directory_name;
    // f_big_a.open(dir_a+= "big_a.csv");
    // 
    // String dir_z = directory_name;
    // f_big_z.open(dir_z+="big_z.csv");
    // 
    // String dir_y = directory_name;
    // f_latent_y.open(dir_y+="latent_y.csv");
    // 
    // String dir_individual_mean = directory_name;
    // f_individual_mean.open(dir_individual_mean+="individual_mean.csv");
    // 
    // String dir_variance_g = directory_name;
    // f_variance_g.open(dir_variance_g+="variance_g.csv");
    // 
    // String dir_pred_y = directory_name;
    // f_pred_y.open(dir_pred_y+="pred_y.csv");
    // 
    // String dir_pred_x = directory_name;
    // f_pred_x.open(dir_pred_x+="pred_x.csv");

  arma::mat retval(1, 1);

    // note that everytime we use this, need to ensure time index is labelled!!

    arma::cube h3n2_response_cube(p, q, n);

    for (int person_index = 0; person_index < n; person_index ++){

      h3n2_response_cube.slice(person_index).cols(0, obs_time_num(person_index) - 1) = (as<arma::mat>(h3n2_response[person_index]));

    }

  for (int iter = start_iter; iter < iters; iter++) {

    // Rcout << "this is which iteration of Gibbs sampling:" << iter <<"\n";

    // Calculate multivariate normal params and take random draws
    // ----------------------------------------------------------

    // Element-wise multiplication
    arma::mat l_temp = big_a % big_z; // p*k

    arma::mat fast_comp = (l_temp.each_col() % phi).t(); // k*p matrix


    // k x k, rather than kq x kq, consisting of all elements needed in sigma_pos
    arma::mat mini_sigma = fast_comp * l_temp;


    if (pred_indicator){

      for (int j = 0; j< n; j++){

        arma::mat sigmay_inv_temp(as<arma::mat>(sigmay_inv[j]));

        arma::mat sigma_pos(sigmay_inv_temp); // initialize the sigma_pos matrix as the matrix sigmay_inv, which is dependent on the common sigmay_inv and the indicator matrix

        arma::mat sigma_posinv(size(sigmay_inv_temp));

        // Expand mini_sigma along the diagonal and add to sigmay_inv
        for (int x = 0; x < k; x++) {
          for (int y = 0; y < k; y++) {
            // Starting point in the output matrix
            int startx = x * obs_time_num(j), starty = y * obs_time_num(j);
            for (int z = 0; z < obs_time_num(j); z++)
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

        arma::mat lower_triangular = chol(sigma_posinv, "lower"); // lower triangular matrix

        // JCC: Random draws from independent MVN
        arma::vec mu_vec(sigma_posinv.n_rows, fill::zeros); // create a zero-vector, length = the dimension of sigma_posinv



        arma::mat cov_matrix(size(sigma_posinv), fill::eye); // create a diagonal matrix with all elements =1 in the main diagonal, dim = the dimension of sigma_posinv



        arma::mat draw_matrix = arma::mvnrnd(mu_vec, cov_matrix);



        // JCC: calculate the common part
        arma::vec common_part(lower_triangular*draw_matrix.col(0));



        // Calculate fast_comp * h3n2_response[j]

        arma::vec mu_pos(k*obs_time_num(j)); // 24*1



        int pos = 0;

        // By iterating over q, we can multiply fast_comp by h3n2_response
        // without expanding or flattening
        for (int ik = 0; ik < k; ik++) {
          for (int iq = 0; iq < obs_time_num(j); iq++) {

            arma::vec h3n2vec = h3n2_response_cube.slice(j).col(iq) - individual_mean.col(j);

            mu_pos[pos++] = as_scalar(fast_comp.row(ik) * h3n2vec);


          }
        }

        mu_pos = sigma_posinv * mu_pos;

        arma::uvec obs_time_vector_temp(as<arma::uvec>(obs_time_index[j]));

        obs_time_vector_temp = obs_time_vector_temp - 1;

        latent_y.slice(j).rows(obs_time_vector_temp) = reshape(mu_pos + common_part, obs_time_num(j), k);

        //////////////////////// for factors without observed gene expressions: including both missing and predicted y ////////////////////////////////////////////////////////////////////////
        arma::mat prod_covnew_covestinv_temp(as<arma::mat>(prod_covnew_covestinv[j]));



        arma::mat cov_cond_dist_temp(as<arma::mat>(cov_cond_dist[j]));



        // calculate mean vector for the conditional distribution
        arma::vec mean_cond_dist(prod_covnew_covestinv_temp*(mu_pos + common_part));



        arma::mat draw = arma::mvnrnd(mean_cond_dist, cov_cond_dist_temp);

        arma::uvec missing_time_vector_temp(as<arma::uvec>(missing_time_index[j]));

        missing_time_vector_temp = missing_time_vector_temp - 1;

        latent_y.slice(j).rows(missing_time_vector_temp) = reshape(draw.col(0), missing_time_num(j), k);

        //////////////////////////////// save predicted y to pred_y ///////////////////////////////////////////////////////////////////////

        pred_y.slice(j) = latent_y.slice(j).rows(q, (q+num_time_test-1)); // num_time_test*k matrix

        //////////////////////////////////// for factors under prediction: no people have observed gene expressions for this /////////////////

      }
    } else {


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

          ///////////////////////////////////////////// for factors without observed gene expressions ////////////////////////////////////////////////////////////////////////
          arma::mat prod_covnew_covestinv_temp(as<arma::mat>(prod_covnew_covestinv[j]));

          arma::mat cov_cond_dist_temp(as<arma::mat>(cov_cond_dist[j]));

          // calculate mean vector for the conditional distribution

          arma::vec mean_cond_dist(prod_covnew_covestinv_temp*draw_obs.col(0));

          arma::mat draw_miss = arma::mvnrnd(mean_cond_dist, cov_cond_dist_temp);

          arma::uvec missing_time_vector_temp(as<arma::uvec>(missing_time_index[j])); // uner r system

          missing_time_vector_temp = missing_time_vector_temp - 1; // under c++ system

          latent_y.slice(person_index).rows(missing_time_vector_temp) = reshape(draw_miss.col(0), missing_time_num(j), k);

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

        }
      }

    }


    // update matrix z: updated version - block update for one row
    int num_z = pow(2,k);

    int final_index = num_z-1;

    IntegerVector all_index = seq(0,final_index);

    arma::cube latent_y_pre_cal(q, k, n);

    for (int person_index = 0; person_index < n; person_index ++){

      arma::uvec obs_time_vector_temp(as<arma::uvec>(obs_time_index[person_index])); // this person's y's index in the full matrix (under R counting system)

      obs_time_vector_temp = obs_time_vector_temp - 1; // transform it to under c++ counting system

      latent_y_pre_cal.slice(person_index).rows(0, obs_time_num(person_index) - 1) = latent_y.slice(person_index).rows(obs_time_vector_temp); // n_i*k

    }

    for (int g = 0; g < p; g++){

      // Rcout << "this is for updating big_z for which gene: " << g <<"\n";

      NumericVector ratio_z(num_z);
      NumericVector prob_z(num_z);

      arma::vec one_vector(k);

      for (int table_index = 1; table_index < num_z; table_index ++){

        arma::rowvec big_l_vector = big_z_table.row(table_index) % big_a.row(g);

        // double first_part = (big_z_table.row(table_index) * log(pai.col(iter-1))) - (big_z_table.row(table_index) * log( one_vector - pai.col(iter-1)));
        double first_part = 0;
        for (int iter_index = 0; iter_index < k; iter_index ++){
          first_part += big_z_table(table_index, iter_index) * (log(pai(iter_index) / (1 - pai(iter_index))));
        }

        double second_part = 0;

        for (int person_index = 0; person_index < n; person_index ++){

          for (int time_index = 0; time_index < obs_time_num(person_index); time_index ++){


            double sum_ly = sum(big_l_vector % (latent_y_pre_cal.slice(person_index)).row(time_index));

              second_part += ((2*((h3n2_response_cube.slice(person_index))(g,time_index)-individual_mean(g,person_index))*sum_ly) - pow(sum_ly,2))*(phi(g)/2.0);

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

    // Update matrix a
    // ---------------

    // This sum can be pre-calculated
    arma::mat sum_1(k, k, fill::zeros);

    for (int l = 0; l < n; l++) {

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

      arma::mat lower_triangular = chol(sigma_aginv, "lower");

      // JCC: Random draws from independent MVN
      arma::vec mu_vec(sigma_aginv.n_rows, fill::zeros); // create a zero-vector, length = the dimension of sigma_posinv
      arma::mat cov_matrix(size(sigma_aginv), fill::eye); // create a diagonal matrix with all elements =1 in the main diagonal, dim = the dimension of sigma_posinv


      arma::mat draw = arma::mvnrnd(mu_vec, cov_matrix);

      sum_2 = sum_2 * phi(g);
      arma::vec mu_ag = sigma_aginv * zdiag * sum_2.t();

      big_a.row(g) = trans(mu_ag + (lower_triangular*draw.col(0)));

    }


    // Remaining params
    // ----------------

    // Update pai
    for (int j = 0; j < k; j++) {

      double sum_j = sum(big_z.col(j));

      double e = e0 + sum_j;

      double f = f0 + (p - sum_j);


      pai(j) = R::rbeta(e, f);
    }



    // Update beta
    double beta_c = c1 + (0.5*p);

    for (int j = 0; j < k; j++) {

      double beta_d = d1 + (sum(big_a.col(j) % big_a.col(j))/2.0);

      beta(j) = R::rgamma(beta_c, 1.0/beta_d);


    }



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



    // update v_{g,i}, which is a p*n matrix, therefore can be treated like big_a and big_z

    if (ind_x){

      double sigma2_c = c2 + (0.5 * n);

      for (int g = 0; g< p; g++){

        // mu_g_pos is specific to each individual

        for (int l = 0; l < n; l++) {

          // variance_g_pos is specific to each individual under irregular time points

          double variance_g_pos_temp = 1.0/((1.0/variance_g(g)) + (obs_time_num(l)*phi(g)));

          arma::mat latent_slice(latent_y_pre_cal.slice(l).rows(0, obs_time_num(l) - 1));

          arma::rowvec tosum = (big_a.row(g) % big_z.row(g)) * latent_slice.t(); // 1*n_i

          tosum = h3n2_response_cube.slice(l)(g, span(0, obs_time_num(l) - 1)) - tosum;

          double mu_g_pos_temp = ((mu_g(g)/variance_g(g)) + (sum(tosum)*phi(g)))*variance_g_pos_temp;

          individual_mean(g,l) = as<NumericVector>(rnorm(1, mu_g_pos_temp, sqrt(variance_g_pos_temp)))[0]; // individual_mean always saves the updated values

        }

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
      }

    }

    //////////////////////////////////////////////////// prediction for x using all updates ////////////////////////////////////////////////////////

    if (pred_indicator){

      if (noise) {

        arma::mat l_temp_updated(big_z % big_a);

        for (int person_index=0; person_index < n; person_index++){

          arma::mat mean_final_matrix(l_temp_updated*pred_y.slice(person_index).t()); // expect this to be a num_time_test*k dimensional matrix

          for (int biomarker_index=0; biomarker_index < p; biomarker_index ++ ){

            arma::rowvec tmp_individual_mean(num_time_test, fill::value(individual_mean(biomarker_index, person_index)));

            arma::mat sigmax(num_time_test,num_time_test,fill::zeros);

            arma::vec vec_temp(num_time_test, fill::ones); // a vector of 1

            vec_temp *= 1/phi(biomarker_index); // a vector of 1/phi(biomarker_index)

            sigmax.diag() = vec_temp; // variance of the normal distribution

            arma::rowvec mean_final_vec_row = mean_final_matrix.row(biomarker_index) + tmp_individual_mean;

            arma::vec mean_final_vec = conv_to<colvec>::from(mean_final_vec_row);

            arma::mat draw = arma::mvnrnd(mean_final_vec, sigmax);

            // the assigned object is a Rcpp column vector; therefore, we need the transformation
            pred_x.slice(person_index).col(biomarker_index) = draw.col(0);
          }
        }

        // Rcout << "this is pred_x:" << pred_x <<"\n";

      } else {

        arma::mat l_temp_updated(big_z % big_a);

        for (int person_index=0; person_index < n; person_index++){

          arma::mat mean_final_matrix(l_temp_updated*pred_y.slice(person_index).t()); // expect this to be a num_time_test*k dimensional matrix

          for (int biomarker_index=0; biomarker_index < p; biomarker_index ++ ){

            pred_x.slice(person_index).col(biomarker_index) = conv_to<colvec>::from(mean_final_matrix.row(biomarker_index));

            // for ind_x = 1, will add; for ind_x = 0, will not change
            arma::rowvec tmp_individual_mean(num_time_test, fill::value(individual_mean(biomarker_index, person_index)));

            arma::vec tmp_individual_mean_col = conv_to<colvec>::from(tmp_individual_mean);

            pred_x.slice(person_index).col(biomarker_index) = pred_x.slice(person_index).col(biomarker_index) + tmp_individual_mean_col;
          }
        }

      }

    }


       ////////////////////////////////////////////////////  save results ////////////////////////////////////////////////////////
      //Rcout << "iter " << iter << endl;
      if (iter % thin_step == 0 && iter > burnin){
        //Rcout << "printing...\n";
        // for vectors
        print_irregular(f_phi, phi);
        print_irregular(f_pai, pai);
        print_irregular(f_beta, beta);

        // for matrices
        print_irregular(f_big_a, big_a);
        print_irregular(f_big_z, big_z);

        // for cubes
        print_irregular(f_latent_y, latent_y);

        //if (pred_indicator){

          print_irregular(f_pred_y, pred_y);

          print_irregular(f_pred_x, pred_x);

       // }

        //if (ind_x){

          print_irregular(f_individual_mean, individual_mean);

          print_irregular(f_variance_g, variance_g);

        //}

      }


  }

  // return retval: a dummy matrix to return;
  return retval;
}

