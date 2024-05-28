#' Generating a table listing all possible combinations of the binary variables for one gene.
#'
#' @param k A numeric scalar. Number of latent factors.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#' vignette("bsfadgp_irregular_data_example", package = "DGP4LCF")
#'
#' @return A table listing all possible combinations of the binary variables for one gene.
#' @export
table_generator<- function(k){
  if (k>=2){
    current_z<- 2^k
    current_table<- matrix(0, nrow=current_z,ncol=k)

    prev_z<- 2^(k-1)
    prev_table<- table_generator(k-1)

    current_table[,(1:(k-1))]<- rbind(prev_table,prev_table)

    current_table[(1:prev_z),k]<- rep(0,times=prev_z)

    current_table[((prev_z+1):current_z),k]<- rep(1,times=prev_z)
  } else if (k==1){
    current_table<- matrix(c(0,1),nrow=2,ncol=1)
  } else {
    message("please input an integer no less than 1")
  }

  return(current_table)
}
