#' @title Displaying significant factor loadings in the heatmap.
#'
#' @description This function is used to visualize results of estimates of factor loadings (in heatmaps).
#'
#' @param factor_loading_matrix A matrix of dimension (p, k), which stores results for factor loadings.
#' @param heatmap_title A character. Title for the heatmap.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#'
#' @return A heatmap presenting posterior median estimates of factor loadings.
#' @export
factor_loading_heatmap<- function(factor_loading_matrix, heatmap_title){

  ######################################### waiting to add appropriate code here #########################################
  # present results of l in heatmap
  p<- nrow(factor_loading_matrix)
  k<- ncol(factor_loading_matrix)

  gene_index_sorted<- rep(0, times = p)

  gene_index<- 0

  for (factor_index in 1:k){

    gene_value_sorted<- sort(abs(factor_loading_matrix[factor_loading_matrix[,factor_index]!=0,factor_index]),decreasing=TRUE)

    for (index_sort in 1:length(gene_value_sorted)){

      temp_index<- which(gene_value_sorted[index_sort]==abs(factor_loading_matrix[,factor_index]))

      if (temp_index %in% gene_index_sorted[1:gene_index]){
        next
      } else {
        gene_index<- gene_index + 1
        gene_index_sorted[gene_index]<- temp_index
      }
    }
  }

  gene_index_sorted_length<- sum(gene_index_sorted!=0)

  for (reorder_index in 1:p){
    if (reorder_index %in% gene_index_sorted){
      next
    } else {
      gene_index<- gene_index + 1
      gene_index_sorted[gene_index]<- reorder_index
    }
  }

  ### reorder
  factor_loading_matrix_reordered<- factor_loading_matrix[gene_index_sorted,]
  factor_loading_matrix_reordered<- factor_loading_matrix_reordered[1:gene_index_sorted_length,]

  rownames(factor_loading_matrix_reordered)[1:gene_index_sorted_length]<-
    paste("Gene", gene_index_sorted[1:gene_index_sorted_length], sep = " ")

  colnames(factor_loading_matrix_reordered)<-
    paste("Factor", 1:k, sep = " ")

  # both positive and negative
  if (range(factor_loading_matrix_reordered)[1]<0 & range(factor_loading_matrix_reordered)[2]>0){

    bk_1<- seq(min(factor_loading_matrix_reordered),-0.01, by=0.01)
    bk_2<- seq(0,max(factor_loading_matrix_reordered), length = 20)
    bk <- c(bk_1, bk_2)

    heatmap<- pheatmap::pheatmap(factor_loading_matrix_reordered,
                                 angle_col = "0",
                                 cluster_rows = FALSE,
                                 cluster_cols = FALSE,
                                 #show_rownames = FALSE,
                                 #show_colnames = FALSE,
                                 main = heatmap_title,
                                 color = c(colorRampPalette(colors = c("light blue","white"))(length(bk_1)),colorRampPalette(colors = c("white","orange","red"))(length(bk_2))),
                                 #color = colorRampPalette(colors = c("white","orange","red"))(length(bk_2)),
                                 display_numbers = ifelse(factor_loading_matrix_reordered!=0, round(factor_loading_matrix_reordered,2)," "),
                                 #breaks = bk_2)
                                 breaks = bk)

  } else if (range(factor_loading_matrix_reordered)[1]>=0){ # only positive

    bk_2<- seq(0,max(factor_loading_matrix_reordered), length = 20)

    heatmap<- pheatmap::pheatmap(factor_loading_matrix_reordered,
                                 angle_col = "0",
                                 cluster_rows = FALSE,
                                 cluster_cols = FALSE,
                                 main = heatmap_title,
                                 color = colorRampPalette(colors = c("white","orange","red"))(length(bk_2)),
                                 display_numbers = ifelse(factor_loading_matrix_reordered!=0, round(factor_loading_matrix_reordered,2)," "),
                                 breaks = bk_2)

  } else if (range(factor_loading_matrix_reordered)[2]<=0){ # only negative

    bk_1<- seq(min(factor_loading_matrix_reordered),0, by=0.01)

    heatmap<- pheatmap::pheatmap(factor_loading_matrix_reordered,
                                 angle_col = "0",
                                 cluster_rows = FALSE,
                                 cluster_cols = FALSE,
                                 main = heatmap_title,
                                 color = colorRampPalette(colors = c("light blue","white"))(length(bk_1)),
                                 display_numbers = ifelse(factor_loading_matrix_reordered!=0, round(factor_loading_matrix_reordered,2)," "),
                                 breaks = bk_1)

  }

  # only positive

  return(heatmap)

}

