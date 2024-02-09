#' Deconvolution results of single-cell versus spatial transcriptome
#'
#' A data frame for each cell type in each spot using the CARD deconvolution method
#'
#' @format A data frame with 428 rows and 15 variables
#'
"spot_celltype"

#' Clustering results of the spatial transcriptome
#'
#' A data frame of clustering result obtained by BayesSpace
#'
#' @format A data frame with 428 rows and 5 variables:
#' \describe{
#'  \item{cluster.init}{cluster.init}
#'  \item{col}{col}
#'  \item{row}{row}
#'  \item{spatial.cluster}{the clustering result}
#'  \item{sizeFactor}{sizeFactor}
#'  ...
#' }
#'
"spot_clusters"
