# CellTopic main
#' @include SpaTopic_function.R

#' @title CellTopic
#' @description Through deconvolution and clustering information, the SpaTopic
#' obtains CellTopics representing the spatial domains.
#'
#'
#' @param spot_celltype A data frame of the result of deconvolution, row
#' is spot, col is celltype.
#' @param spot_clusters A data frame contains clustering information for spots,
#' row is spots.
#' @param cluster A character. Use the first column in spot_clusters as the
#' column name if \code{NULL}. Or provide your own column names that represent
#' clustering information in spot_clusters.
#' @param num_topics A integer or a vector of integer. The number of topics.
#' default is \code{NULL}, find the best number of topics by default. if a vector of
#' integer, find the best number of topics by the given range.
#' @param percent A numeric from \code{0} to \code{1}. The percent of topics. Default is \code{0.6}.
#' @param Binarization Logical indicating if to choose one topic for each
#' CellTopic. Default is \code{FALSE}.
#' @param meta.cell Logical indicating if return MetaTopic, which is a data frame of the
#' cluster result of CellTopic.
#' @param k A integer of how much MetaTopic to choose. If \code{meta.cell = TRUE}, \code{k}
#' is Requested to given. Default is \code{4}.
#' @param method See \code{\link[stats]{hclust}}.
#'
#' @return A list with three data frame and one vector. "CellTopic" is a data
#' frame which can be add to a Seurat object. The "domain_topic" is a data frame,
#' row is CellTopic and col is domain. The "celltype_topic" is a data frame, row is
#' celltype and col is CellTopic. The "Cell_topic" is a vector of which topic be chosen
#' in each CellTopic. If meta.cell = TRUE, the "MetaTopic" will be given, which
#' is a data frame of the cluster result of CellTopic.
#' @export
#'
#' @examples
#' options (warn = -1)
#' result_list <- CellTopic(spot_celltype,
#' spot_clusters,
#' cluster = "spatial.cluster",
#' num_topics = 13,
#' percent = 0.7,
#' Binarization = FALSE,
#' meta.cell = FALSE,
#' k = NULL)
#' head(result_list[["CellTopic"]])
#' print(result_list[["domain_topic"]])
#' print(result_list[["celltype_topic"]])
#' print(result_list[["Cell_topic"]])
CellTopic <- function(
    spot_celltype,
    spot_clusters,
    cluster = NULL,
    num_topics = NULL,
    percent = NULL,
    Binarization = FALSE,
    meta.cell = FALSE,
    k = NULL,
    method = "complete"
  ){
  domain_cellytpe <- ksTest(
    spot_celltype = spot_celltype,
    spot_clusters = spot_clusters,
    cluster = cluster
  )
  topic_list <- MatrixFactorization(domain_cellytpe = domain_cellytpe, num_topics = num_topics)
  dt_topic_data <- topic_list[["domain_topic"]]
  ct_topic_data <- topic_list[["celltype_topic"]]
  result_list <- FindCellTopic(
    dt_topic_data = dt_topic_data,
    ct_topic_data = ct_topic_data,
    spot_clusters = spot_clusters,
    cluster = cluster,
    percent = percent,
    Binarization = Binarization
  )
  if(meta.cell){
    ct_data <- result_list[["celltype_topic"]]
    meta_tree <- MetaTopic(ct_topic_data = ct_data, k = k, method = method)
    result_list[["MetaTopic"]] <- meta_tree
  }
  return(result_list)
}


