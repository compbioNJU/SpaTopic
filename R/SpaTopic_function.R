# CellTopic Functions

character()
#' @title ksTest
#' @description
#' KS test is performed on the deconvolution results and the overall level of
#' each cluster.
#'
#' @param spot_celltype A data frame of the result of deconvolution, row
#' is spot, col is celltype.
#' @param spot_clusters A data frame contains clustering information for spots,
#' row is spots.
#' @param cluster A character. Use the first column in spot_clusters as the
#' column name if \code{NULL}. Or provide your own column names that represent
#' clustering information in spot_clusters.
#'
#' @return A data frame, row is spot_domain, col is celltype score.
#' @importFrom stats ks.test
#' @export
#'
ksTest <- function(
    spot_celltype,
    spot_clusters,
    cluster = NULL
    ) {
  if (!is.data.frame(spot_clusters)) {
    spot_clusters <- as.data.frame(spot_clusters)
  }
  if (!is.data.frame(spot_celltype)) {
    spot_celltype <- as.data.frame(spot_celltype)
  }
  if (nrow(spot_clusters) != nrow(spot_celltype)) {
    stop("Data frame dimensions are not uniform!")
  }
  if (is.null(cluster)) {
    cluster <- colnames(spot_clusters)[1]
  } else if (!cluster %in% colnames(spot_clusters)) {
    cat("You can use these clusters:\n", colnames(spot_clusters), "\n")
    stop("The selected cluster is not in the data!")
  }
  cluster_sort <- sort(unique(spot_clusters[[cluster]]))
  celltype <- colnames(spot_celltype)
  data <- merge(spot_celltype, spot_clusters, by = 0, all.x = TRUE, )
  rownames(data) <- data[["Row.names"]]
  data[["Row.names"]] <- NULL
  domain_cellytpe <- data.frame(row.names = cluster_sort)
  # Do KS test of each celltype and each cluster
  for (i in celltype) {
    celltype_ks <- c()
    domain_all <- data[, i]
    for (j in cluster_sort){
      domain <- data[which(data[[cluster]] == j), i]
      ks <- stats::ks.test(domain, domain_all, alternative = "two.sided")$statistic
      ks <- as.numeric(ks)
      if (mean(domain) - mean(domain_all) < 0) {
        # If the average is low, a negative value should be taken
        ks <- -ks
      }
      celltype_ks <- c(celltype_ks, ks)
    }
    domain_cellytpe <- cbind(domain_cellytpe, celltype_ks)
  }
  colnames(domain_cellytpe) <- celltype
  rownames(domain_cellytpe) <- paste0("spot_domain_", rownames(domain_cellytpe))
  return(domain_cellytpe)
}


#' @title MatrixFactorization
#' @description
#' Perform matrix factorization of the KS test results, see the matrix
#' factorization method at \code{\link[topicmodels]{LDA}}.
#'
#' @param domain_cellytpe A data frame. The result of function \code{\link[SpaTopic]{ksTest}}.
#' @param num_topics A integer or a vector of integer. The number of topics. To
#' find the best number of topics if \code{NULL}. Or to find the best number of
#' topics by the given range if a vector of integer. Or provide your own number
#' if a integer.
#'
#' @return A list with two data frame. The "domain_topic" is a data frame,
#' row is topic and col is domain. The "celltype_topic" is a data frame, row is
#' celltype and col is topic.
#' @importFrom slam as.simple_triplet_matrix
#' @importFrom topicmodels LDA
#' @importFrom modeltools posterior
#' @export
#'
MatrixFactorization <- function(domain_cellytpe, num_topics = NULL) {
  if (!is.data.frame(domain_cellytpe)) {
    domain_cellytpe <- as.data.frame(domain_cellytpe)
  }
  # Convert to an integer
  domain_cellytpe_integer <- (domain_cellytpe - min(domain_cellytpe)) / (max(domain_cellytpe) - min(domain_cellytpe)) * 100
  domain_cellytpe_integer <- round(domain_cellytpe_integer, digits = 0)
  # Generate corpus
  corpus <- slam::as.simple_triplet_matrix(as.matrix(domain_cellytpe_integer))
  if (is.null(num_topics)) {
    num_topics <- tryTopic(corpus)
  }else if(length(num_topics) > 1){
    num_topics <- tryTopic(corpus, num_topics)
  }
  # Run the LDA model
  lda_model <- topicmodels::LDA(corpus, k = num_topics, control = list(seed = 1234))
  # domain topic
  dt_topic_data <- modeltools::posterior(lda_model)$topics
  colnames(dt_topic_data) <- paste0("topic", 1:ncol(dt_topic_data))
  dt_topic_data <- as.data.frame(t(dt_topic_data))
  # celltype topic
  ct_topic_data <- modeltools::posterior(lda_model)$terms
  colnames(ct_topic_data) <- colnames(spot_celltype)
  row.names(ct_topic_data) <- paste0("topic", 1:nrow(ct_topic_data))
  ct_topic_data <- as.data.frame(t(ct_topic_data))
  # result list
  topic_list <- list(
    domain_topic = dt_topic_data,
    celltype_topic = ct_topic_data
  )
  return(topic_list)
}


#' @title TryTopic
#' @description
#' Probably find the best number of topics.
#'
#' @param corpus See \code{\link[topicmodels]{LDA}}.
#' @param num_topics_to_try A vector of numbers. Default is \code{c(7:20)}.
#'
#' @return Probably the best number of topics.
#'
#' @noRd
TryTopic <- function(corpus, num_topics_to_try = NULL) {
  if (is.null(num_topics_to_try)) {
    num_topics_to_try <- 7:20
  }
  perplexities <- sapply(num_topics_to_try, perplexity_calculation, corpus)
  best_num_topics <- num_topics_to_try[which.min(perplexities)]
  return(best_num_topics)
}


#' Calculate the perplexity of the specified number of topics.
#'
#' @param num_topics The number of topics.
#' @param corpus See \code{\link[topicmodels]{LDA}}.
#'
#' @return The inverse of the model perplexity.
#' @importFrom topicmodels LDA
#' @importFrom stats logLik
#'
#' @noRd
Perplexity_calculation <- function(num_topics, corpus) {
  lda_model <- topicmodels::LDA(corpus, k = num_topics)
  perplexity <- stats::logLik(lda_model)[1]
  return(-perplexity)
}




#' @title FindCellTopic
#' @description
#' Based on the results of the \code{\link[SpaTopic]{MatrixFactorization}} and
#' the clustering information, the CellTopics representing each domain are
#' calculated.
#'
#' @param dt_topic_data A data frame, row is topic and col is domain. The result
#' of function \code{\link[SpaTopic]{MatrixFactorization}}.
#' @param ct_topic_data A data frame, row is celltype and col is topic. The
#' result of function \code{\link[SpaTopic]{MatrixFactorization}}.
#' @param spot_clusters A data frame contains clustering information for spots,
#' row is spots.
#' @param cluster A character. Use the first column in spot_clusters as the
#' column name if \code{NULL}. Or provide your own column names that represent
#' clustering information in spot_clusters.
#' @param percent A numeric from \code{0} to \code{1}. The percent of topics. Default is \code{0.6}.
#' @param Binarization Logical indicating if to choose one topic for each
#' CellTopic. Default is \code{FALSE}.
#'
#' @return A list with three data frame and one vector. "CellTopic" is a data
#' frame which can be add to a Seurat object. The "domain_topic" is a data frame,
#' row is CellTopic. and col is domain. The "celltype_topic" is a data frame, row is
#' celltype and col is CellTopic. The "Cell_topic" is a vector of which topic be chosen
#' in each CellTopic.
#' @export
#'
FindCellTopic <- function(
    dt_topic_data,
    ct_topic_data,
    spot_clusters,
    cluster = NULL,
    percent = NULL,
    Binarization = FALSE
    ) {
  if (is.null(cluster)) {
    cluster <- colnames(spot_clusters)[1]
  } else if (!cluster %in% colnames(spot_clusters)) {
    cat("You can use these clusters:\n", colnames(spot_clusters), "\n")
    stop("The selected cluster is not in the data!")
  } else if (!identical(as.character(sort(unique(spot_clusters[[cluster]]))),
                        sort(gsub("spot_domain_", "", colnames(dt_topic_data))))) {
    stop("Be sure to use the same cluster in function ksTest")
  }
  if (is.null(percent)) {
    percent <- 0.6
  }
  # Choose topic
  if (Binarization == FALSE) {
    tpc_choose <- lapply(dt_topic_data, TopnTopics, percent)
  } else {
    tpc_choose <- TopTopic(dt_topic_data)
  }
  # Create metadata, cell topic data and domain topic data
  MetaT <- spot_clusters[cluster]
  clst <- as.character(unique(sort(MetaT[[cluster]])))
  MetaT[[cluster]] <- as.character(MetaT[[cluster]])
  for (i in 1:length(clst)) {
    MetaT[which(MetaT[[cluster]] == clst[i]), ] <- as.character(tpc_choose)[i]
  }
  ct_data <- data.frame(matrix(ncol = 0, nrow = nrow(ct_topic_data)))
  rownames(ct_data) <- rownames(ct_topic_data)
  dt_data <- data.frame(matrix(ncol = ncol(dt_topic_data), nrow = 0))
  colnames(dt_data) <- colnames(dt_topic_data)
  Cell_topic <- c()
  for (i in 1:length(unique(tpc_choose))) {
    ct_add <- rep(0, nrow(ct_topic_data))
    dt_add <- rep(0, ncol(dt_topic_data))
    for (topic_add in unique(tpc_choose)[[i]]) {
      ct_add <- ct_add + ct_topic_data[, topic_add] * dt_topic_data[topic_add, i]
      dt_add <- dt_add + dt_topic_data[topic_add, ]
    }
    CellTopic <- paste("CellTopic", i, sep = "")
    ct_data <- cbind(ct_data, ct_add)
    colnames(ct_data)[ncol(ct_data)] <- CellTopic
    dt_data <- rbind(dt_data, dt_add)
    rownames(dt_data)[nrow(dt_data)] <- CellTopic
    MetaT[which(MetaT == unique(as.character(lapply(tpc_choose, unname)))[[i]]), ] <- CellTopic
    Cell_topic <- c(Cell_topic, paste(unique(tpc_choose)[[i]], collapse = "_"))
    names(Cell_topic)[i] <- CellTopic
  }
  colnames(MetaT) <- "CellTopic"
  # Normalization
  ct_data <- as.data.frame(apply(ct_data, 2, function(x) {return(x / sum(x))}))
  # Add each CellTopic to metadata
  for (i in 1:nrow(dt_data)) {
    MetaTs <- spot_clusters[cluster]
    clst <- as.character(unique(sort(MetaTs[[cluster]])))
    MetaTs[[cluster]] <- as.character(MetaTs[[cluster]])
    for (j in 1:length(clst)) {
      MetaTs[which(MetaTs == clst[j]), ] <- dt_data[i, j]
    }
    MetaT <- cbind(MetaT, MetaTs)
    colnames(MetaT)[ncol(MetaT)] <- rownames(dt_data)[i]
  }
  # result
  result_list <- list(
    CellTopic = MetaT,
    domain_topic = dt_data,
    celltype_topic = ct_data,
    Cell_topic = Cell_topic
  )
  return(result_list)
}


#' Select precent of topics for each domain.
#'
#' @param x A vector of integer. A column of date frame domain_topic.
#' @param percent The percent of topics.
#'
#' @return A vector of integer, which indicates the selected topics.
#'
#' @noRd
TopnTopics <- function(x, percent) {
  top <- c(which(x == max(x)))
  sumtpc <- max(x)
  while (sumtpc <= percent) {
    x[which(x == max(x))] <- 0
    top <- c(top, which(x == max(x)))
    sumtpc <- sumtpc + max(x)
  }
  return(top)
}


#' Choose the topic for each domain.
#'
#' @param dt_topic_data A data frame of domain_topic.
#'
#' @return A list with the topic chosen for each domain.
#'
#' @noRd
TopTopic <- function(dt_topic_data){
  tpc_rank <- as.data.frame(t(apply(dt_topic_data, 1, rank)))
  tpc_choose_n <- lapply(tpc_rank, function(x) {which(x == max(x))})
  for (i in 1:length(tpc_choose_n)) {
    if (length(tpc_choose_n[[i]]) > 1) {
      tpc_choose_n[[i]] <- tpc_choose_n[[i]][which(dt_topic_data[tpc_choose_n[[i]], i] == max(dt_topic_data[tpc_choose_n[[i]], i]))]
    }
    tpc_choose_n[[i]] <- unname(tpc_choose_n[[i]])
  }
  return(tpc_choose_n)
}


#' @title MetaTopic
#' @description
#' Make a MetaTopic from CellTopic by celltype scores.
#'
#' @param ct_topic_data A data frame, row is celltype and col is CellTopic.
#' @param k A integer of how much MetaTopic to choose.
#' @param method See \code{\link[stats]{hclust}}.
#'
#' @return A data frame of the cluster result of CellTopic.
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @export
#'
MetaTopic <- function(ct_topic_data, k, method) {
  if(is.null(k)){
    k <- 4
  }
  MetaCellT <- t(scale(ct_topic_data))
  Meta <- stats::hclust(stats::dist(MetaCellT), method = method)
  tree <- as.data.frame(stats::cutree(Meta, k = k))
  colnames(tree) <- "MetaCluster"
  return(tree)
}
