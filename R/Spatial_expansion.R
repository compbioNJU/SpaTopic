#' @title Spatial_expansion
#' @description
#' Buffers and angles for dividing the specified spatial domain.
#'
#' @param st_obj A Seurat object.
#' @param CellTopic A character of which CellTopic to be chosen.
#' @param meta.col A character of which meta.data to be chosen. Default is \code{"CellTopic"}.
#' @param distance A numeric. The distance of expansion. Default is \code{NULL},
#' use the nearest spot distance.
#' @param alpha A numeric used for alpha-shape calculation. Default is \code{NULL},
#' use nearest spot distance based on data type. See also \code{\link[alphahull]{ashape}}.
#' @param type By default \code{type = "Hexagon"} for data with hexagonal spot arrangement.
#' When the spot arrangement of data is square, set \code{tpye = "Square"}.
#' @param edge.plot Logical indicating if to display the regional division plot
#' of CellTopic. Default is \code{TRUE}.
#' @param image.factor A character of the resolution scale factor. Default is \code{"lowres"}.
#'
#' @return A Seurat object. The expansion of the data is added to the \code{Seurat@meta.data}
#' and the message is added to \code{Seurat@misc}.
#' @import Seurat
#' @import sf
#' @import ggplot2
#' @importFrom alphahull ashape
#' @export
#'

Spatial_expansion <- function(st_obj, CellTopic,
                              meta.col = "CellTopic",
                              distance = NULL,
                              alpha = NULL,
                              type = c("Hexagon", "Square"),
                              edge.plot = TRUE,
                              image.factor = "lowres") {
  coordinate <- SeuratObject::GetTissueCoordinates(st_obj, scale = NULL, cols = c("imagerow", "imagecol"))
  cat("Use the low resolution image\n")
  scale.factor <- Seurat::ScaleFactors(st_obj[["image"]])[[image.factor]]
  df <- base::cbind(coordinate, st_obj@meta.data[meta.col])
  df$expansion <- NA
  df_select <- df[df[[meta.col]] %in% CellTopic, ]
  type <- base::match.arg(type)
  cat("The type of data is", type, "\n")
  if (is.null(alpha)) {
    switch(type,
           "Hexagon" = alpha <- min(stats::dist(df_select[, c("imagecol", "imagerow")])) * 1.16,
           "Square" = alpha <- min(stats::dist(df_select[, c("imagecol", "imagerow")])) * 1.41
    )
  }
  cat("The alpha value of alpha-shape method is", alpha, "\n")
  sf_data <- sf::st_as_sf(df, coords = c("imagecol", "imagerow"))
  sf_data_select <- sf::st_as_sf(df_select, coords = c("imagecol", "imagerow"))
  ashape_obj <- alphahull::ashape(df_select[, c("imagecol", "imagerow")], alpha = alpha)
  df_ashape <- ashape_obj$edges[, c("x1", "y1", "x2", "y2")]
  line_list <- base::apply(df_ashape, 1, simplify = F, function(x) {
    sf::st_linestring(base::matrix(x, ncol = 2, byrow = T))
  })
  line <- sf::st_line_merge(sf::st_sfc(sf::st_multilinestring(line_list)))
  polygon <- sf::st_cast(sf::st_polygonize(line))
  single_polygon <- base::apply(sf::st_within(sf_data_select, polygon, sparse = F), 2, any)
  polygon <- polygon[single_polygon]
  line_fix <- sf::st_cast(polygon, "LINESTRING")
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = polygon) +
    ggplot2::geom_point(data = df_select, ggplot2::aes(x = imagecol, y = imagerow), color = "red") +
    ggplot2::theme_void()
  if (edge.plot) {
    print(p)
  }
  single_spot_T <- base::apply(sf::st_intersects(polygon, sf_data_select, sparse = F), 2, any)
  single_spot <- row.names(sf_data_select)[!single_spot_T]
  sf_data[single_spot, meta.col] <- "SingleSpot"
  in_spot <- base::apply(sf::st_intersects(polygon, sf_data, sparse = F), 2, any)
  sf_data[in_spot, meta.col] <- "InSpot"
  if (is.null(distance)) {
    distance <- base::min(stats::dist(df[, c("imagecol", "imagerow")]))
  }
  cat("The distance of expansion is", distance, "\n")
  i <- 1
  while (any(is.na(sf_data$expansion))) {
    area <- sf::st_buffer(line_fix, dist = distance * i)
    sf_data$covers <- sf::st_covers(area, sf_data, sparse = FALSE)[1, ]
    sf_data$expansion[sf_data$covers & sf_data[[meta.col]] != "InSpot" & is.na(sf_data$expansion)] <- distance * scale.factor * i
    sf_data$expansion[sf_data$covers & sf_data[[meta.col]] == "InSpot" & is.na(sf_data$expansion)] <- -distance * scale.factor * i
    i <- i + 1
  }
  sf_data[[base::paste0(c("expansion", CellTopic), collapse = "_")]] <- sf_data$expansion
  st_obj <- SeuratObject::AddMetaData(st_obj, as.data.frame(sf_data)[base::paste0(c("expansion", CellTopic), collapse = "_")])
  center <- sf::st_coordinates(sf::st_centroid(line_fix))
  angle <- apply(sf_data, 1, function(x, center) {
    point <- sf::st_coordinates(x$geometry)
    base::atan2(point[2] - center[2], point[1] - center[1])
  }, center)
  angle <- angle * 180 / pi + 180
  st_obj <- SeuratObject::AddMetaData(st_obj, angle, col.name = base::paste0(c("angle", CellTopic), collapse = "_"))
  st_obj@misc$expansion[[paste0(CellTopic, collapse = "_")]] <- list(
    distance = distance,
    alpha = alpha,
    type = type,
    image = image.factor,
    image.factor = scale.factor,
    edge.plot = p
  )
  cat("The expansion of the data is added to the Seurat@meta.data and the message is added to Seurat@misc\n")
  return(st_obj)
}



#' @title Expansion_gene_plot
#' @description
#' Show smooth plots of genes as a function of buffer distance at different angles.
#'
#' @param st_obj A Seurat object.
#' @param CellTopic A character of which CellTopic to be chosen.
#' @param gene A vector contains one or more genes for display
#' @param meta.col A character of which meta.data to be chosen. Default is \code{"CellTopic"}.
#' @param angle A numeric indicates the starting angle. Default is \code{90}, which means directly above.
#' @param parts A number indicates that the angle needs to be divided into several equal parts.
#' @param plot Logical indicating if to display the result plot or just data for
#' plotting. Default is \code{TRUE}.
#' @param cols A vector of colors used for plot.
#'
#' @return If \code{plot = TRUE}, return smooth plots of genes as a function of
#' buffer distance at different angles. If \code{plot = FALSE}, return the data
#' for plotting.
#' @import ggplot2
#' @import cowplot
#' @importFrom SeuratObject GetAssayData
#' @importFrom tidyr gather
#' @export
#'

Expansion_gene_plot <- function(st_obj, CellTopic, gene,
                                meta.col = "CellTopic",
                                angle = 90,
                                parts = 1,
                                plot = TRUE,
                                cols = NULL) {
  expansion_name <- paste0(c("expansion", CellTopic), collapse = "_")
  angle_name <- paste0(c("angle", CellTopic), collapse = "_")
  if (is.null(st_obj@misc$expansion[[paste0(CellTopic, collapse = "_")]])) {
    stop("The expansion is missing, please run CellTopic_expansion first.")
  }
  if (parts == 1) {
    data <- SeuratObject::GetAssayData(st_obj, layer = "data")
    meta <- st_obj@meta.data[, c(expansion_name, angle_name)]
    colnames(data) <- meta[[expansion_name]]
    data <- data[, order(as.numeric(colnames(data)))]
    gene_expansion <- t(base::as.matrix(data[gene, , drop = F]))
    gene_expansion_df <- data.frame(gene_expansion, expansion = as.numeric(row.names(gene_expansion)))
    gene_expansion_long <- tidyr::gather(gene_expansion_df, key = "gene", value = "expression", -expansion)
    gene_expansion_long$angle <- "360"
    if (plot) {
      p <- ggplot2::ggplot(gene_expansion_long) +
        ggplot2::geom_smooth(aes(x = expansion, y = expression, color = gene), se = FALSE) +
        ggplot2::theme_minimal()
      if (!is.null(cols)) {
        p <- p + ggplot2::scale_color_manual(values = cols)
      }
      return(p)
    }
    return(gene_expansion_long)
  }
  angle_seq <- seq(0, 360, by = 360 / parts)
  angle_seq <- base::Map(function(x) {
    if (x > 360) {
      x <- x - 360
    }
    x
  }, (angle_seq + angle))
  data <- SeuratObject::GetAssayData(st_obj, layer = "data")
  meta <- st_obj@meta.data[, c(expansion_name, angle_name)]
  colnames(data) <- meta[[expansion_name]]
  gene_expansion_long_all <- data.frame()
  if (plot) {
    p_list <- list()
    sincos <- data.frame(row.names = c("sin", "cos"))
    label <- c()
  }
  for (i in 1:(length(angle_seq) - 1)) {
    if (angle_seq[[i]] <= angle_seq[[i + 1]]) {
      data_select <- data[, meta[[angle_name]] >= angle_seq[[i]] & meta[[angle_name]] < angle_seq[[i + 1]]]
    } else {
      data_select <- data[, meta[[angle_name]] >= angle_seq[[i]] | meta[[angle_name]] < angle_seq[[i + 1]]]
    }
    data_select <- data_select[, order(as.numeric(colnames(data_select)))]
    gene_expansion <- t(base::as.matrix(data_select[gene, , drop = F]))
    gene_expansion_df <- data.frame(gene_expansion, expansion = as.numeric(row.names(gene_expansion)))
    gene_expansion_long <- tidyr::gather(gene_expansion_df, key = "gene", value = "expression", -expansion)
    gene_expansion_long$angle <- paste0(angle_seq[i], "-", angle_seq[i + 1])
    label <- c(label, base::paste0(angle_seq[i], "-", angle_seq[i + 1]))
    gene_expansion_long_all <- base::rbind(gene_expansion_long_all, gene_expansion_long)
    if (plot) {
      mid_angle <- ifelse(angle_seq[[i]] <= angle_seq[[i + 1]],
                          (angle_seq[[i + 1]] - angle_seq[[i]]) / 2 + angle_seq[[i]],
                          (angle_seq[[i + 1]] - angle_seq[[i]] + 360) / 2 + angle_seq[[i]]
      )
      sincos[, i] <- c(sinpi((mid_angle - 90) / 180), cospi((mid_angle - 90) / 180))
      p_list[[i]] <- ggplot2::ggplot(gene_expansion_long) +
        ggplot2::geom_smooth(aes(x = expansion, y = expression, color = gene), se = FALSE) +
        ggplot2::theme_minimal()
      if (!is.null(cols)) {
        p_list[[i]] <- p_list[[i]] + ggplot2::scale_color_manual(values = cols)
      }
    }
  }
  if (plot) {
    sincos <- t(sincos)
    sincos[, 1][sincos[, 1] < 0 | (sincos[, 1] == 0 & sincos[, 2] < 0)] <- -1
    sincos[, 1][sincos[, 1] > 0 | (sincos[, 1] == 0 & sincos[, 2] > 0)] <- 1
    p_list <- p_list[base::order(sincos[, 1], -sincos[, 2])]
    label <- label[base::order(sincos[, 1], -sincos[, 2])]
    sincos <- sincos[base::order(sincos[, 1], -sincos[, 2]), ]
    legend <- cowplot::get_legend(p_list[[1]])
    min_x <- max_x <- 0
    min_y <- max_y <- 0
    for (i in 1:length(p_list)) {
      p_list[[i]] <- p_list[[i]] +
        ggplot2::ggtitle(label[i]) +
        ggplot2::theme(plot.title = element_text(size = 10, hjust = 0.5))
      if (min_x > ggplot2::layer_scales(p_list[[i]])$x$range$range[1]) min_x <- ggplot2::layer_scales(p_list[[i]])$x$range$range[1]
      if (max_x < ggplot2::layer_scales(p_list[[i]])$x$range$range[2]) max_x <- ggplot2::layer_scales(p_list[[i]])$x$range$range[2]
      if (min_y > ggplot2::layer_scales(p_list[[i]])$y$range$range[1]) min_y <- ggplot2::layer_scales(p_list[[i]])$y$range$range[1]
      if (max_y < ggplot2::layer_scales(p_list[[i]])$y$range$range[2]) max_y <- ggplot2::layer_scales(p_list[[i]])$y$range$range[2]
    }
    p_list <- base::lapply(p_list, function(x) {
      x + ggplot2::coord_cartesian(ylim = c(min_y, max_y)) +
        ggplot2::theme(legend.position = "none")
    })
    for(i in 1:length(p_list)){
      if(i < parts/2){
        p_list[[i]] <- p_list[[i]] +
          ggplot2::scale_x_reverse(limits = c(max_x, min_x)) +
          ggplot2::theme(axis.title.x = element_blank(), axis.text.x = element_blank())
      }else if(i == ceiling(parts/2)){
        p_list[[i]] <- p_list[[i]] +
          ggplot2::scale_x_reverse(limits = c(max_x, min_x))
      }else if(i == parts){
        p_list[[i]] <- p_list[[i]] +
          ggplot2::scale_x_continuous(limits = c(min_x, max_x)) +
          ggplot2::theme(axis.title.y = element_blank(), axis.text.y = element_blank())
      }else{
        p_list[[i]] <- p_list[[i]] +
          ggplot2::scale_x_continuous(limits = c(min_x, max_x)) +
          ggplot2::theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                         axis.title.x = element_blank(), axis.text.x = element_blank())
      }
    }
    p <- cowplot::plot_grid(plotlist = p_list, ncol = 2, byrow = FALSE, align = 'h')
    p <- cowplot::plot_grid(p, legend, rel_widths = c(4, 1))
    return(p)
  }
  return(gene_expansion_long_all)
}
