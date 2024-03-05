# CellTopic Plot

#' CellTopic Plot
#' @description
#' A drawing function that outputs a CellTopic spatial distribution plot and corresponding celltype barpolt.
#'
#' @param st_obj A seurat object with CellTopic meta data.
#' @param celltype_topic A data frame, row is celltype and col is CellTopic.
#' @param celltopic A character of which CellTopic to plot.
#' @param highlight.circle A character of which meta.data to high light. Default is \code{NULL}.
#' @param cols.circle A character of circle's color. Default is \code{black}.
#' @param cols.highlight A vector with two character of highlight spots' colors.
#' @param cols.celltype A vector with character of celltypes' colors. Default is \code{NULL}.
#' @param plot.size A numeric of plot's size. Default is \code{7}.
#' @param nrow A numeric of rows in the plot grid.
#' @param ncol A numeric of cols in the plot grid.
#' @param ... Arguments passed to other methods.
#'
#' @return A polt with SpatialDimPlot and celltype barpolt.
#' @import Seurat
#' @import ggplot2
#' @importFrom stats reorder
#' @importFrom ggforce geom_circle
#' @importFrom cowplot plot_grid
#' @export
#'

CellTopic_plot <- function(st_obj,
                           celltype_topic,
                           celltopic = colnames(celltype_topic),
                           highlight.circle = NULL,
                           cols.circle = "black",
                           cols.highlight = c("#DE2D26", "grey50"),
                           cols.celltype = NULL,
                           plot.size = 7,
                           nrow = NULL,
                           ncol = NULL,
                           ...){
  Seurat::Idents(st_obj) <- st_obj$CellTopic
  bar_plot_data <- celltype_topic
  bar_plot_data$CellType <- rownames(celltype_topic)
  plot_list <- list()
  for(i in celltopic){
    p <- Seurat::SpatialDimPlot(st_obj, image.alpha = 0,
                                cells.highlight = Seurat::CellsByIdentities(object = st_obj, idents = i),
                                cols.highlight = cols.highlight, facet.highlight = TRUE, ...)
    if(!is.null(highlight.circle)){
      TLS_location <- Seurat::GetTissueCoordinates(st_obj)[st_obj@meta.data[[highlight.circle]],]
      middle <- sum(layer_scales(p)$y$range$range)
      dim_plot <- p + ggforce::geom_circle(aes(x0 = .data$imagecol, y0 = !!middle-.data$imagerow, r = 3),
                                           data = TLS_location, color = cols.circle, inherit.aes = F, ...)
    }
    else{
      dim_plot <- p
    }

    bar_plot <- ggplot2::ggplot(bar_plot_data,
                                aes(x = stats::reorder(.data$CellType, .data[[i]], decreasing = TRUE),
                                    y = .data[[i]])) +
      ggplot2::geom_bar(aes(fill = .data$CellType), stat='identity', width = 0.7) +
      ggplot2::xlab('CellType') +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 8),
            axis.title.y = element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    if(!is.null(cols.celltype)){
      bar_plot <- bar_plot + ggplot2::scale_fill_manual(values = cols.celltype)
    }
    plot_list[[i]] <- cowplot::plot_grid(plotlist = list(dim_plot, bar_plot), ncol = 1)
  }
  plot <- cowplot::plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol)
  return(plot)
}

