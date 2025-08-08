#' Plot method for quantile.bspcov objects
#'
#' Create visualization plots for posterior quantiles of covariance matrices.
#' Produces heatmap visualizations showing uncertainty across different quantile levels.
#'
#' @param x an object of class \code{quantile.bspcov} from \code{quantile()} function.
#' @param type character string specifying plot type. Options are "heatmap" (default), "uncertainty", or "comparison".
#' @param titles character vector of titles for each quantile plot. If NULL, auto-generated titles are used.
#' @param ncol number of columns in the plot layout. Default is 3.
#' @param color_limits optional vector of length 2 specifying color scale limits. If NULL, computed from data.
#' @param color_low color for low values in heatmap. Default is "black".
#' @param color_high color for high values in heatmap. Default is "white".
#' @param base_size base font size for plot theme. Default is 14.
#' @param legend_position position of legend. Default is "bottom".
#' @param x_label label for x-axis. Default is "Variable".
#' @param y_label label for y-axis. Default is "Variable".
#' @param width plot width when saving. Default is calculated based on number of quantiles.
#' @param height plot height when saving. Default is 6.
#' @param ... additional arguments passed to plotting functions.
#'
#' @return A ggplot object (single quantile) or patchwork object (multiple quantiles) showing heatmap visualizations.
#' @author Kyeongwon Lee
#' @seealso quantile, plot.bspcov, plot.postmean.bspcov
#'
#' @importFrom ggplot2 ggplot aes geom_raster labs theme_bw theme element_text element_blank scale_fill_gradient guide_colourbar
#' @importFrom patchwork plot_layout
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#' 
#' n <- 25
#' p <- 50
#' Sigma0 <- diag(1, p)
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)
#' res <- bspcov::bandPPP(X, 2, 0.01, nsample = 100)
#' quant <- quantile(res)
#' 
#' # Plot uncertainty visualization
#' plot(quant)
#' 
#' # Plot with custom titles and labels
#' plot(quant, titles = c("Lower Bound", "Median", "Upper Bound"),
#'      x_label = "Variable", y_label = "Variable")
#' 
#' # Plot with gene-specific labels
#' plot(quant, x_label = "Gene", y_label = "Gene")
#' 
#' # Plot comparison layout
#' plot(quant, type = "comparison", ncol = 3)
#'
plot.quantile.bspcov <- function(x, 
                                type = "heatmap",
                                titles = NULL,
                                ncol = 3,
                                color_limits = NULL,
                                color_low = "black", 
                                color_high = "white",
                                base_size = 14,
                                legend_position = "bottom",
                                x_label = "Variable",
                                y_label = "Variable",
                                width = NULL,
                                height = 6,
                                ...) {
  
  # Validate inputs
  if (!inherits(x, "quantile.bspcov")) {
    stop("Input must be of class 'quantile.bspcov'")
  }
  
  if (!type %in% c("heatmap", "uncertainty", "comparison")) {
    stop("type must be one of 'heatmap', 'uncertainty', or 'comparison'")
  }
  
  # Get quantile names and probabilities
  quant_names <- names(x)
  n_quantiles <- length(x)
  
  # Auto-generate titles if not provided
  if (is.null(titles)) {
    probs <- gsub("q", "", quant_names)
    titles <- paste0("Quantile ", as.numeric(probs) * 100, "%")
  }
  
  # Function to create individual heatmap
  vis_quantile <- function(quant_matrix, title) {
    p <- ncol(quant_matrix)
    quant_mat <- matrix(as.numeric(quant_matrix), nrow = p, ncol = p)
    df <- reshape2::melt(quant_mat, varnames = c("x", "y"), value.name = "cov")
    
    ggplot2::ggplot(df, ggplot2::aes(x, y, fill = cov)) +
      ggplot2::geom_raster(interpolate = FALSE) +
      ggplot2::labs(title = title, x = x_label, y = y_label) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = legend_position,
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      )
  }
  
  # Compute color scale limits if not provided
  if (is.null(color_limits)) {
    all_values <- unlist(lapply(x, as.numeric))
    color_limits <- c(min(all_values, na.rm = TRUE), max(all_values, na.rm = TRUE))
  }
  
  # Create color scale
  cont_scale <- ggplot2::scale_fill_gradient(
    name = "",
    low = color_low,
    high = color_high,
    limits = color_limits,
    guide = ggplot2::guide_colourbar(
      barwidth = ggplot2::unit(8, "cm"), 
      barheight = ggplot2::unit(0.5, "cm")
    )
  )
  
  # Create plots for each quantile
  plots <- list()
  for (i in 1:n_quantiles) {
    plots[[i]] <- vis_quantile(x[[i]], titles[i])
  }
  
  # Handle single quantile case
  if (n_quantiles == 1) {
    return(plots[[1]] + cont_scale)
  }
  
  # Combine plots using patchwork
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- Reduce(`+`, plots) + 
      patchwork::plot_layout(ncol = ncol, guides = "collect") &
      cont_scale &
      ggplot2::theme(legend.position = legend_position)
    
    return(combined_plot)
  } else {
    warning("patchwork package not available. Returning list of individual plots.")
    return(plots)
  }
}

#' Save quantile plot to file
#'
#' Convenience function to save quantile plots with appropriate dimensions.
#'
#' @param x an object of class \code{quantile.bspcov}.
#' @param filename filename to save the plot.
#' @param width plot width. If NULL, calculated based on number of quantiles.
#' @param height plot height. Default is 6.
#' @param ... additional arguments passed to \code{plot.quantile.bspcov} and \code{ggsave}.
#'
#' @export
#' @examples
#' \dontrun{
#' n <- 25
#' p <- 50
#' Sigma0 <- diag(1, p)
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)
#' res <- bspcov::bandPPP(X, 2, 0.01, nsample = 100)
#' quant <- quantile(res)
#' 
#' # Save uncertainty plot
#' save_quantile_plot(quant, "uncertainty_plot.png")
#' }
save_quantile_plot <- function(x, filename, width = NULL, height = 6, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for saving plots")
  }
  
  # Calculate width if not provided
  if (is.null(width)) {
    n_quantiles <- length(x)
    width <- max(8, n_quantiles * 4)
  }
  
  # Create plot
  p <- plot.quantile.bspcov(x, ...)
  
  # Save plot
  ggplot2::ggsave(filename, p, width = width, height = height, ...)
  
  message("Plot saved to: ", filename)
}
