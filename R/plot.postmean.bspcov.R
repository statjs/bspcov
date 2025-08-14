#' Plot method for postmean.bspcov objects
#'
#' Create heatmap visualization for posterior mean estimate of sparse covariance matrix.
#' Provides flexible visualization options with customizable aesthetics and labeling.
#'
#' @param x an object of class \code{postmean.bspcov} from \code{estimate()} function.
#' @param title character string for plot title. If NULL, auto-generated title is used.
#' @param color_limits optional vector of length 2 specifying color scale limits. If NULL, computed from data.
#' @param color_low color for low values in heatmap. Default is "black".
#' @param color_high color for high values in heatmap. Default is "white".
#' @param base_size base font size for plot theme. Default is 14.
#' @param legend_position position of legend. Default is "bottom".
#' @param x_label label for x-axis. Default is "Variable".
#' @param y_label label for y-axis. Default is "Variable".
#' @param show_values logical indicating whether to display values on tiles. Default is FALSE.
#' @param ... additional arguments passed to plotting functions.
#'
#' @return A ggplot object showing heatmap visualization of the posterior mean covariance matrix.
#' @author Seongil Jo, Kyeongwon Lee
#' @seealso estimate, plot.bspcov, plot.quantile.bspcov
#'
#' @importFrom ggplot2 ggplot aes geom_raster labs theme_bw theme element_text element_blank scale_fill_gradient guide_colourbar geom_text
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#' \donttest{
#' # Example with simulated data
#' n <- 25
#' p <- 50
#' Sigma0 <- diag(1, p)
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)
#' res <- bspcov::thresPPP(X, eps=0.01, thres=list(value=0.5,fun='hard'), nsample=100)
#' est <- bspcov::estimate(res)
#' 
#' # Basic plot
#' plot(est)
#' 
#' # Plot with custom color scheme
#' plot(est, color_low = "blue", color_high = "red")
#' }
#'
plot.postmean.bspcov <- function(x, 
                                title = NULL,
                                color_limits = NULL,
                                color_low = "black", 
                                color_high = "white",
                                base_size = 14,
                                legend_position = "bottom",
                                x_label = "Variable",
                                y_label = "Variable",
                                show_values = FALSE,
                                ...) {
  
  # Validate inputs
  if (!inherits(x, "postmean.bspcov")) {
    stop("Input must be of class 'postmean.bspcov'")
  }
  
  # Auto-generate title if not provided
  if (is.null(title)) {
    title <- "Posterior Mean"
  }
  
  # Convert matrix to data frame for plotting
  p <- ncol(x)
  cov_mat <- matrix(as.numeric(x), nrow = p, ncol = p)
  df <- reshape2::melt(cov_mat, varnames = c("x", "y"), value.name = "cov")
  
  # Compute color scale limits if not provided
  if (is.null(color_limits)) {
    color_limits <- c(min(df$cov, na.rm = TRUE), max(df$cov, na.rm = TRUE))
  }
  
  # Create base plot
  p_plot <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = cov)) +
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
    ) +
    ggplot2::scale_fill_gradient(
      name = "",
      low = color_low,
      high = color_high,
      limits = color_limits,
      guide = ggplot2::guide_colourbar(
        barwidth = ggplot2::unit(8, "cm"), 
        barheight = ggplot2::unit(0.5, "cm")
      )
    )
  
  # Add values on tiles if requested (useful for small matrices)
  if (show_values) {
    p_plot <- p_plot + 
      ggplot2::geom_text(ggplot2::aes(label = round(cov, 3)), 
                        color = "gray50", size = base_size/4)
  }
  
  return(p_plot)
}
