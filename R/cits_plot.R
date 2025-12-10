#' CITS Automatic Analysis + Visualization
#'
#' Fits a `cits()` model, calculates fitted values with 95% confidence intervals,
#' and visualizes the results using ggplot2. Optionally, a vertical line can be drawn
#' at the intervention time if `intervention_time` is provided.
#'
#' @param res Result of CITS analysis: a list returned by `cits()`.
#' @param y_col Name of the outcome variable (string). This corresponds to `y` in the `cits()` model.
#' @param T_col Name of the time index variable (string). This corresponds to `T` in the `cits()` model.
#' @param E_col Name of the group indicator variable (string). This corresponds to `E` in the `cits()` model.
#' @param intervention_time Optional numeric. Draws a vertical line at this time if provided. Corresponds to intervention point marked by `I` in the data.
#' @return A ggplot object with observed points, fitted lines, 95% confidence ribbon, and optional intervention line.
#' @importFrom stats qnorm
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon geom_vline annotate labs scale_color_manual scale_fill_manual theme_minimal
#' @export
plot_cits_result <- function(res, y_col = "y", T_col = "T", E_col = "E", intervention_time = NULL) {

  df <- as.data.frame(res$data)
  df$fitted <- as.numeric(df$fitted)
  df$se <- as.numeric(df$se)

  # Compute 95% confidence intervals
  alpha <- 0.05
  z <- qnorm(1 - alpha/2)
  df$lwr <- df$fitted - z * df$se
  df$upr <- df$fitted + z * df$se

  # Ensure group variable is a factor
  df[[E_col]] <- factor(df[[E_col]])

  # Base plot: points, fitted line, confidence ribbon
  plt <- ggplot(df, aes(x = .data[[T_col]], y = .data[[y_col]], color = .data[[E_col]])) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted), linewidth = 0.8) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = .data[[E_col]]), alpha = 0.2, color = NA) +
    labs(x = "Time", y = "Outcome", color = "Group", fill = "Group") +
    scale_color_manual(values = c("0" = "darkred", "1" = "blue")) +
    scale_fill_manual(values = c("0" = "pink", "1" = "skyblue")) +
    theme_minimal()

  # Optional vertical intervention line
  if (!is.null(intervention_time)) {
    plt <- plt +
      geom_vline(xintercept = intervention_time, linetype = "dashed", color = "black") +
      annotate("text", x = intervention_time, y = max(df[[y_col]]) * 1.02,
               label = "Intervention", angle = 0, vjust = 0, hjust = 0.5, fontface = "bold")
  }

  return(plt)
}
