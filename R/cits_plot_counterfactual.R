#' CITS Automatic Analysis + Counterfactual Visualization
#'
#' Fits a `cits()` model, calculates fitted values with 95% confidence intervals,
#' and visualizes the counterfactual for the treatment group (E = 1) after the intervention
#' as a dashed line. A vertical line marks the intervention time.
#'
#' @param res Result of CITS analysis: a list returned by `cits()`.
#' @param y_col Name of the outcome variable (string). This corresponds to `y` in the `cits()` model.
#' @param T_col Name of the time index variable (string). This corresponds to `T` in the `cits()` model.
#' @param I_col Name of the intervention indicator variable (string). This corresponds to `I` in the `cits()` model.
#' @param E_col Name of the group indicator variable (string). This corresponds to `E` in the `cits()` model.
#' @param intervention_time Numeric. Time point at which the intervention occurs. Required for plotting counterfactual.
#' @return A ggplot object with observed points, fitted line, counterfactual line, confidence ribbon, and intervention vertical line.
#' @importFrom stats qnorm predict formula fitted
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon geom_vline annotate labs scale_color_manual scale_fill_manual theme_minimal
#' @export
plot_cits_result_cf <- function(res, y_col = "y", T_col = "T",
                                I_col = "I", E_col = "E", intervention_time) {

  df <- as.data.frame(res$data)
  df$fitted <- as.numeric(df$fitted)
  df$se <- as.numeric(df$se)

  # 95% confidence intervals
  z <- qnorm(0.975)
  df$lwr <- df$fitted - z * df$se
  df$upr <- df$fitted + z * df$se

  # Ensure group variable is factor
  df[[E_col]] <- factor(df[[E_col]])

  # Create counterfactual: after intervention, treatment group's intervention set to 0
  df_cf <- df
  mask <- df_cf[[E_col]] == 1 & df_cf[[T_col]] >= intervention_time
  df_cf[[I_col]][mask] <- 0
  df_cf[[I_col]][!mask & df_cf[[E_col]] == 1] <- df_cf[[I_col]][!mask & df_cf[[E_col]] == 1]

  # Recompute interaction terms
  df_cf$TI  <- df_cf[[T_col]] * df_cf[[I_col]]
  df_cf$ET  <- as.numeric(df_cf[[E_col]]) * df_cf[[T_col]]
  df_cf$EI  <- as.numeric(df_cf[[E_col]]) * df_cf[[I_col]]
  df_cf$ETI <- as.numeric(df_cf[[E_col]]) * df_cf[[T_col]] * df_cf[[I_col]]

  # Prepare newdata for prediction
  model_vars <- all.vars(formula(res$model))[-1] # exclude outcome y
  df_cf_pred <- df_cf[, model_vars, drop = FALSE]
  for(v in model_vars) df_cf_pred[[v]] <- as.numeric(df_cf_pred[[v]])

  # Predict counterfactual
  df_cf$fitted_cf <- predict(res$model, newdata = df_cf_pred)

  # Plot
  plt <- ggplot(df, aes(x = .data[[T_col]], y = .data[[y_col]], color = .data[[E_col]])) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted), linewidth = 0.8) +
    geom_line(data = df_cf[mask, ], aes(y = fitted_cf),
              color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = .data[[E_col]]), alpha = 0.2, color = NA) +
    geom_vline(xintercept = intervention_time, linetype = "dashed", color = "black") +
    annotate("text", x = intervention_time, y = max(df[[y_col]]) * 1.02,
             label = "Intervention", angle = 0, vjust = 0, hjust = 0.5, fontface = "bold") +
    labs(x = "Time", y = "Outcome", color = "Group", fill = "Group") +
    scale_color_manual(values = c("0" = "darkred", "1" = "blue")) +
    scale_fill_manual(values = c("0" = "pink", "1" = "skyblue")) +
    theme_minimal()

  return(plt)
}
