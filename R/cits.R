#' Controlled Interrupted Time Series (CITS) Estimation
#'
#' Fit a GLS-based Controlled Interrupted Time Series (CITS) model with optional ARMA(p,q) autocorrelation.
#' Robust standard errors (CR2) are computed using the \code{clubSandwich} package.
#' Interaction terms are automatically created if not provided.
#'
#' @param data A data frame containing the variables for CITS analysis.
#' @param y_col Outcome variable column name (string).
#' @param T_col Time index column name (string).
#' @param I_col Intervention indicator column name (string).
#'        Numeric: 1 indicates the intervention is applied at that time, 0 otherwise.
#' @param E_col Group indicator column name (string).
#'        Numeric: 1 indicates the treatment/experimental group, 0 indicates the control group.
#' @param TI_col Optional: Column name for the T × I interaction (default = NULL). Will be computed if NULL.
#' @param ET_col Optional: Column name for the E × T interaction (default = NULL). Will be computed if NULL.
#' @param EI_col Optional: Column name for the E × I interaction (default = NULL). Will be computed if NULL.
#' @param ETI_col Optional: Column name for the E × T × I interaction (default = NULL). Will be computed if NULL.
#' @param p_range Range of autoregressive (AR) terms to search (default 0:3).
#' @param q_range Range of moving average (MA) terms to search (default 0:3).
#'
#' @details
#' This function fits a CITS model using generalized least squares (GLS). It automatically calculates
#' interaction terms if they are not provided in the input data. If ARMA fitting fails or produces
#' non-stationary parameters, the function falls back to standard GLS without correlation.
#'
#' The treatment group (`E_col = 1`) is the group that receives the intervention, and the control group
#' (`E_col = 0`) does not receive the intervention. The intervention indicator (`I_col`) marks whether
#' the intervention is applied at a given time point (`1` = applied, `0` = not applied).
#'
#' @return A list containing:
#'   \describe{
#'     \item{model}{The fitted GLS model object.}
#'     \item{robust_se}{CR2 robust covariance matrix from \code{clubSandwich}.}
#'     \item{data}{Data frame with fitted values and standard errors.}
#'     \item{best_p}{Selected AR order based on AIC.}
#'     \item{best_q}{Selected MA order based on AIC.}
#'     \item{arma_used}{Logical indicating whether ARMA correlation was used (TRUE) or GLS fallback was used (FALSE).}
#'   }
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   T = 1:100,
#'   E = rep(c(0,1), each=100),
#'   I = c(rep(0,50), rep(1,50), rep(0,50), rep(1,50)),
#'   y = rnorm(200)
#' )
#' res <- cits(df, y_col="y", T_col="T", I_col="I", E_col="E")
#' summary(res$model)
#' }
#' @export
cits <- function(data,
                 y_col,
                 T_col,
                 I_col,
                 E_col,
                 TI_col = NULL,
                 ET_col = NULL,
                 EI_col = NULL,
                 ETI_col = NULL,
                 p_range = 0:3,
                 q_range = 0:3) {

  df <- data

  if (is.null(TI_col))  df$TI <- df[[T_col]] * df[[I_col]]
  if (is.null(ET_col))  df$ET <- df[[E_col]] * df[[T_col]]
  if (is.null(EI_col))  df$EI <- df[[E_col]] * df[[I_col]]
  if (is.null(ETI_col)) df$ETI <- df[[E_col]] * df[[T_col]] * df[[I_col]]

  df$y <- df[[y_col]]
  df$T <- df[[T_col]]
  df$I <- df[[I_col]]
  df$E <- df[[E_col]]

  mod_formula <- y ~ T + I + TI + E + ET + EI + ETI

  calc_aic <- function(p, q) {
    tryCatch({
      m <- nlme::gls(
        mod_formula,
        data = df,
        correlation = nlme::corARMA(p = p, q = q),
        method = "ML"
      )
      stats::AIC(m)
    }, error = function(e) NA)
  }

  arma_grid <- expand.grid(p = p_range, q = q_range)
  arma_grid <- arma_grid[!(arma_grid$p == 0 & arma_grid$q == 0), ]
  arma_grid$AIC <- mapply(calc_aic, arma_grid$p, arma_grid$q)
  arma_valid <- arma_grid[!is.na(arma_grid$AIC), ]

  if (nrow(arma_valid) == 0) {
    message("All ARMA models failed; fitting GLS without correlation structure")
    model <- nlme::gls(mod_formula, data = df, method = "ML")
    robust_se <- clubSandwich::vcovCR(model, cluster = df$E, type = "CR2")
    preds <- AICcmodavg::predictSE.gls(model, df, se.fit = TRUE)
    df <- dplyr::mutate(df, fitted = as.numeric(preds$fit), se = as.numeric(preds$se))
    return(list(model = model,
                robust_se = robust_se,
                data = df,
                best_p = NA,
                best_q = NA,
                arma_used = FALSE))
  }

  best <- arma_valid[which.min(arma_valid$AIC), ]
  best_p <- best$p
  best_q <- best$q

  model <- nlme::gls(
    mod_formula,
    data = df,
    method = "ML",
    correlation = nlme::corARMA(p = best_p, q = best_q, form = ~ T | E)
  )

  robust_se <- clubSandwich::vcovCR(model, cluster = df$E, type = "CR2")
  preds <- AICcmodavg::predictSE.gls(model, df, se.fit = TRUE)
  df <- dplyr::mutate(df, fitted = as.numeric(preds$fit), se = as.numeric(preds$se))

  return(list(model = model,
              robust_se = robust_se,
              data = df,
              best_p = best_p,
              best_q = best_q,
              arma_used = TRUE))
}
