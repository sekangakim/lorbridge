#' Fit a DONSCA Model
#'
#' Fits a Doubly-Ordered Nonsymmetric Correspondence Analysis (DONSCA) model
#' via `CAvariants`, using all available dimensions.
#'
#' @param tab A numeric matrix or table. Both rows and columns must represent
#'   ordered categories.
#'
#' @return A `CAvariants` fit object containing, among others,
#'   `Rprinccoord` (row principal coordinates) and `Cstdcoord` (column
#'   standard coordinates).
#'
#' @importFrom CAvariants CAvariants
#' @export
donsca_fit <- function(tab) {
  D <- min(nrow(tab) - 1, ncol(tab) - 1)
  CAvariants(tab, catype = "DONSCA", firstaxis = 1, lastaxis = D)
}


#' Doubly-Anchored Cosine Theta for DONSCA
#'
#' Computes doubly-anchored cosine theta values for all non-anchor row and
#' column contrasts in a DONSCA solution. Each cosine theta quantifies the
#' geometric alignment between the direction (row_i - row_anchor) and
#' (col_j - col_anchor) in the full multivariate CA space.
#'
#' @param fit A fitted DONSCA object from `donsca_fit()`.
#' @param col_anchor_idx Integer. Column index of the anchor (reference) column.
#' @param row_anchor_idx Integer. Row index of the anchor (reference) row.
#' @param dims Integer or `"all"`. Number of dimensions to use (default `"all"`).
#'
#' @return A `data.frame` with columns: Row, Col, cos_theta.
#'
#' @export
donsca_cosines <- function(fit, col_anchor_idx, row_anchor_idx, dims = "all") {
  F <- fit$Rprinccoord; G <- fit$Cstdcoord
  if (is.null(F) || is.null(G)) stop("fit must contain Rprinccoord and Cstdcoord.")
  K_shared <- min(ncol(F), ncol(G))
  K <- if (identical(dims, "all")) K_shared else min(K_shared, as.integer(dims))
  F <- F[, seq_len(K), drop = FALSE]; G <- G[, seq_len(K), drop = FALSE]
  fi0 <- as.numeric(F[row_anchor_idx, , drop = FALSE])
  gj0 <- as.numeric(G[col_anchor_idx, , drop = FALSE])
  I <- nrow(F); J <- nrow(G)
  out <- list(); k <- 0L
  for (i in seq_len(I)) {
    if (i == row_anchor_idx) next
    dfi <- as.numeric(F[i, ]) - fi0
    nfi <- sqrt(sum(dfi^2)); if (nfi == 0) next
    for (j in seq_len(J)) {
      if (j == col_anchor_idx) next
      dgj <- as.numeric(G[j, ]) - gj0
      ngj <- sqrt(sum(dgj^2)); if (ngj == 0) next
      k <- k + 1L
      out[[k]] <- data.frame(
        Row = rownames(F)[i], Col = rownames(G)[j],
        cos_theta = sum(dfi * dgj) / (nfi * ngj),
        stringsAsFactors = FALSE)
    }
  }
  if (k == 0L) return(data.frame(Row = character(), Col = character(),
                                  cos_theta = numeric()))
  do.call(rbind, out[seq_len(k)])
}


#' Multinomial Logistic Regression with CCMs (DONSCA Bridge)
#'
#' Fits a multinomial logistic regression model with a single standardised
#' continuous predictor (per 1 SD), using a specified baseline outcome
#' category. Returns log-odds ratios, odds ratios, Wald CIs, and the full
#' set of CCMs for each non-baseline outcome level.
#'
#' @param outcome Factor or character vector. Polytomous outcome variable.
#' @param predictor Numeric vector. Continuous predictor. Standardised
#'   internally to unit SD before fitting.
#' @param baseline Character. Baseline/reference level of the outcome
#'   (default: first level of `outcome` as a factor).
#' @param alpha Numeric. Significance level for CIs (default 0.05).
#'
#' @return A `data.frame` with one row per non-baseline outcome level,
#'   containing LOR, OR, 95% CI, and CCMs (YuleQ, YuleY, r_meta).
#'
#' @importFrom nnet multinom
#' @importFrom stats coef sd
#' @export
mlr_ccm <- function(outcome, predictor, baseline = NULL, alpha = 0.05) {
  z975     <- qnorm(1 - alpha / 2)
  pred_z   <- (predictor - mean(predictor, na.rm = TRUE)) /
               sd(predictor, na.rm = TRUE)
  out_fac  <- factor(outcome)
  if (!is.null(baseline)) out_fac <- relevel(out_fac, ref = baseline)

  df_fit   <- data.frame(outcome = out_fac, pred_z = pred_z)
  fit      <- nnet::multinom(outcome ~ pred_z, data = df_fit, trace = FALSE)
  sm       <- summary(fit)
  coef_mat <- sm$coefficients
  se_mat   <- sm$standard.errors
  nb_levs  <- rownames(coef_mat)

  slope    <- coef_mat[, "pred_z"]
  slope_se <- se_mat[,  "pred_z"]
  LOR      <- slope
  LOR_lo   <- slope - z975 * slope_se
  LOR_hi   <- slope + z975 * slope_se
  OR       <- exp(LOR); OR_lo <- exp(LOR_lo); OR_hi <- exp(LOR_hi)

  do.call(rbind, lapply(nb_levs, function(lv) {
    ccm <- ccm_row(OR[lv], OR_lo[lv], OR_hi[lv],
                   LOR[lv], LOR_lo[lv], LOR_hi[lv])
    data.frame(Level = lv,
               LOR = LOR[lv], LOR_lo = LOR_lo[lv], LOR_hi = LOR_hi[lv],
               OR  = OR[lv],  OR_lo  = OR_lo[lv],  OR_hi  = OR_hi[lv],
               YuleQ = ccm$YuleQ, Q_lo = ccm$Q_lo, Q_hi = ccm$Q_hi,
               YuleY = ccm$YuleY, Y_lo = ccm$Y_lo, Y_hi = ccm$Y_hi,
               r_meta = ccm$r_meta, r_lo = ccm$r_lo, r_hi = ccm$r_hi,
               row.names = NULL)
  }))
}
