#' Binary Logistic Regression with a Continuous Predictor
#'
#' Fits a binary logistic regression model with a single continuous predictor,
#' standardised to unit standard deviation (per 1 SD). Reports the log-odds
#' ratio (LOR), odds ratio (OR), Wald confidence interval, p-value,
#' Nagelkerke R-squared, and the full set of closeness-of-concordance
#' measures (CCMs) on the \[-1, +1\] scale.
#'
#' @param outcome Integer or numeric vector. Binary outcome (0/1).
#' @param predictor Numeric vector. Continuous predictor variable. Will be
#'   standardised internally (mean = 0, SD = 1) before fitting.
#' @param alpha Numeric. Significance level for confidence intervals
#'   (default 0.05).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{model}{The fitted `glm` object.}
#'     \item{summary_table}{A `data.frame` with LOR, OR, 95\% CI, p-value,
#'       Nagelkerke R-squared, and CCMs (YuleQ, YuleY, r_meta) with CIs.}
#'     \item{predictor_mean}{Mean of the predictor used for standardisation.}
#'     \item{predictor_sd}{SD of the predictor used for standardisation.}
#'   }
#'
#' @details
#' The predictor is standardised as `(x - mean(x)) / sd(x)`, so the
#' reported OR and LOR correspond to a one-standard-deviation increase.
#' Nagelkerke R-squared is computed as:
#' (1 - exp((2/n)(LL_null - LL_fit))) / (1 - exp((2/n) * LL_null)).
#'
#' @references
#' Kim, S.-K., & Grochowalski, J. H. (2019). Gaining from discretization of
#' continuous data: The correspondence analysis biplot approach.
#' *Behavior Research Methods*, 51(2), 589-601.
#' \doi{10.3758/s13428-018-1161-1}
#'
#' @examples
#' data(lorbridge_data)
#' res <- blr_continuous(outcome  = lorbridge_data$minority,
#'                       predictor = lorbridge_data$VM)
#' print(res$summary_table)
#'
#' @importFrom stats glm binomial coef qnorm logLik
#' @export
blr_continuous <- function(outcome, predictor, alpha = 0.05) {
  stopifnot(length(outcome) == length(predictor))
  z975     <- qnorm(1 - alpha / 2)
  pred_mean <- mean(predictor, na.rm = TRUE)
  pred_sd   <- sd(predictor,   na.rm = TRUE)
  pred_z    <- (predictor - pred_mean) / pred_sd

  df_fit <- data.frame(outcome = outcome, pred_z = pred_z)
  fit    <- glm(outcome ~ pred_z, data = df_fit, family = binomial(link = "logit"))
  sm     <- summary(fit)

  lor    <- coef(fit)["pred_z"]
  se     <- coef(sm)["pred_z", "Std. Error"]
  lor_lo <- lor - z975 * se
  lor_hi <- lor + z975 * se
  pval   <- coef(sm)["pred_z", "Pr(>|z|)"]

  OR <- exp(lor); OR_lo <- exp(lor_lo); OR_hi <- exp(lor_hi)

  # Nagelkerke R-squared
  n       <- nrow(df_fit)
  null_ll <- as.numeric(logLik(glm(outcome ~ 1, data = df_fit,
                                   family = binomial)))
  fit_ll  <- as.numeric(logLik(fit))
  nagel   <- (1 - exp((2/n) * (null_ll - fit_ll))) /
             (1 - exp((2/n) * null_ll))

  ccm <- ccm_row(OR, OR_lo, OR_hi, lor, lor_lo, lor_hi)
  tbl <- data.frame(
    Predictor      = "Continuous (per 1 SD)",
    LOR            = lor,    LOR_lo = lor_lo, LOR_hi = lor_hi,
    OR             = OR,     OR_lo  = OR_lo,  OR_hi  = OR_hi,
    p              = pval,
    Nagelkerke_R2  = nagel,
    YuleQ          = ccm$YuleQ, Q_lo   = ccm$Q_lo,   Q_hi   = ccm$Q_hi,
    YuleY          = ccm$YuleY, Y_lo   = ccm$Y_lo,   Y_hi   = ccm$Y_hi,
    r_meta         = ccm$r_meta, r_lo  = ccm$r_lo,   r_hi   = ccm$r_hi,
    row.names      = NULL
  )
  list(model = fit, summary_table = tbl,
       predictor_mean = pred_mean, predictor_sd = pred_sd)
}


#' Binary Logistic Regression with a Categorical (Ordinal) Predictor
#'
#' Fits a binary logistic regression model with a categorical predictor
#' (treated as an unordered factor with a user-specified reference level).
#' For each non-reference category, reports the pairwise LOR, OR, Wald
#' confidence interval, p-value, CCMs, and cosine theta from a 2 x J simple
#' correspondence analysis (Kim & Grochowalski, 2019 bridge).
#'
#' @param outcome Integer or numeric vector. Binary outcome (0/1).
#' @param predictor Factor or character vector. Categorical predictor. Will be
#'   converted to a factor internally.
#' @param ref_level Character. Reference category label (default: first level).
#' @param alpha Numeric. Significance level for confidence intervals
#'   (default 0.05).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{model}{The fitted glm object.}
#'     \item{results}{A data.frame with one row per non-reference category,
#'       containing LOR, OR, 95 percent CI, p-value, CCMs, and cosine theta.}
#'     \item{tab_2xJ}{The 2 x J contingency table used for cosine theta.}
#'   }
#'
#' @details
#' Cosine theta is computed via a 1D singular value decomposition (SVD) of
#' the standardised residual matrix of the 2 x J correspondence table,
#' bypassing `CAvariants` (which requires at least 2 dimensions). In a
#' 2-row table the 1D CA solution is exact, and cosine theta equals +1 or -1
#' for each non-reference column, with the sign indicating the direction of
#' association relative to the reference category and majority group.
#'
#' @references
#' Kim, S.-K., & Grochowalski, J. H. (2019). Gaining from discretization of
#' continuous data: The correspondence analysis biplot approach.
#' *Behavior Research Methods*, 51(2), 589-601.
#' \doi{10.3758/s13428-018-1161-1}
#'
#' @examples
#' data(lorbridge_data)
#' res <- blr_categorical(outcome   = lorbridge_data$minority,
#'                        predictor = lorbridge_data$VMbin,
#'                        ref_level = "VM4")
#' print(res$results)
#'
#' @importFrom stats glm binomial coef relevel qnorm setNames
#' @export
blr_categorical <- function(outcome, predictor, ref_level = NULL, alpha = 0.05) {
  stopifnot(length(outcome) == length(predictor))
  z975 <- qnorm(1 - alpha / 2)

  pred_fac <- factor(predictor)
  if (!is.null(ref_level)) pred_fac <- relevel(pred_fac, ref = ref_level)
  ref <- levels(pred_fac)[1]

  df_fit <- data.frame(outcome = outcome, predictor = pred_fac)
  fit    <- glm(outcome ~ predictor, data = df_fit,
                family = binomial(link = "logit"))
  sm     <- summary(fit)

  # 2 x J contingency table (row 0 = majority, row 1 = minority)
  tab_2xJ <- table(outcome, pred_fac)
  rownames(tab_2xJ) <- ifelse(rownames(tab_2xJ) == "0", "Majority", "Minority")

  non_ref <- setdiff(levels(pred_fac), ref)

  # 1D CA via SVD for cosine theta
  P      <- tab_2xJ / sum(tab_2xJ)
  r_mass <- rowSums(P); c_mass <- colSums(P)
  S      <- diag(1 / sqrt(r_mass)) %*%
            (P - outer(r_mass, c_mass)) %*%
            diag(1 / sqrt(c_mass))
  sv     <- svd(S, nu = 1, nv = 1)
  F_ca   <- diag(1 / sqrt(r_mass)) %*% sv$u * sv$d[1]
  G_ca   <- diag(1 / sqrt(c_mass)) %*% sv$v * sv$d[1]
  rownames(F_ca) <- rownames(tab_2xJ)
  rownames(G_ca) <- colnames(tab_2xJ)
  g_vec  <- stats::setNames(as.numeric(G_ca[, 1]), rownames(G_ca))
  f_vec  <- stats::setNames(as.numeric(F_ca[, 1]), rownames(F_ca))
  g0     <- g_vec[ref]
  dfi    <- f_vec["Minority"] - f_vec["Majority"]

  results <- do.call(rbind, lapply(non_ref, function(bn) {
    a <- tab_2xJ["Minority", bn]; b <- tab_2xJ["Minority", ref]
    c <- tab_2xJ["Majority", bn]; d <- tab_2xJ["Majority", ref]
    lc   <- lor_ci_2x2(a, b, c, d, alpha)
    OR   <- exp(lc$lor); OR_lo <- exp(lc$lo); OR_hi <- exp(lc$hi)
    pv   <- coef(sm)[paste0("predictor", bn), "Pr(>|z|)"]
    ccm  <- ccm_row(OR, OR_lo, OR_hi, lc$lor, lc$lo, lc$hi)

    dgj <- g_vec[bn] - g0
    cos_t <- if (abs(dfi) < .Machine$double.eps || abs(dgj) < .Machine$double.eps)
      NA_real_ else (dfi * dgj) / (abs(dfi) * abs(dgj))

    data.frame(Category = bn,
               LOR = lc$lor, LOR_lo = lc$lo, LOR_hi = lc$hi,
               OR = OR, OR_lo = OR_lo, OR_hi = OR_hi,
               p = pv,
               YuleQ = ccm$YuleQ, Q_lo = ccm$Q_lo, Q_hi = ccm$Q_hi,
               YuleY = ccm$YuleY, Y_lo = ccm$Y_lo, Y_hi = ccm$Y_hi,
               r_meta = ccm$r_meta, r_lo = ccm$r_lo, r_hi = ccm$r_hi,
               cos_theta = cos_t,
               row.names = NULL)
  }))
  list(model = fit, results = results, tab_2xJ = tab_2xJ)
}


#' Cosine Theta from a 2-Row Correspondence Analysis (SVD)
#'
#' Computes anchored cosine theta values from a 2 x J contingency table via a
#' direct 1D SVD of the standardised residual matrix, bypassing the
#' `CAvariants` function (which requires at least 2 dimensions). This
#' implements the Kim and Grochowalski (2019) log-odds ratio to cosine theta
#' bridge for 2-row tables.
#'
#' @param tab_2xJ A 2 x J matrix or table. Row 1 = majority/reference group;
#'   row 2 = minority/focal group.
#' @param ref_col Character. Name of the reference (anchor) column.
#'
#' @return A named numeric vector of cosine theta values, one per non-reference
#'   column. Values are +1 or -1 in the 1D case (sign carries the direction).
#'
#' @references
#' Kim, S.-K., & Grochowalski, J. H. (2019). Gaining from discretization of
#' continuous data: The correspondence analysis biplot approach.
#' *Behavior Research Methods*, 51(2), 589-601.
#' \doi{10.3758/s13428-018-1161-1}
#'
#' @examples
#' data(lorbridge_data)
#' tab <- table(lorbridge_data$minority, lorbridge_data$VMbin)
#' rownames(tab) <- c("Majority", "Minority")
#' cosine_theta_2row(tab, ref_col = "VM4")
#' @importFrom stats setNames
#' @export
cosine_theta_2row <- function(tab_2xJ, ref_col) {
  stopifnot(nrow(tab_2xJ) == 2, ref_col %in% colnames(tab_2xJ))
  P      <- tab_2xJ / sum(tab_2xJ)
  r_mass <- rowSums(P); c_mass <- colSums(P)
  S      <- diag(1 / sqrt(r_mass)) %*%
            (P - outer(r_mass, c_mass)) %*%
            diag(1 / sqrt(c_mass))
  sv     <- svd(S, nu = 1, nv = 1)
  F_ca   <- diag(1 / sqrt(r_mass)) %*% sv$u * sv$d[1]
  G_ca   <- diag(1 / sqrt(c_mass)) %*% sv$v * sv$d[1]
  rownames(F_ca) <- rownames(tab_2xJ)
  rownames(G_ca) <- colnames(tab_2xJ)
  g_vec  <- stats::setNames(as.numeric(G_ca[, 1]), rownames(G_ca))
  f_vec  <- stats::setNames(as.numeric(F_ca[, 1]), rownames(F_ca))
  g0  <- g_vec[ref_col]
  dfi <- f_vec[2] - f_vec[1]
  non_ref <- setdiff(colnames(tab_2xJ), ref_col)
  sapply(non_ref, function(bn) {
    dgj <- g_vec[bn] - g0
    if (abs(dfi) < .Machine$double.eps || abs(dgj) < .Machine$double.eps)
      return(NA_real_)
    (dfi * dgj) / (abs(dfi) * abs(dgj))
  })
}
