#' Log-Odds Ratio with Haldane-Anscombe Correction
#'
#' Computes the log-odds ratio (LOR) and its Wald confidence interval for a
#' 2x2 contingency table, applying the Haldane-Anscombe continuity correction
#' (adding 0.5 to all cells) when any cell count is zero.
#'
#' @param a Numeric. Cell count: row 1 (focal), column 1 (focal).
#' @param b Numeric. Cell count: row 1 (focal), column 2 (reference/anchor).
#' @param c Numeric. Cell count: row 2 (reference/anchor), column 1 (focal).
#' @param d Numeric. Cell count: row 2 (reference/anchor), column 2 (reference/anchor).
#' @param alpha Numeric. Significance level for the confidence interval (default 0.05).
#' @param cc Numeric. Continuity correction added to all cells when any count
#'   is zero (default 0.5, Haldane-Anscombe).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{lor}{Point estimate of the log-odds ratio.}
#'     \item{se}{Standard error of the log-odds ratio.}
#'     \item{lo}{Lower bound of the Wald confidence interval.}
#'     \item{hi}{Upper bound of the Wald confidence interval.}
#'   }
#'
#' @references
#' Haldane, J. B. S. (1956). The estimation and significance of the logarithm
#' of a ratio of frequencies. *Annals of Human Genetics*, 20(4), 309-311.
#'
#' @examples
#' # Compare Race1 vs Race2 at IQ bin 1 vs IQ bin 4 (reference)
#' lor_ci_2x2(a = 30, b = 25, c = 28, d = 24)
#'
#' @export
lor_ci_2x2 <- function(a, b, c, d, alpha = 0.05, cc = 0.5) {
  if (any(c(a, b, c, d) == 0)) {
    a <- a + cc; b <- b + cc; c <- c + cc; d <- d + cc
  }
  lor <- log((a * d) / (b * c))
  se  <- sqrt(1/a + 1/b + 1/c + 1/d)
  z   <- qnorm(1 - alpha / 2)
  list(lor = lor, se = se, lo = lor - z * se, hi = lor + z * se)
}


#' Closeness-of-Concordance Measures from OR and CI Endpoints
#'
#' Computes a full row of closeness-of-concordance measures (CCMs) from an
#' odds ratio (OR), its confidence interval endpoints, and the corresponding
#' log-odds ratio (LOR). CCMs include Yule's Q, Yule's Y, and the
#' meta-analytic correlation r_meta (probit transformation of LOR), all on
#' the \[-1, +1\] scale introduced by Kim and Grochowalski (2019).
#'
#' @param OR Numeric. Odds ratio point estimate.
#' @param OR_lo Numeric. Lower confidence limit of the odds ratio.
#' @param OR_hi Numeric. Upper confidence limit of the odds ratio.
#' @param LOR Numeric. Log-odds ratio point estimate.
#' @param LOR_lo Numeric. Lower confidence limit of the log-odds ratio.
#' @param LOR_hi Numeric. Upper confidence limit of the log-odds ratio.
#'
#' @return A one-row `data.frame` with columns:
#'   OR, OR_lo, OR_hi, LOR, LOR_lo, LOR_hi,
#'   YuleQ, Q_lo, Q_hi, YuleY, Y_lo, Y_hi,
#'   r_meta, r_lo, r_hi.
#'
#' @details
#' **Yule's Q** = (OR - 1) / (OR + 1). Ranges from -1 to +1; equals the
#' Pearson correlation for 2x2 tables under a tetrachoric model.
#'
#' **Yule's Y** = (sqrt(OR) - 1) / (sqrt(OR) + 1). A shrunken version of Q
#' with better sampling properties for sparse tables.
#'
#' **r_meta** converts LOR to Cohen's d via d = LOR * sqrt(3) / pi, then to
#' a correlation-like metric via d / sqrt(d^2 + 4). Equivalent to the
#' biserial correlation used in meta-analysis.
#'
#' @references
#' Kim, S.-K., & Grochowalski, J. H. (2019). Gaining from discretization of
#' continuous data: The correspondence analysis biplot approach.
#' *Behavior Research Methods*, 51(2), 589-601.
#' \doi{10.3758/s13428-018-1161-1}
#'
#' @examples
#' lc <- lor_ci_2x2(30, 25, 28, 24)
#' ccm_row(exp(lc$lor), exp(lc$lo), exp(lc$hi), lc$lor, lc$lo, lc$hi)
#'
#' @export
ccm_row <- function(OR, OR_lo, OR_hi, LOR, LOR_lo, LOR_hi) {
  yules_Q   <- function(x) (x - 1) / (x + 1)
  yules_Y   <- function(x) { s <- sqrt(x); (s - 1) / (s + 1) }
  r_meta_fn <- function(l) { d <- l * sqrt(3) / pi; d / sqrt(d^2 + 4) }

  data.frame(
    OR     = OR,     OR_lo  = OR_lo,  OR_hi  = OR_hi,
    LOR    = LOR,    LOR_lo = LOR_lo, LOR_hi = LOR_hi,
    YuleQ  = yules_Q(OR),   Q_lo = yules_Q(OR_lo),  Q_hi = yules_Q(OR_hi),
    YuleY  = yules_Y(OR),   Y_lo = yules_Y(OR_lo),  Y_hi = yules_Y(OR_hi),
    r_meta = r_meta_fn(LOR), r_lo = r_meta_fn(LOR_lo), r_hi = r_meta_fn(LOR_hi)
  )
}
