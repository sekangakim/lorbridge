#' SONSCA Column-Isometric Coordinates
#'
#' Extracts row standard coordinates and column principal coordinates from a
#' Singly-Ordered Nonsymmetric Correspondence Analysis (SONSCA) fit via
#' `CAvariants`. This is the column-isometric scaling recommended for
#' cosine theta computation.
#'
#' @param tab A numeric matrix or table. Rows are nominal categories (e.g.,
#'   racial groups); columns are ordered categories (e.g., score bins).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{row_coords}{Row standard coordinate matrix (rows = row categories).}
#'     \item{col_coords}{Column principal coordinate matrix (rows = column categories).}
#'   }
#'
#' @importFrom CAvariants CAvariants
#' @export
sonsca_coords <- function(tab) {
  s <- CAvariants(tab, catype = "SONSCA")
  row_coords <- as.matrix(if (!is.null(s$Rstdcoord))  s$Rstdcoord  else s$Rprinccoord)
  col_coords <- as.matrix(if (!is.null(s$Cprinccoord)) s$Cprinccoord else s$Cstdcoord)
  if (is.null(rownames(row_coords)) || any(rownames(row_coords) == "")) rownames(row_coords) <- rownames(tab)
  if (is.null(rownames(col_coords)) || any(rownames(col_coords) == "")) rownames(col_coords) <- colnames(tab)
  list(row_coords = row_coords, col_coords = col_coords)
}


#' Doubly-Anchored Cosine Theta for SONSCA
#'
#' Computes the matrix of doubly-anchored cosine theta values between all
#' non-anchor row and column pairs in SONSCA coordinate space.
#'
#' @param row_coords Row coordinate matrix (from `sonsca_coords()`).
#' @param col_coords Column coordinate matrix (from `sonsca_coords()`).
#' @param row_anchor Character. Row label used as the anchor (reference).
#' @param col_anchor Character. Column label used as the anchor (reference).
#'
#' @return A matrix of cosine theta values with rows = non-anchor row
#'   categories and columns = non-anchor column categories. Anchor rows/
#'   columns yield `NA` (zero displacement).
#'
#' @export
sonsca_cosines <- function(row_coords, col_coords, row_anchor, col_anchor) {
  stopifnot(row_anchor %in% rownames(row_coords), col_anchor %in% rownames(col_coords))
  F0 <- sweep(row_coords, 2, row_coords[row_anchor, ], "-")
  G0 <- sweep(col_coords, 2, col_coords[col_anchor, ], "-")
  fn <- sqrt(rowSums(F0^2)); gn <- sqrt(rowSums(G0^2))
  COS <- (F0 %*% t(G0)) / outer(fn, gn, "*")
  COS[!is.finite(COS)] <- NA_real_
  dimnames(COS) <- list(rownames(row_coords), rownames(col_coords))
  COS
}


#' Bootstrap Confidence Intervals for SONSCA Cosine Theta
#'
#' Generates bias-corrected and accelerated (BCa-style percentile) bootstrap
#' confidence intervals for doubly-anchored cosine theta estimates from
#' SONSCA. Resamples the contingency table under a multinomial model.
#'
#' @param tab A numeric matrix or table for SONSCA.
#' @param row_anchor Character. Row anchor label.
#' @param col_anchor Character. Column anchor label.
#' @param row_groups Character vector. Non-anchor row labels to include.
#' @param col_groups Character vector. Non-anchor column labels to include.
#' @param B Integer. Number of bootstrap replications (default 2000).
#' @param alpha Numeric. Significance level for CIs (default 0.05).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{point}{Numeric vector of point estimates (flattened row x col).}
#'     \item{lo}{Lower CI bounds.}
#'     \item{hi}{Upper CI bounds.}
#'     \item{boot_n}{Number of successful bootstrap replications.}
#'   }
#'
#' @importFrom stats rmultinom quantile
#' @export
sonsca_bootstrap <- function(tab, row_anchor, col_anchor,
                             row_groups, col_groups,
                             B = 2000, alpha = 0.05) {
  .align_sign <- function(Fb, Gb, Fref, Gref, K) {
    Fb <- Fb[, 1:K, drop = FALSE]; Gb <- Gb[, 1:K, drop = FALSE]
    for (k in seq_len(K)) {
      sgn <- sign(sum(Fb[, k] * Fref[, k], na.rm = TRUE))
      if (is.na(sgn) || sgn == 0) sgn <- 1
      Fb[, k] <- sgn * Fb[, k]; Gb[, k] <- sgn * Gb[, k]
    }
    list(row_coords = Fb, col_coords = Gb)
  }

  base <- sonsca_coords(tab)
  K0   <- min(ncol(base$row_coords), ncol(base$col_coords))
  COS0 <- sonsca_cosines(base$row_coords[, 1:K0, drop = FALSE],
                          base$col_coords[, 1:K0, drop = FALSE],
                          row_anchor, col_anchor)
  t0   <- as.numeric(COS0[row_groups, col_groups, drop = FALSE])
  n    <- sum(tab); p <- as.vector(tab) / n
  I    <- nrow(tab)
  store <- matrix(NA_real_, nrow = length(t0), ncol = B)
  got <- 0; tries <- 0; max_iter <- 100 * B

  while (got < B && tries < max_iter) {
    tries <- tries + 1
    tb <- matrix(rmultinom(1, size = n, prob = p), nrow = I,
                 dimnames = list(rownames(tab), colnames(tab)))
    if (!all(rowSums(tb) > 0) || !all(colSums(tb) > 0)) next
    bc <- tryCatch(sonsca_coords(tb), error = function(e) NULL)
    if (is.null(bc)) next
    if (min(ncol(bc$row_coords), ncol(bc$col_coords)) < K0) next
    al   <- .align_sign(bc$row_coords, bc$col_coords, base$row_coords, base$col_coords, K0)
    COSb <- tryCatch(
      sonsca_cosines(al$row_coords, al$col_coords, row_anchor, col_anchor),
      error = function(e) NULL)
    if (is.null(COSb)) next
    v <- as.numeric(COSb[row_groups, col_groups, drop = FALSE])
    if (!all(is.finite(v))) next
    got <- got + 1; store[, got] <- v
  }
  if (got == 0) stop("No valid bootstrap replicates obtained.")
  store <- store[, seq_len(got), drop = FALSE]
  lo <- apply(store, 1, quantile, probs = alpha / 2,     na.rm = TRUE)
  hi <- apply(store, 1, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  list(point = t0, lo = lo, hi = hi, boot_n = got)
}


#' Pairwise CCMs for SONSCA
#'
#' Computes pairwise closeness-of-concordance measures (CCMs) for a single
#' (row, column) contrast against the anchor (row_anchor, col_anchor) in a
#' SONSCA contingency table.
#'
#' @param tab A numeric matrix or table.
#' @param row_k Character. Focal row label.
#' @param bin_j Character. Focal column label.
#' @param row_anchor Character. Row anchor label.
#' @param col_anchor Character. Column anchor label.
#' @param alpha Numeric. Significance level (default 0.05).
#'
#' @return A one-row `data.frame` with columns: Race, Bin, and all CCM
#'   columns from `ccm_row()`.
#'
#' @export
sonsca_ccm <- function(tab, row_k, bin_j, row_anchor, col_anchor,
                       alpha = 0.05) {
  a <- tab[row_k, bin_j];         b <- tab[row_k, col_anchor]
  c <- tab[row_anchor, bin_j];    d <- tab[row_anchor, col_anchor]
  lc    <- lor_ci_2x2(a, b, c, d, alpha)
  OR_pt <- exp(lc$lor); OR_lo <- exp(lc$lo); OR_hi <- exp(lc$hi)
  data.frame(Race = row_k, Bin = bin_j,
             ccm_row(OR_pt, OR_lo, OR_hi, lc$lor, lc$lo, lc$hi),
             row.names = NULL)
}


#' Percent Inertia from SONSCA
#'
#' Computes the percentage of total inertia explained by each dimension in
#' a SONSCA solution, using a direct SVD of the standardised residual matrix.
#'
#' @param tab A numeric matrix or table.
#'
#' @return Numeric vector of percent inertia values (summing to 100).
#'
#' @export
inertia_pct <- function(tab) {
  P  <- tab / sum(tab); r <- rowSums(P); cv <- colSums(P)
  Dc <- diag(1 / cv, nrow = length(cv))
  Q  <- as.matrix(P) %*% Dc
  Z  <- Q - outer(r, cv)
  d  <- svd(Z)$d; eig <- d^2; 100 * (eig / sum(eig))
}
