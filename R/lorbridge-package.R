#' lorbridge: Bridging Log-Odds Ratios and Correspondence Analysis
#'
#' @description
#' `lorbridge` provides a unified analytical workflow that connects
#' conventional binary and multinomial logistic regression with singly-ordered
#' (SONSCA) and doubly-ordered (DONSCA) nonsymmetric correspondence analysis.
#'
#' Log-odds ratios (LORs) from logistic regression are re-expressed as
#' cosine theta estimates and closeness-of-concordance measures (CCMs) --
#' Yule's Q, Yule's Y, and r_meta -- all on the familiar \[-1, +1\] scale.
#' Bootstrap confidence intervals for cosine theta are provided throughout.
#'
#' @section Main functions:
#' \describe{
#'   \item{`lor_ci_2x2()`}{Haldane-Anscombe-corrected 2x2 LOR with Wald CI.}
#'   \item{`ccm_row()`}{Full CCM row (OR, LOR, YuleQ, YuleY, r_meta) with CIs.}
#'   \item{`blr_continuous()`}{Binary logistic regression, continuous predictor.}
#'   \item{`blr_categorical()`}{Binary logistic regression, categorical predictor.}
#'   \item{`cosine_theta_2row()`}{Cosine theta from 2-row CA via SVD.}
#'   \item{`sonsca_coords()`}{SONSCA column-isometric coordinates.}
#'   \item{`sonsca_cosines()`}{Doubly-anchored cosine theta matrix (SONSCA).}
#'   \item{`sonsca_bootstrap()`}{Bootstrap CIs for SONSCA cosine theta.}
#'   \item{`sonsca_ccm()`}{Pairwise CCMs for SONSCA contrasts.}
#'   \item{`inertia_pct()`}{Percent inertia per dimension.}
#'   \item{`donsca_fit()`}{Fit a DONSCA model.}
#'   \item{`donsca_cosines()`}{Doubly-anchored cosine theta (DONSCA).}
#'   \item{`mlr_ccm()`}{Multinomial logistic regression with CCMs.}
#' }
#'
#' @section Data:
#' \describe{
#'   \item{`lorbridge_data`}{Individual-level dataset (N = 900) with VM
#'     scores, VM bins, and binary minority/majority group membership.}
#'   \item{`tab_IQ`}{4 x 6 contingency table: 4 racial groups x 6 IQ bins.}
#'   \item{`tab_VM`}{4 x 6 contingency table: 4 racial groups x 6 VM bins.}
#'   \item{`tab_IQ_VM`}{6 x 6 contingency table: 6 IQ bins x 6 VM bins
#'     (for DONSCA).}
#' }
#'
#' @references
#' Kim, S.-K., & Grochowalski, J. H. (2019). Gaining from discretization of
#' continuous data: The correspondence analysis biplot approach.
#' *Behavior Research Methods*, 51(2), 589-601.
#' \doi{10.3758/s13428-018-1161-1}
#'
#' Kim, S.-K. (2024). Factorization of person response profiles to identify
#' summative profiles carrying central response patterns.
#' *Psychological Methods*, 29(4), 723-730.
#' \doi{10.1037/met0000568}
#'
#' @keywords internal
"_PACKAGE"


#' Individual-Level VM and Minority Group Dataset
#'
#' An individual-level dataset with N = 900 observations, reconstructed by
#' row expansion from the `vm_raw` grouped data (SONSCA Analysis C). Each row
#' represents one individual with their Vocabulary Meaning (VM) test score,
#' discretised VM bin, binary minority/majority group indicator, and race
#' group label.
#'
#' @format A `data.frame` with 900 rows and 4 columns:
#' \describe{
#'   \item{VM}{Numeric. Raw Vocabulary Meaning score (range 54-149).}
#'   \item{VMbin}{Factor. Discretised VM bin (VM1-VM6), with VM4 as the
#'     reference level. Breakpoints: <=64.28, <=81, <=100, <=121, <=138.36,
#'     >138.36.}
#'   \item{minority}{Integer. Binary outcome: 1 = minority
#'     (Race1 + Race2 + Race3), 0 = majority (Race4).}
#'   \item{Race}{Character. Race group label (Race1, Race2, Race3, Race4).}
#' }
#'
#' @source Reconstructed from the `vm_raw` table in Kim, S.-K. (2026),
#'   unified SONSCA/DONSCA analysis script.
"lorbridge_data"


#' IQ Contingency Table (4 Races x 6 IQ Bins)
#'
#' A 4 x 6 contingency table cross-classifying four racial groups by six
#' discretised IQ score bins. Used in SONSCA Analysis A.
#'
#' @format A numeric matrix with 4 rows (Race1, Race2, Race3, Race4) and
#'   6 columns (IQ1-IQ6).
#'
#' @source Kim, S.-K. (2026), unified SONSCA/DONSCA analysis script,
#'   Analysis A.
"tab_IQ"


#' VM Contingency Table (4 Races x 6 VM Bins)
#'
#' A 4 x 6 contingency table cross-classifying four racial groups by six
#' discretised Vocabulary Meaning (VM) score bins. Used in SONSCA Analysis B.
#'
#' @format A numeric matrix with 4 rows (Race1, Race2, Race3, Race4) and
#'   6 columns (VM1-VM6).
#'
#' @source Kim, S.-K. (2026), unified SONSCA/DONSCA analysis script,
#'   Analysis B.
"tab_VM"


#' IQ x VM Contingency Table (6 IQ Bins x 6 VM Bins)
#'
#' A 6 x 6 contingency table cross-classifying six discretised IQ bins by
#' six discretised Vocabulary Meaning (VM) bins. Both rows and columns are
#' ordered, making this suitable for DONSCA.
#'
#' @format A numeric matrix with 6 rows (IQ1-IQ6) and 6 columns (VM1-VM6).
#'
#' @source Kim, S.-K. (2026), unified SONSCA/DONSCA analysis script,
#'   Analysis 3a.
"tab_IQ_VM"
