#' Simulation results
#'
#' A dataset containing the published simuation results.
#'
#' @format A data frame with 252000 rows and 12 variables:
#' \describe{
#'   \item{set}{Setting; integer from 1 to 36}
#'   \item{sp}{Spacing; character, one of "regular", "partially random", "totally random"}
#'   \item{Ng}{Number of subjects per group; character either "50:50:50" or "40:50:60"}
#'   \item{ni}{Number of observations per subject; integer either 50 or 100}
#'   \item{sigma}{Standard deviation of noise; integer from 1 to 3}
#'   \item{it}{Iteration; integer from 1 to 500}
#'   \item{meth}{Method name; character}
#'   \item{vc}{Variance-covariance structure used in proximity; character, one of "II", "ED", "MD", or NA}
#'   \item{qc}{Cluster quality criteria}
#'   \item{khat}{Estimated number of clusters}
#'   \item{ari}{Adjusted Rand Index}
#'   \item{csa}{Cluster Specific Accuracy}
#' }
"sim_res"
