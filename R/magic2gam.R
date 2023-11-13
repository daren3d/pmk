#' Update `mgcv::gam` with results from `mgcv::magic`
#'
#' @param gam_F Output of `mgcv::gam`.
#' @param magic Output of `mgcv::magic`.
#' @param w Variance-covariance matrix that necessitated the use of `mgcv::magic`.
#'
#' @return A `mgcv::gamObject`.
#' @export
#'
#' @examples
#' # to do.
magic2gam <- function(gam_F, magic, w){
  mod_gam <- gam(G = gam_F)
  #
  mod_gam$coefficients <- magic$b
  mod_gam$sig2 <- magic$scale
  mod_gam$gcv.ubre <- magic$score
  mod_gam$sp <- magic$sp
  # ? <- magic$sp.full
  mod_gam$rV <- magic$rV  # root of Bayesian covariance matrix
  # ? <- magic$gcv.info  # a long list...
  mod_gam$R <- magic$R
  #
  magic2 <- magic.post.proc(gam_F$X, magic, w)
  mod_gam$Vp <- magic2$Vb  # Bayesian covariance matrix.  Used by plot
  mod_gam$Ve <- magic2$Ve
  mod_gam$hat <- magic2$hat
  mod_gam$edf <- magic2$edf
  #
  magic_fitted <- as.numeric(gam_F$X %*% magic$b)
  mod_gam$fitted.values <- magic_fitted
  mod_gam$residuals <- magic_fitted - mod_gam$y
  return(mod_gam)
}
