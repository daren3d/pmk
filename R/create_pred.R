#' Create matix of predictions
#'
#' @param dat `data.frame` with columns `id`, `time`, `y`.
#' @param lo number of predicted points.
#'
#' @return `list`.
#' @export
#'
#' @examples
#' set.seed(808)
#' N <- 50
#' n <- 50
#' dat <- data.frame(id = rep(1:N, each = n),
#' time = rep(seq(0, 1, length.out = n), times = N),
#' y = rnorm(25))
#' cp <- create.pred(dat, 3)
create.pred <- function(dat, lo){
  idd <- unique(dat$id)
  N <- length(idd)

  # Sequence of time points
  t_first <- dat %>%
    group_by(id) %>%
    filter(row_number() == 1) %>%
    pull(time)
  t_last <- dat %>%
    group_by(id) %>%
    filter(row_number() == n()) %>%
    pull(time)
  tps <- seq(max(t_first), min(t_last), length.out = lo)
  s <- data.frame(time = tps)
  # Fit models
  val <- matrix(NA, nrow = N, ncol = length(tps))
  vcm <- vector("list", length = N)
  sigma_est <- rep(NA, N)
  mod0 <- y ~ s(time, bs = "bs")
  for(i in 1:N){
    dat_sub <- dat %>%
      filter(id == idd[i])

    # Assume independence.
    argh4 <- diag(nrow(dat_sub))
    # Use `gam` to set up `magic`.
    mod_gam_F <- gam(mod0, data = dat_sub, fit = FALSE)
    # Fit final model with `magic`.
    mod_magic <- magic(mod_gam_F$y, mod_gam_F$X, mod_gam_F$sp, mod_gam_F$S,
                       mod_gam_F$off, rank = mod_gam_F$rank, C = mod_gam_F$C,
                       w = argh4)
    # It might be necessary to "update" a `gamObject` to inherit methods.
    mod <- magic2gam(mod_gam_F, mod_magic, argh4)
    mod_Xs <- predict(mod, newdata = s, type = "lpmatrix")

    # Save output
    val[i, ] <- c(mod_Xs %*% coefficients(mod))
    vcm[[i]] <- mod_Xs %*% vcov(mod) %*% t(mod_Xs)
    sigma_est[i] <- sqrt(mod$sig2)
  }
  out <- list(val = val, vcm = vcm, N = N,
              tps = tps, sigma_est = mean(sigma_est))
  return(out)
}
