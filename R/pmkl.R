#' Prediction model-based k-medoids for longitudinal data
#'
#' @param dat A `data.frame` of longitudinal data in long form  with columns `id`, `time`, `response`.
#' @param lo Integer for the number of points to predict at.
#' @param vc A string specifying the proximity.  Choices are, `II` (individual information), `ED` (Euclidean distance) and `MD` (Mahalanobis distance).
#' @param qc The quality criterion used to estimate the number of clusters.  Choices are `CH` (Calinski-Harabasz), `ASW` (Average Silhouette Width) and `Gap` (Gap statistic).
#' @param max_k The maximum number of clusters to consider.  The method considers 2:max_k.
#' @param clusGap_boot Argument passed to `cluster::clusGap` for the number of boot strap iterations.
#'
#' @return The results from `cluster::pam`; a `pam` object.
#' @export
#' @import dplyr
#' @importFrom magrittr `%>%`
#' @import stats
#'
#' @examples
#' set.seed(808)
#' n <- 50
#' dat <- data.frame(id = rep(1:9, each = n),
#'                   time = c(sort(runif(n)), sort(runif(n)), sort(runif(n)),
#'                            sort(runif(n)), sort(runif(n)), sort(runif(n)),
#'                            sort(runif(n)), sort(runif(n)), sort(runif(n))),
#'                   response = c(rnorm(3 * n, mean = 1, sd = 0.1),
#'                                rnorm(3 * n, mean = 2, sd = 0.1),
#'                                rnorm(3 * n, mean = 3, sd = 0.1)))
#' mod <- pmkl(dat, vc = "ED")
#' table(mod$clustering, c(1, 1, 1, 2, 2, 2, 3, 3, 3))
pmkl <- function (dat, lo = 10, vc = "II", qc = "CH", max_k = 6, clusGap_boot = 100) {
  co <- create.pred(dat, lo)
  prox_mat <- create.prox_mat(co, vc)
  khat <- cal.khat(prox_mat, qc, max_k, clusGap_boot)
  mod <- cluster::pam(as.dist(prox_mat), k = khat, nstart = 10)
  return(mod)
}

create.pred <- function(dat, lo = 10){
  idd <- unique(dat$id)
  N <- length(idd)
  # Sequence of time points
  t_first <- dat %>%
    group_by(id) %>%
    dplyr::filter(row_number() == 1) %>%
    pull(time)
  t_last <- dat %>%
    group_by(id) %>%
    dplyr::filter(row_number() == n()) %>%
    pull(time)
  tps <- seq(max(t_first), min(t_last), length.out = lo)
  s <- data.frame(time = tps)
  # Fit models
  val <- matrix(NA, nrow = N, ncol = length(tps))
  vcm <- vector("list", length = N)
  for(i in 1:N){
    dat_sub <- dat %>%
      dplyr::filter(id == idd[i])
    mod <- mgcv::gam(response ~ s(time, bs = "bs"), data = dat_sub)
    mod_Xs <- predict(mod, newdata = s, type = "lpmatrix")
    val[i, ] <- c(mod_Xs %*% coefficients(mod))
    vcm[[i]] <- mod_Xs %*% vcov(mod) %*% t(mod_Xs)
  }
  out <- list(val = val, vcm = vcm, N = N,
              tps = tps)
  return(out)
}
create.prox_mat <- function(co, vc = "II"){
  vc <- match.arg(vc, c("II", "ED", "MD"))
  if (vc == "ED") {
    prox_mat <- as.matrix(dist(co$val))
  } else if (vc == "MD") {
    scv <- solve(chol(lmf::nearPD(var(co$val))))
    prox_mat <- as.matrix(dist(co$val %*% scv))
  } else {
    val <- co$val
    vcm <- co$vcm
    N <- co$N
    prox_mat <- matrix(0, nrow = N, ncol = N)
    for(i1 in 1:N){
      for(i2 in 1:N){
        if(i1 < i2){
          dd1 <- val[i1, , drop = FALSE] - val[i2, , drop = FALSE]
          dd2 <- solve(lmf::nearPD(vcm[[i1]] + vcm[[i2]]))
          dd3 <- sqrt(dd1 %*% dd2 %*% t(dd1))
          prox_mat[i1, i2] <- dd3
          prox_mat[i2, i1] <- dd3
        }
      }
    }
  }
  return(prox_mat)
}
flm <- function(x){
  tmp1 <- diff(x) <= 0
  if(all(!tmp1)){
    return(length(tmp1))
  }else{
    return(which.max(tmp1))
  }
}
cal.khat <- function(prox_mat, qc = "CH", max_k = 6, clusGap_boot = 100){
  qc <- match.arg(qc, c("ASW", "CH", "Gap"))
  if(qc == "ASW"){
    khat <- rep(NA, max_k)
    for(k in 2:max_k){
      mod <- cluster::pam(as.dist(prox_mat), k = k, nstart = 10)
      khat[k] <- mod$silinfo$avg.width
    }
    khat <- flm(khat[-1]) + 1
  }else if(qc == "CH"){
    khat <- rep(NA, max_k)
    N <- nrow(prox_mat)
    R <- sum(prox_mat^2) / N / 2  # total sum of squares (TSS)
    for(k in 2:max_k){
      mod <- cluster::pam(as.dist(prox_mat), k = k, nstart = 10)
      Rg <- rep(NA, k)
      for(gp in 1:k){
        id_gp <- mod$clustering == gp
        pm_gp <- prox_mat[id_gp, id_gp]
        Rg[gp] <- sum(pm_gp^2) / sum(id_gp) / 2
      }
      W <- sum(Rg)  # total within sum of squares (TWSS)
      B <- R - W  # between cluster sum of squares (BSS)
      khat[k] <- (B / (k - 1)) / (W / (N - k))
    }
    khat <- flm(khat[-1]) + 1
  }else{
    cg <- cluster::clusGap(prox_mat, cluster::pam, K.max = max_k,
                           B = clusGap_boot, verbose = FALSE, diss = TRUE)
    khat <- cluster::maxSE(cg$Tab[, "gap"], cg$Tab[, "SE.sim"])
  }
  return(khat)
}
