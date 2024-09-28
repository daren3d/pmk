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
#' @import cluster
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr n
#' @importFrom dplyr pull
#' @importFrom dplyr row_number
#' @import lmf
#' @importFrom magrittr `%>%`
#' @import mgcv
#' @importFrom stats as.dist coef dist var vcov
#'
#' @examples
#' set.seed(808)
#' N <- 50
#' n <- 50
#' dat <- data.frame(id = rep(1:N, each = n),
#' time = rep(seq(0, 1, length.out = n), times = N),
#' response = rnorm(25))
#' mod <- pmkl(dat)
pmkl <- function (dat, lo = 10, vc = "II", qc = "CH", max_k = 6, clusGap_boot = 100) {
  co <- create.pred(dat, lo)
  prox <- create.pro_mat(co, vc)
  khat <- cal.khat(prox, qc, max_k, clusGap_boot)
  mod <- cluster::pam(as.dist(prox), k = khat, nstart = 10)
  return(mod)
}

create.pred <- function(dat, lo = 10){
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
  for(i in 1:N){
    dat_sub <- dat %>%
      filter(id == idd[i])
    mod <- mgcv::gam(response ~ s(time, bs = "bs"), data = dat_sub)
    mod_Xs <- predict(mod, newdata = s, type = "lpmatrix")
    val[i, ] <- c(mod_Xs %*% coefficients(mod))
    vcm[[i]] <- mod_Xs %*% vcov(mod) %*% t(mod_Xs)
  }
  out <- list(val = val, vcm = vcm, N = N,
              tps = tps)
  return(out)
}
create.pro_mat <- function(co, vc = "II"){
  vc <- match.arg(vc, c("II", "ED", "MD"))
  if (vc == "ED") {
    pro_mat <- as.matrix(dist(co$val))
  } else if (vc == "MD") {
    scv <- solve(chol(lmf::nearPD(var(co$prox))))
    pro_mat <- as.matrix(dist(co$prox %*% scv))  # Mahalanobis
  } else {
    val <- co$val
    vcm <- co$vcm
    N <- co$N
    pro_mat <- matrix(0, nrow = N, ncol = N)
    for(i1 in 1:N){
      for(i2 in 1:N){
        if(i1 < i2){
          dd1 <- val[i1, , drop = FALSE] - val[i2, , drop = FALSE]
          dd2 <- solve(lmf::nearPD(vcm[[i1]] + vcm[[i2]]))
          dd3 <- sqrt(dd1 %*% dd2 %*% t(dd1))
          pro_mat[i1, i2] <- dd3
          pro_mat[i2, i1] <- dd3
        }
      }
    }
  }
  return(pro_mat)
}
flm <- function(x){
  which.max(diff(x) <= 0)
}
cal.khat <- function(pro_mat, qc = "CH", max_k = 6, clusGap_boot = 100){
  qc <- match.arg(qc, c("ASW", "CH", "Gap"))
  if(qc == "ASW"){
    khat <- rep(NA, max_k)
    for(k in 2:max_k){
      mod <- cluster::pam(as.dist(pro_mat), k = k, nstart = 10)
      khat[k] <- mod$silinfo$avg.width
    }
    khat <- flm(khat[-1]) + 1
  }else if(qc == "CH"){
    khat <- rep(NA, max_k)
    N <- nrow(pro_mat)
    R <- sum(pro_mat^2) / N / 2  # total sum of squares (TSS)
    for(k in 2:max_k){
      mod <- cluster::pam(as.dist(pro_mat), k = k, nstart = 10)
      Rg <- rep(NA, k)
      for(gp in 1:k){
        id_gp <- mod$clustering == gp
        pm_gp <- pro_mat[id_gp, id_gp]
        Rg[gp] <- sum(pm_gp^2) / sum(id_gp) / 2
      }
      W <- sum(Rg)  # total within sum of squares (TWSS)
      B <- R - W  # between cluster sum of squares (BSS)
      khat[k] <- (B / (k - 1)) / (W / (N - k))
    }
    khat <- flm(khat[-1]) + 1
  }else{
    cg <- cluster::clusGap(pro_mat, pam, K.max = max_k, B = clusGap_boot,
                           verbose = FALSE, diss = TRUE)
    khat <- cluster::maxSE(cg$Tab[, "gap"], cg$Tab[, "SE.sim"])
  }
  return(khat)
}
