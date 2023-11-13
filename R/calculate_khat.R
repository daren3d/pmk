#' Calculate K hat
#'
#' @param pro_mat A N by N matrix of pairwise proximities.
#' @param qc A string specifying the cluster quality criteria.  Choices are "ASW", "CH" and "Gap", for average silhouette width, Calinski criteria and Gap statistic.
#' @param max_k An integer for the largest K.
#'
#' @return An integer for the estimated number of clusters.
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
#' pro_mat <- create.pro_mat(cp, "ED")
#' khat <- calculate.khat(pro_mat, "CH", 6)
calculate.khat <- function(pro_mat, qc, max_k){
  if(qc == "ASW"){
    khat <- rep(NA, max_k)
    for(k in 2:max_k){
      mod <- pam(as.dist(pro_mat), k = k, nstart = 10)
      khat[k] <- mod$silinfo$avg.width
    }
    khat <- which.max(diff(khat[-1]) <= 0) + 1
  }else if(qc == "CH"){
    khat <- rep(NA, max_k)
    N <- nrow(pro_mat)
    R <- sum(pro_mat^2) / N / 2  # total sum of squares (TSS)
    for(k in 2:max_k){
      mod <- pam(as.dist(pro_mat), k = k, nstart = 10)
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
    khat <- which.max(diff(khat[-1]) <= 0) + 1
  }else if(qc == "Gap"){
    cg <- clusGap(pro_mat, pam, K.max = max_k, verbose = FALSE, diss = TRUE)
    khat <- maxSE(cg$Tab[, "gap"], cg$Tab[, "SE.sim"])
  }else{
    stop("Improper `qc` in `calculate.khat`.")
  }
  return(khat)
}
