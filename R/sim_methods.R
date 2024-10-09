#' An iteration of a simulation run
#'
#' @return A `data.frame`.
#' @export
#' @import dplyr
#' @importFrom magrittr `%>%`
#'
#' @examples
#' \dontrun{
#' clusGap_boot <- 500
#' set.seed(808)
#' sp <- "regular"
#' Ng <- "50:50:50"
#' ni <- 50
#' sigma <- 1
#' iteration()}
iteration <- function(){
  cd <- create.data(sp, Ng, ni, sigma)
  G <- cd$G
  ##
  pe1 <- sim.pmkl(cd, G, max_k, vc = "II", qc = "ASW")
  pe2 <- sim.pmkl(cd, G, max_k, vc = "II", qc = "CH")
  pe3 <- sim.pmkl(cd, G, max_k, vc = "II", qc = "Gap", clusGap_boot)
  pe4 <- sim.pmkl(cd, G, max_k, vc = "ED", qc = "ASW")
  pe5 <- sim.pmkl(cd, G, max_k, vc = "ED", qc = "CH")
  pe6 <- sim.pmkl(cd, G, max_k, vc = "ED", qc = "Gap", clusGap_boot)
  pe7 <- sim.pmkl(cd, G, max_k, vc = "MD", qc = "ASW")
  pe8 <- sim.pmkl(cd, G, max_k, vc = "MD", qc = "CH")
  pe9 <- sim.pmkl(cd, G, max_k, vc = "MD", qc = "Gap", clusGap_boot)
  abe <- sim.abe(cd, G, max_k, "CH")
  mld <- sim.mld(cd, G, "Gapb")
  kml <- sim.kml(cd, G, max_k, "ED")
  mcl <- sim.mclust(cd, G, max_k)
  ##
  out <- rbind(pe1, pe2, pe3, pe4, pe5, pe6, pe7, pe8, pe9,
               abe, mld, kml, mcl) %>%
    as.data.frame() %>%
    mutate(set = set, sp = sp, Ng = Ng, ni = ni, sigma = sigma, it = sim)
  return(out)
}

cal.csa <- function(clus, oracle){
  con_mat <- table(clus, oracle)
  pcM <- combinat::permn(nrow(con_mat))
  con_mat_perm <- lapply(pcM, function(x) con_mat[x, ])
  csa <- sapply(con_mat_perm, function(x) sum(diag(x)) / length(oracle))
  return(max(csa))
}

# prediction model-based k-medoids
sim.pmkl <- function(cd, G, max_k = 6, vc = "II", qc = "CH",
                     clusGap_boot = 100, lo = 10){
  co <- create.pred(cd$dat, lo)
  # create proximities
  prox_mat <- create.prox_mat(co, vc)
  # Outcome 1: number of clusters
  khat <- cal.khat(prox_mat, qc, max_k, clusGap_boot)
  # Outcome 2: adjusted rand index
  mod2 <- cluster::pam(as.dist(prox_mat), k = khat, nstart = 10)
  ari <- mclust::adjustedRandIndex(mod2$clustering, cd$oracle)
  # Outcome 3: cluster specific accuracy
  mod3 <- cluster::pam(as.dist(prox_mat), k = G, nstart = 10)
  csa <- cal.csa(mod3$clustering, cd$oracle)
  #
  return(data.frame(meth = "pred", vc = vc, qc = qc,
                    khat = khat, ari = ari, csa = csa))
}

create.beta <- function(dat){
  idd <- unique(dat$id)
  N <- length(idd)

  # Fit models
  val <- matrix(NA, nrow = N, ncol = 6)
  vcm <- vector("list", length = N)
  for(i in 1:N){
    dat_sub <- dat %>%
      dplyr::filter(id == idd[i])
    y_i <- dat_sub$response
    B_i <- splines::bs(dat_sub$time, knots = c(10/3, 20/3), intercept = FALSE)
    mod <- lm(y_i ~ B_i)
    # Save output
    val[i,] <- coef(mod)
    vcm[[i]] <- vcov(mod)
  }
  out <- list(val = val, vcm = vcm, N = N)
  return(out)
}

# Abraham, regression spline coefficients
sim.abe <- function(cd, G, max_k = 6, qc = "CH", clusGap_boot = 100){
  qc <- match.arg(qc, c("ASW", "CH", "Gap"))
  co <- create.beta(cd$dat)
  val <- co$val
  # Outcome 1: number of clusters
  if(qc == "ASW"){
    khat <- rep(NA, max_k)
    for(k in 2:max_k){
      mod <- kmeans(val, centers = k, nstart = 10)
      khat[k] <- mean(cluster::silhouette(x = mod$cluster, dist = dist(val))[, 3])
    }
    khat <- flm(khat[-1]) + 1
  }else if(qc == "CH"){
    khat <- rep(NA, max_k)
    for(k in 2:max_k){
      mod <- kmeans(val, centers = k, nstart = 10)
      khat[k] <- mod$betweenss / mod$tot.withinss * (cd$cb$N - k) / (k - 1)
    }
    khat <- flm(khat[-1]) + 1
  }else {  # qc == "Gap"
    cg <- cluster::clusGap(val, kmeans, K.max = max_k, B = clusGap_boot, verbose = FALSE)
    khat <- cluster::maxSE(cg$Tab[, "gap"], cg$Tab[, "SE.sim"])
  }
  # Outcome 2: adjusted rand index
  mod2 <- kmeans(val, khat, nstart = 10)
  ari <- mclust::adjustedRandIndex(mod2$cluster, cd$oracle)
  # Outcome 3: cluster specific accuracy
  mod3 <- kmeans(val, G, nstart = 10)
  csa <- cal.csa(mod3$cluster, cd$oracle)
  #
  return(data.frame(meth = "abe", vc = NA, qc = qc,
                    khat = khat, ari = ari, csa = csa))
}

# clusterMLD, long form
sim.mld <- function(cd, G, qc = "Gapb"){
  mod <- with(cd$dat, clusterMLD::LongDataCluster(time, response, id))
  qc <- match.arg(qc, c("CH", "Gapb"))
  if(qc == "CH"){
    # Outcome 1: number of clusters
    khat <- mod$No.CH
    # Outcome 2: adjusted rand index
    clus <- mod$Dat.label %>%
      group_by(id) %>%
      dplyr::filter(row_number() == 1) %>%
      pull(label.CH)
    ari <- mclust::adjustedRandIndex(clus, cd$oracle)
  }else{
    # Outcome 1: number of clusters
    khat <- mod$No.Gapb
    # Outcome 2: adjusted rand index
    clus <- mod$Dat.label %>%
      group_by(id) %>%
      dplyr::filter(row_number() == 1) %>%
      pull(label.Gapb)
    ari <- mclust::adjustedRandIndex(clus, cd$oracle)
  }
  # Outcome 3: cluster specific accuracy
  l2 <- lapply(1:G, function(x){data.frame(id = mod$Cluster.Lists[[G]][[x]],
                                           clus = x)}) |>
    do.call(what = rbind) |>
    arrange(id)
  csa <- cal.csa(l2$clus, cd$oracle)
  #
  return(data.frame(meth = "mld", vc = NA, qc = qc,
                    khat = khat, ari = ari, csa = csa))
}

# kml, wide form
sim.kml <- function(cd, G, max_k = 6, vc = "ED"){
  dat_wide <- cd$dat_wide
  t <- as.numeric(colnames(dat_wide))
  id <- rownames(dat_wide)
  vc <- match.arg(vc, c("ED", "MD"))
  if (vc == "MD") {
    csv <- chol(solve(lmf::nearPD(var(dat_wide))))
    dat_wide <- as.matrix(dat_wide) %*% t(csv)
  }
  tid <- 1:ncol(dat_wide)
  mod <- kml::cld(dat_wide, idAll = id, time = t, timeInData = tid)
  invisible(utils::capture.output(kml::kml(mod, nbClusters = 2:max_k,
                                           nbRedrawing = 1,
                                           parAlgo = kml::parALGO(saveFreq = Inf))))
  # Outcome 1: number of clusters
  ch <- unlist(mod["criterionValuesAsMatrix"])
  names(ch) <- NULL
  khat <- flm(ch) + 1
  # Outcome 2: adjusted rand index
  ari <- mclust::adjustedRandIndex(kml::getClusters(mod, khat), cd$oracle)
  # Outcome 3: cluster specific accuracy
  csa <- cal.csa(kml::getClusters(mod, G), cd$oracle)
  #
  return(data.frame(meth = "kml", vc = vc, qc = "CH",
                    khat = khat, ari = ari, csa = csa))
}

# mclust, regression spline coefficients
sim.mclust <- function(cd, G, max_k = 6){
  co <- create.beta(cd$dat)
  val <- co$val
  # Outcome 1: number of clusters
  invisible(utils::capture.output(mod <- mclust::Mclust(val, G = 1:max_k)))
  khat <- mod$G
  # Outcome 2: adjusted rand index
  ari <- mclust::adjustedRandIndex(mod$classification, cd$oracle)
  # Outcome 3: cluster specific accuracy
  invisible(utils::capture.output(mod3 <- mclust::Mclust(val, G = G)))
  csa <- cal.csa(mod3$classification, cd$oracle)
  #
  return(data.frame(meth = "mclust", vc = NA, qc = "BIC",
                    khat = khat, ari = ari, csa = csa))
}

tab.pii <- function(set){
  res2 <- sim_res %>%
    dplyr::filter(set == set, meth == "pred") %>%
    mutate(vc = factor(vc, levels = c("II", "ED", "MD")),
           qc = factor(qc, levels = c("CH", "ASW", "Gap"))) %>%
    select(-it) %>%
    group_by(vc, qc) %>%
    summarise(khat_mean = sprintf("%.3f", round(mean(khat), 3)),
              khat_sd = sprintf("%.3f", round(sd(khat), 3)),
              ari_mean = sprintf("%.3f", round(mean(ari), 3)),
              ari_sd = sprintf("%.3f", round(sd(ari), 3)),
              csa_mean = sprintf("%.3f", round(mean(csa), 3)),
              csa_sd = sprintf("%.3f", round(sd(csa), 3)),
              .groups = "drop") %>%
    mutate(khat = paste0(khat_mean, " (", khat_sd, ")"),
           ari = paste0(ari_mean, " (", ari_sd, ")"),
           csa = paste0(csa_mean, " (", csa_sd, ")"),
           .keep = "unused") %>%
    select(vc, qc, khat, ari, csa)
  return(res2)
}

tab.rival <- function(set){
  res2_pii <- sim_res %>%
    dplyr::filter(set == set, meth == "pred", vc == "II", qc == "CH")
  res2_riv <- sim_res %>%
    dplyr::filter(set == set, meth != "pred")
  res3 <- rbind(res2_pii, res2_riv) %>%
    mutate(meth = factor(meth, levels = c("pred", "mld", "abe", "mclust", "kml", "fun"))) %>%
    select(-it, -vc, -qc) %>%
    group_by(meth) %>%
    summarise(khat_mean = sprintf("%.3f", round(mean(khat), 3)),
              khat_sd = sprintf("%.3f", round(sd(khat), 3)),
              ari_mean = sprintf("%.3f", round(mean(ari), 3)),
              ari_sd = sprintf("%.3f", round(sd(ari), 3)),
              csa_mean = sprintf("%.3f", round(mean(csa), 3)),
              csa_sd = sprintf("%.3f", round(sd(csa), 3)),
              .groups = "drop") %>%
    mutate(khat = paste0(khat_mean, " (", khat_sd, ")"),
           ari = paste0(ari_mean, " (", ari_sd, ")"),
           csa = paste0(csa_mean, " (", csa_sd, ")"),
           .keep = "unused") %>%
    select(meth, khat, ari, csa)
  return(res3)
}
