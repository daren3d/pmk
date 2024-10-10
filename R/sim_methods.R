#' An iteration of a simulation run
#'
#' @return A `data.frame`.
#' @export
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
  ##
  out <- rbind(pe1, pe2, pe3, pe4, pe5, pe6, pe7, pe8, pe9) %>%
    as.data.frame() %>%
    dplyr::mutate(set = set, sp = sp, Ng = Ng, ni = ni, sigma = sigma, it = sim)
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


tab.pii <- function(set){
  res2 <- sim_res %>%
    dplyr::filter(set == set, meth == "pred") %>%
    dplyr::mutate(vc = factor(vc, levels = c("II", "ED", "MD")),
                  qc = factor(qc, levels = c("CH", "ASW", "Gap"))) %>%
    dplyr::select(-it) %>%
    dplyr::group_by(vc, qc) %>%
    dplyr::summarise(khat_mean = sprintf("%.3f", round(mean(khat), 3)),
                     khat_sd = sprintf("%.3f", round(sd(khat), 3)),
                     ari_mean = sprintf("%.3f", round(mean(ari), 3)),
                     ari_sd = sprintf("%.3f", round(sd(ari), 3)),
                     csa_mean = sprintf("%.3f", round(mean(csa), 3)),
                     csa_sd = sprintf("%.3f", round(sd(csa), 3)),
                     .groups = "drop") %>%
    dplyr::mutate(khat = paste0(khat_mean, " (", khat_sd, ")"),
                  ari = paste0(ari_mean, " (", ari_sd, ")"),
                  csa = paste0(csa_mean, " (", csa_sd, ")"),
                  .keep = "unused") %>%
    dplyr::select(vc, qc, khat, ari, csa)
  return(res2)
}
