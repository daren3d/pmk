#' Create data to be used in simulation
#'
#' @param sp Spacing between a subject's observations.  Options include, `regular`, `partially random` and `totally random`.
#' @param Ng A string with the number of subjects in each group, separated by colons (:).
#' @param ni Number of observations per subject.
#' @param sigma Standard deviation of the normal error term.
#'
#' @return A `list` with data in long and wide form, true cluster assignments and true number of clusters.
#' @export
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 facet_wrap
#' @importFrom tidyr pivot_wider
#'
#' @examples
#' set.seed(808)
#' cd <- create.data(sigma = 0.1)
#' ggplot(cd$dat, aes(x = time, y = response, group = id)) +
#'   geom_line() +
#'   facet_wrap(~ oracle)
create.data <- function(sp = "regular", Ng = "50:50:50", ni = 50, sigma = 1){
  Ng <- Ng |>
    strsplit(":") |>
    do.call(what = c) |>
    as.numeric()
  G <- length(Ng)
  uid <- 1:sum(Ng)
  N <- length(uid)
  or <- lapply(1:G, function(x){rep(LETTERS[x], each = Ng[x])}) |>
    do.call(what = c)
  t2 <- runif(ni * N, 0, 10 / ni)
  t3 <- c(replicate(N, sort(runif(ni, min = 0, max = 10))))
  dat_oracle <- data.frame(id = uid,
                           oracle = or)
  dat_time <- data.frame(id = rep(uid, each = ni),
                         t0 = rep(seq.int(ni), times = N),
                         t1 = rep(seq(0, 10, length.out = ni), times = N),
                         t2 = t2,
                         t3 = t3)
  sp <- match.arg(sp, c("regular", "partially random", "totally random"))
  if(sp == "regular"){
    dat_time$time <- dat_time$t1
  }else if(sp == "partially random"){
    dat_time$time <- dat_time$t1 + dat_time$t2
  }else{
    dat_time$time <- dat_time$t3
  }
  dat <- dat_time %>%
    left_join(dat_oracle, by = "id")
  n <- nrow(dat_time)
  dat <- dat |>
    lapply(X = 1:G, FUN = create.fun) |>
    do.call(what = rbind) %>%
    mutate(er = rnorm(n, sd = sigma),
           response = f + er)
  dat_wide <- dat %>%
    select(id, t0, response) %>%
    pivot_wider(names_from = t0, values_from = response) %>%
    select(-id)
  return(list(dat = dat, dat_wide = dat_wide, oracle = or, G = G))
}
create.A <- function(t){
  y <- log(t + 0.1, base = 10) + 1
  return(y)
}
create.B <- function(t){
  y <- cos(2 * pi / 10 * t - pi) + 1
  return(y)
}
create.C <- function(t){
  y <- -((t - 5) / 5)^3 + 1
  return(y)
}
create.fun <- function(x, dat){
  g <- LETTERS[x]
  dat <- dat %>%
    filter(oracle == g) %>%
    mutate(f = eval(parse(text = paste0("create.", g, "(time)"))))
  return(dat)
}
