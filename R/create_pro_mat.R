#' Create prediction matrix
#'
#' @param co A list from `create.pred`.
#' @param vc A string specifying the proximity.  Choices are "II", "ED" and "MD" for individual information, euclidean distance and mahalanobis distance.
#'
#' @return A N by N matrix of pairwise proximity.
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
create.pro_mat <- function(co, vc){
  val <- co$val
  vcm <- co$vcm
  N <- co$N
  if(vc == "II"){
    pro_mat <- matrix(0, nrow = N, ncol = N)
    for(i1 in 1:N){
      val1 <- val[i1, , drop = FALSE]
      vcm1 <- vcm[[i1]]
      for(i2 in 1:N){
        if(i1 < i2){
          dd1 <- val1 - val[i2, , drop = FALSE]
          dd2 <- solve(nearPD(as.matrix(vcm1 + vcm[[i2]])))
          dd3 <- sqrt(dd1 %*% dd2 %*% t(dd1))  # / sum(!is.na(dd1))
          pro_mat[i1, i2] <- dd3
          pro_mat[i2, i1] <- dd3
        }
      }
    }
  }else if(vc == "ED"){
    pro_mat <- as.matrix(dist(val))  # Euclidean
  }else if(vc == "MD"){  # Mahalanobis
    csv <- chol(solve(nearPD(var(val))))  # invert then decompose
    pro_mat <- as.matrix(dist(val %*% t(csv)))
    # scv <- solve(chol(lmf::nearPD(var(val))))  # or decompose then invert?
    # pro_mat <- as.matrix(dist(val %*% scv))
  }else{
    stop("Improper `vc` in `create.pro_mat`.")
  }
  return(pro_mat)
}
