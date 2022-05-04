#Helpers
firstup <- function(x) {
  #Capitalize first letter
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#Port of StatMatch::gower.dist() with default opts but slightly faster
gowerdist <- function(data.x) {
  nx <- nrow(data.x)

  num <- den <- matrix(0, nx, nx)

  for (k in seq_len(ncol(data.x))) {
    x <- data.x[, k]
    na.x <- is.na(x)
    delta <- matrix(1, nx, nx)
    dist <- matrix(0, nx, nx)

    if (is.logical(x)) {
      dist[!x,] <- 1
      dist[,!x] <- 1
      delta[!x, !x] <- 0
    }
    else if (is.character(x) || (is.factor(x) && !is.ordered(x))) {
      dist[outer(x, x, FUN = "!=")] <- 1
    }
    else if (is.ordered(x)) {
      x <- as.numeric(x)
      rng <- max(x[!na.x]) - 1
      if (rng != 0) {
        zx <- (x - 1)/rng
        dist <- abs(outer(zx, zx, FUN = "-"))/(max(zx[!na.x]) - min(zx[!na.x]))
      }
    }
    else {
      rng <- max(x[!na.x]) - min(x[!na.x])
      if (rng != 0) dist <- abs(outer(x, x, FUN = "-"))/rng
    }

    delta[na.x,] <- 0
    delta[,na.x] <- 0

    n <- dist * delta
    n[is.na(n)] <- 0
    num <- num + n

    den <- den + delta
  }

  out <- num/den
  is.na(out[!is.finite(out)]) <- TRUE

  rownames(out) <- colnames(out) <- rownames(data.x)

  return(out)
}

safe_as.numeric <- function(x) {
  l <- TRUE
  if (!is.list(x)) {
    l <- FALSE
    x <- list(x)
  }

  non_num <- FALSE
  for (i in seq_along(x)) {
    nas <- is.na(x[[i]])
    suppressWarnings(xn <- as.numeric(x[[i]]))
    if (!non_num && anyNA(xn) && !identical(nas, is.na(xn))) non_num <- TRUE
    x[[i]] <- xn
  }

  if (non_num) {
    warning("NAs introduced by coercion; some values could not be converted to numbers.", call. = FALSE)
  }

  if (!l) x[[1]]
  else x

}
