#' Compute Gower distances between characters
#'
#' Computes Gower distance between characters from a phylogenetic data matrix.
#'
#' @param x A phylogenetic data matrix in Nexus (.nex) format, or in any other
#' data frame or matrix format with a column for each character and terminal
#' taxa as rows, which will be read using [ape::read.nexus.data()]. The data cannot include polymorphisms.
#' @param numeric Whether to treat the values contained in the `x` as
#' numeric or categorical. If `FALSE` (default), features will be
#' considered categorical; if `TRUE`, they will be considered numeric.
#'
#' @returns
#' The Gower distance matrix.
#'
#' @author
#' This function uses code adapted from `StatMatch::gower.dist()`
#' written by Marcello D'Orazio.
#'
#' @seealso
#' `vignette("char-part")` for the use of this function as part
#' of an analysis pipeline.
#'
#' @examples
#' # Load example phylogenetic data matrix
#' data("characters")
#'
#' # Create distance matrix
#' Dmatrix <- get_gower_dist(characters)
#'
#' # Reading data matrix as numeric data
#' Dmatrix <- get_gower_dist(characters, numeric = TRUE)
#'

#' @export
get_gower_dist <- function(x, numeric = FALSE) {
  if (is.matrix(x) || is.data.frame(x) ||
      (is.list(x) && all(lengths(x) == length(x[[1L]])))) {
    data <- as.data.frame(as.matrix(x))
  }
  else if (length(x) == 1L && is.character(x) && endsWith(x, ".nex")) {
    data <- as.data.frame(ape::read.nexus.data(x))
  }
  else {
    stop("'x' must be a data frame, matrix, or file path to a .nex file.", call. = FALSE)
  }

  if (numeric) {
    data[] <- safe_as.numeric(data)
  }

  gowerdist(data)
}

#Port of StatMatch::gower.dist() with default opts but slightly faster
gowerdist <- function(data.x) {
  nx <- nrow(data.x)

  num <- den <- matrix(0, nrow = nx, ncol = nx)

  for (k in seq_len(ncol(data.x))) {
    x <- data.x[, k]
    na.x <- is.na(x)
    delta <- matrix(1, nrow = nx, ncol = nx)
    dist <- matrix(0, nrow = nx, ncol = nx)

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
        zx <- (x - 1) / rng
        dist <- abs(outer(zx, zx, FUN = "-")) / (max(zx[!na.x]) - min(zx[!na.x]))
      }
    }
    else {
      rng <- max(x[!na.x]) - min(x[!na.x])
      if (rng != 0) {
        dist <- abs(outer(x, x, FUN = "-")) / rng
      }
    }

    delta[na.x,] <- 0
    delta[,na.x] <- 0

    n <- dist * delta
    n[is.na(n)] <- 0
    num <- num + n

    den <- den + delta
  }

  out <- num / den
  is.na(out[!is.finite(out)]) <- TRUE

  rownames(out) <- colnames(out) <- rownames(data.x)

  out
}
