#Helpers

oneSummary <- function(x, digits = 3) {
  qq <- unname(quantile(x, c(0, .25, .5, .75, 1)))
  d <- data.frame(
    n = length(x),
    mean = round(mean(x), digits),
    sd = round(sd(x), digits),
    min = round(qq[1], digits),
    Q1 = round(qq[2], digits),
    median = round(qq[3], digits),
    Q3 = round(qq[4], digits),
    max = round(qq[5], digits)
  )
  rownames(d) <- NULL

  d
}

firstup <- function(x) {
  #Capitalize first letter
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#Reshape AllRuns from wide to long with Time_bins as time and parameters as varying
FBD_reshape <- function(samples) {
  if (!is.data.frame(samples) ||
      !any(startsWith(names(samples), "net_speciation_")) ||
      !any(startsWith(names(samples), "relative_extinction_")) ||
      !any(startsWith(names(samples), "relative_fossilization_"))) {
    stop("'samples' must be a data frame with columns for net_speciation, relative_extinction, and relative_fossilization.", call. = FALSE)
  }

  samples_long <- reshape(samples, direction = "long",
                          varying = list(names(samples)[startsWith(names(samples), "net_speciation_")],
                                         names(samples)[startsWith(names(samples), "relative_extinction_")],
                                         names(samples)[startsWith(names(samples), "relative_fossilization_")]),
                          v.names = c("net_speciation", "relative_extinction", "relative_fossilization"),
                          timevar = "Time_bin",
                          sep = "_",
                          idvar = "Gen", ids = samples[["Gen"]])
  samples_long[["Time_bin"]] <- factor(samples_long[["Time_bin"]])
  attr(samples_long, "reshapeLong") <- NULL
  samples_long
}

#Reshape rate_table with a column for each rate to long
clock_reshape <- function(rate_table) {
  rate_table_long <- reshape(rate_table, direction = "long",
                             varying = which(startsWith(names(rate_table), "rates")),
                             v.names = "rate",
                             timevar = "clock",
                             idvar = "nodes",
                             sep = "_")
  rate_table_long[["clock"]] <- factor(rate_table_long[["clock"]])
  rownames(rate_table_long) <- NULL
  attr(rate_table_long, "reshapeLong") <- NULL

  rate_table_long
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
