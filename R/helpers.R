#Helpers
firstup <- function(x) {
  #Capitalize first letter
  substr(x, 1L, 1L) <- toupper(substr(x, 1L, 1L))
  x
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

  if (l) x
  else x[[1L]]
}
