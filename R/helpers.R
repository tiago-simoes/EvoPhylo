#Helpers

oneSummary <- function(x, digits = 3) {
  d <- data.frame(
    n = length(x),
    mean = round(mean(x), digits),
    sd = round(sd(x), digits),
    min = round(min(x), digits),
    Q1 = round(quantile(x, .25), digits),
    median = round(median(x), digits),
    Q3 = round(quantile(x, .75), digits),
    max = round(max(x), digits)
  )
  rownames(d) <- NULL

  d
}

firstup <- function(x) {
  #Capitalize first letter
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
