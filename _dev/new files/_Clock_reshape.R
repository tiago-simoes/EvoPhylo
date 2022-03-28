clock_reshape <- function(RatesByClade) {
  RatesByClade_long <- reshape(RatesByClade, direction = "long",
                          varying = names(RatesByClade)[startsWith(names(RatesByClade), "rates")],
                          v.names = c("rates"),
                          timevar = "clock",
                          idvar = "nodes",
                          sep = "_")
  RatesByClade_long[["clock"]] <- factor(RatesByClade_long[["clock"]])
  RatesByClade_long
}