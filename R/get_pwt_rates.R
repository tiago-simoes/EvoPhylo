#' Conduct pairwise t-tests between node rates and clock base rate from a
#' Mr.Bayes or BEAST2 output.
#'
#' Produces a data frame containing the results of 1-sample t-tests for the
#' mean of posterior clock rates against each node's absolute clock rate. `get_pwt_rates_MrBayes()` is to be used with Mr. Bayes output, and `get_pwt_rates_BEAST2()` is to be used with BEAST2 output.
#'
#' @name get_pwt_rates
#'
#' @details
#' These functions first transform relative clock rates to
#' absolute rate values for each node and each clock by multiplying these by
#' the mean posterior clock rate base value. Then, for each node and clock, a
#' one-sample t-test is performed with the null hypothesis that the mean of the
#' posterior clockrates is equal to that node and clock's absolute clock rate.
#'
#' @param rate_table A data frame containing a single "value" column (for all rate values) and one column for the "clock" variable (indicating to which clock partition each rate values refers to), such as from the output of [get_clockrate_table_MrBayes()] or [get_clockrate_table_BEAST2()] with an extra `clade` column added, and followed by [clock_reshape()].
#' @param posterior A data frame of posterior parameter estimates including a
#' "clockrate" column indicating the base of the clock rate estimate for each
#' generation that will be used for pairwise t-tests. Such data frame can be
#' imported using [combine_log()] (no need to reshape from wide to
#' long). See the [`posterior1p`] or [`posterior3p`]
#' datasets for an examples of how the input file should look.
#'
#' @returns
#' A long data frame with one row per node per clock and the following
#' columns:
#' \item{clade}{The name of the clade, taken from the "clade" column of `rate_table`}
#' \item{nodes}{The node number, taken from the "node" column of `rate_table`}
#' \item{clock}{The clock partition number}
#' \item{relative.rate}{The relative mean clock rate per node, taken from the "rates" columns of `rate_table`}
#' \item{absolute.rate(mean)}{The absolute mean clock rate per node; the relative clock rate multiplied by the mean of the posterior clock rates}
#' \item{null}{The absolute clock rate used as the null value in the t-test}
#' \item{p.value}{The p-value of the test comparing the mean of the posterior clockrates to each absolute clockrate}
#'
#' @seealso
#' `vignette("rates-selection")` for the use of this function as
#' part of an analysis pipeline.
#'
#' [combine_log()]
#'
#' @examples
#' # Load example rate table and posterior data sets
#' data("RateTable_Means_3p_Clades")
#' data("posterior3p")
#'
#' get_pwt_rates_MrBayes(RateTable_Means_3p_Clades, posterior3p)
#'
#'  \dontrun{
#' # Load example rate table and posterior data sets
#' RateTable_Means_Clades <- system.file("extdata", "RateTable_Means_Clades.csv", package = "EvoPhylo")
#' RateTable_Means_Clades <- read.csv(RateTable_Means_Clades, header = TRUE)
#'
#' posterior <- system.file("extdata", "Penguins_log.log", package = "EvoPhylo")
#' posterior <- read.table(posterior, header = TRUE)
#'
#' get_pwt_rates_BEAST2(RateTable_Means_Clades, posterior)
#' }
#'

#' @export
#' @rdname get_pwt_rates
get_pwt_rates_MrBayes <- function(rate_table, posterior) {
  if (missing(rate_table) || !is.data.frame(rate_table)) {
    stop("'rate_table' must be a data frame.", call. = FALSE)
  }
  if (!any(startsWith(names(rate_table), "rates"))) {
    stop("'rate_table' must contain \"rates\" columns containing clockrate summaries.", call. = FALSE)
  }
  if (!hasName(rate_table, "clade")) {
    stop("A 'clade' column must be present in 'rate_table'.", call. = FALSE)
  }
  if (missing(posterior) || !is.data.frame(posterior)) {
    stop("'posterior' must be a data frame.", call. = FALSE)
  }
  if (hasName(posterior, "clockrate.all.")) {
    names(posterior)[which(names(posterior) == "clockrate.all.")] <- "clockrate"
  }
  if (!hasName(posterior, "clockrate")) {
    stop("A 'clockrate' column must be present in 'posterior'.", call. = FALSE)
  }

  posterior.clockrate <- posterior$clockrate
  post.df <- length(posterior.clockrate) - 1
  post.mean <- mean(posterior.clockrate)

  rate_table_long <- clock_reshape(rate_table)
  rate_table_long$abs_rate <- rate_table_long$rate * post.mean

  post.se <- sd(posterior.clockrate)/sqrt(length(posterior.clockrate))
  post.ts <- abs(post.mean - rate_table_long$abs_rate) / post.se
  pvals <- 2 * pt(post.ts, df = post.df, lower.tail = FALSE)

  out <- data.frame(rate_table_long$clade,
                    rate_table_long$nodes,
                    rate_table_long$clock,
                    rate_table_long$rate,
                    rate_table_long$abs_rate,
                    pvals)
  names(out) <- c("clade", "nodes", "clock", "relative.rate.mean", "absolute.rate.mean", "p.value")
  return(out)
}

#' @export
#' @rdname get_pwt_rates
get_pwt_rates_BEAST2 <- function(rate_table, posterior) {
  if (missing(rate_table) || !is.data.frame(rate_table)) {
    stop("'rate_table' must be a data frame.", call. = FALSE)
  }
  if (!any(startsWith(names(rate_table), "rates"))) {
    stop("'rate_table' must contain \"rates\" columns containing clockrate summaries.", call. = FALSE)
  }
  if (!hasName(rate_table, "clade")) {
    stop("A 'clade' column must be present in 'rate_table'.", call. = FALSE)
  }
  if (missing(posterior) || !is.data.frame(posterior)) {
    stop("'posterior' must be a data frame.", call. = FALSE)
  }
  if (!any(startsWith(names(posterior), "rate"))) {
    stop("At least one clock 'rate' column must be present in 'posterior'.", call. = FALSE)
  }

  clade <- nodes <- clock <- rates <- mean.rates.post <- NULL

  posterior.clockrate <- get_clockrate_posterior(posterior)       #get rate table from posterior sample using 'get_clockrate_posterior' helper

  post.df <- nrow(posterior.clockrate) - 1

  #Get mean of posterior means
  post.mean <- as.data.frame(lapply(posterior.clockrate, mean))
  rownames(post.mean)<- "rates"
  post.mean.long <- post_mean_reshape(post.mean)                 #convert posterior mean rates table to long using 'post_mean_reshape' helper

  #Get standard error of posterior means
  post.se <- as.data.frame(lapply(posterior.clockrate, function(x) {
    sd(x) / sqrt(length(x))
  }))
  rownames(post.se) <- "se"
  post.se.long <- post_se_reshape(post.se)                       #convert post.se table to long using 'post_se_reshape' helper

  #Get branch rates in long format
  rate.table.long <- clock_reshape(rate_table)

  #Combine all tables by background clock partition
  comb_rates <- Reduce(function(x, y) merge(x, y, all = FALSE, by = "clock"),
                       list(post.mean.long, post.se.long, rate.table.long))
  comb_rates <- subset(comb_rates, select = c(clade, nodes, clock, rates, mean.rates.post, post.se))

  #Calculate relative branch rates
  comb_rates$rel_rate <- comb_rates$rates / comb_rates$mean.rates.post

  #Calculate pairwise t-tests (post.ts) and p-values (pvals)
  post.ts <- abs(comb_rates$mean.rates.post - comb_rates$rate)/comb_rates$post.se
  pvals <- 2 * pt(post.ts, df = post.df, lower.tail = FALSE)

  #Output all to a new table
  out <- data.frame(comb_rates$clade,
                    comb_rates$nodes,
                    comb_rates$clock,
                    comb_rates$mean.rates.post,
                    comb_rates$rates,
                    comb_rates$rel_rate,
                    pvals)
  names(out) <- c("clade", "node", "clock", "background.rate.mean", "absolute.rate.mean", "relative.rate.mean", "p.value")
  return(out)
}
