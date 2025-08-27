#' Plots distribution of background rates extracted from posterior log files.
#'
#' Plots The distribution and mean of background rates extracted from the
#' posterior log files from Mr. Bayes or BEAST2, as well as the distribution of
#' background rates if log transformed to test for normality of data
#' distribution.
#'
#' @details
#' Plots The distribution and mean (red dotted line) of background rates
#' extracted from the posterior log files from Mr. Bayes or BEAST2, as well as
#' the distribution of background rates if log transformed.  Background rates
#' should be normally distributed for meeting the assumptions of t-tests and
#' other tests passed on by downstream functions, including
#' [get_pwt_rates_MrBayes()], [get_pwt_rates_BEAST2()], and
#' [plot_treerates_sgn()].
#'
#' @param type Whether to use data output from "Mr.Bayes" or "BEAST2".
#' @param posterior A data frame of posterior parameter estimates (log file).
#' From Mr.Bayes, it includes a "clockrate" column indicating the mean
#' (background) clock rate estimate for each generation that will be used for
#' pairwise t-tests. Such data frame can be imported using
#' [combine_log()] (no need to reshape from wide to long). See the
#' [`posterior1p`] or [`posterior3p`] datasets for an
#' examples of how the input file should look.  From BEAST2, it will include at
#' least one "rate<filename>.mean" column indicating the mean (background)
#' clock rate estimate for each generation. If there are "P" unlinked clock
#' partitions in BEAST2, there will be P x "rate<filename>.mean" columns (one
#' for each partition) in the posterior log file.
#' @param clock The clock partition number to calculate selection mode. Ignored
#' if only one clock is available.
#' @param trans Type of data transformation to perform on background rates
#' extracted from the posterior log file from Mr. Bayes or BEAST2. Options
#' include "none" (if rates are normally distributed), natural log
#' transformation "log", and log of base 10 transformation "log10".
#' @param size Font size for title of plot
#' @param quantile Upper limit for X axis (passed on to [ggplot2::xlim()]) to remove
#' outliers from histogram. The quantile can be any value between "0" and "1",
#' but values equal or above "0.95" provide good results in most cases in which
#' the data distribution is right skewed.
#'
#' @returns
#' A `ggplot` object that can be manipulated using
#' \pkg{ggplot2} syntax (e.g., to change the theme or labels).
#'
#' @seealso
#' `vignette("rates-selection")` for the use of this function as
#' part of an analysis pipeline.
#'
#' @examples
#' ## MrBayes example
#' # Load example tree and posterior
#'
#' data("posterior3p")
#'
#' P <- plot_back_rates(type = "MrBayes", posterior3p, clock = 1,
#'                      trans = "log10", size = 10, quantile = 0.95)
#' P

#' @export
plot_back_rates <- function(type = c("MrBayes", "BEAST2"),
                            posterior,
                            clock = 1,
                            trans = c("none", "log", "log10"),
                            size = 12, quantile = 0.95) {

  if (!is.character(type) || length(type) != 1L || !type %in% c("MrBayes", "BEAST2")) {
    stop("Bad type call")
  }

  rates.post <- ..density.. <- clockrate <- Rel.Back.Rate <- NULL

  if (type == "BEAST2") {
    if (missing(posterior) || !is.data.frame(posterior)) {
      stop("'posterior' must be a data frame.", call. = FALSE)
    }

    #get BEAST2 relative background clock rate (for the desired clock partition) and data transform
    posterior.clockrate <- get_clockrate_posterior(posterior)                     #get rate table from posterior sample using 'get_clockrate_posterior' helper
    posterior.clockrate.long <- posterior_clockrate_reshape(posterior.clockrate)  #convert posterior mean rates table to long using 'posterior_clockrate_reshape' helper
    posterior.final <- posterior.clockrate.long[posterior.clockrate.long$clock == clock,]                       #keep values for desired clock

    #check original absolute data distribution
    P1 <- ggplot2::ggplot(posterior.final, aes(x = .data$rates.post)) +
      geom_histogram(aes(y = ..density..), bins = 50, colour="black", fill="white") +
      geom_density(alpha = .2, fill="cyan") +
      geom_vline(aes(xintercept = mean(.data$rates.post)),color="red", linetype="dashed", size=1) +
      labs(x = "Absolute rates", y = "Density") +
      xlim(c(0, quantile(posterior.final$rates.post, quantile))) +
      theme_bw() +
      #labs(title = paste0("Clock partition ", clock), call. = FALSE) +
      theme(plot.title = element_text(size = size, face = "bold", hjust = 0.5))

    if (trans == "none") {
      posterior.rel.clockrate <- posterior.final$clockrate / mean(posterior.final$rates.post)
    }
    else if (trans == "log10"){
      # log10 transform data
      posterior.final$clockrate.log <- log10(posterior.final$rates.post)
    }
    else {
      # ln transform data
      posterior.final$clockrate.log <- log(posterior.final$rates.post)
    }
  }
  else {
    if (missing(posterior) || !is.data.frame(posterior)) {
      stop("'posterior' must be a data frame.", call. = FALSE)
    }

    if (utils::hasName(posterior, "clockrate.all.")) {
      names(posterior)[which(names(posterior) == "clockrate.all.")] <- "clockrate"
    }

    if (!utils::hasName(posterior, "clockrate")) {
      stop("A 'clockrate' column must be present in 'posterior'.", call. = FALSE)
    }

    #get Mr. Bayes relative background clock rate (shared among all partitions)
    ## check original abosulte background rate distribution
    posterior.final <- posterior
    P1 <- ggplot2::ggplot(posterior.final, aes(x = .data$clockrate)) +
      geom_histogram(aes(y = ..density..), bins = 50, colour="black", fill="white")+
      geom_density(alpha=.2, fill="cyan")+
      geom_vline(aes(xintercept = mean(.data$clockrate)), color="red", linetype="dashed", size=1)+
      labs(x = "Absolute rates", y = "Density") +
      xlim(c(0,quantile(posterior.final$clockrate, quantile)))+
      theme_bw() +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

    if (trans == "none") {
      posterior.rel.clockrate <- posterior.final$clockrate / mean(posterior.final$clockrate)
    }
    else if (trans == "log10"){
      # log10 transform data
      posterior.final$clockrate.log <- log10(posterior.final$clockrate)
    }
    else {
      # ln transform data
      posterior.final$clockrate.log <- log(posterior.final$clockrate)
    }
  }

  #get relative background clock rate
  posterior.rel.clockrate <- as.data.frame(posterior.final$clockrate.log / mean(posterior.final$clockrate.log))
  names(posterior.rel.clockrate) <- "Rel.Back.Rate"

  #Plot final rel background rate dist
  if (trans == "none") {
    ## get final relative background rate distribution
    P3 <- ggplot2::ggplot(posterior.rel.clockrate, aes(x = .data$Rel.Back.Rate)) +
      geom_histogram(aes(y = ..density..), bins = 50, colour="black", fill="white")+
      geom_density(alpha = .2, fill="cyan")+
      geom_vline(aes(xintercept=mean(.data$Rel.Back.Rate)),color="red", linetype="dashed", size=1)+
      #xlim(c(0,quantile(posterior.rel.clockrate$Rel.Back.Rate, quantile)))+
      labs(x = "Relative rates", y = "Density") +
      theme_bw() +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

    All_P <- P1 + P3
  }
  else{
    # plot transformed background rate values
    ## transformed absolute rates
    P2 <- ggplot2::ggplot(posterior.final, aes(x = .data$clockrate.log)) +
      geom_histogram(aes(y = ..density..), bins = 50, colour="black", fill="white")+
      geom_density(alpha = .2, fill="cyan")+
      geom_vline(aes(xintercept = mean(.data$clockrate.log)), color="red", linetype="dashed", size=1)+
      labs(x = "Absolute rates (transformed)", y = "Density") +
      theme_bw() +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

    ## transformed relative rates
    P3 <- ggplot2::ggplot(posterior.rel.clockrate, aes(x = .data$Rel.Back.Rate)) +
      geom_histogram(aes(y = ..density..), bins = 50, colour="black", fill="white")+
      geom_density(alpha=.2, fill="cyan")+
      geom_vline(aes(xintercept = mean(.data$Rel.Back.Rate)), color="red", linetype = "dashed", size=1)+
      labs(x = "Relative rates (transformed)", y = "Density") +
      theme_bw() +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

    All_P <- P1 + P2 + P3
  }

  if (type == "BEAST2") {
    Plot <- All_P + patchwork::plot_annotation(
      title = paste0("Background rate distribution for clock partition ", clock)) &
      theme(plot.title = element_text(size = 12, hjust = 0.54, face = "bold"))
  }
  else {
    #Type = "Mr.Bayes" (shared background rates)
    Plot <- All_P + patchwork::plot_annotation(
      title = 'Background rate distribution') &
      theme(plot.title = element_text(size = 12, hjust = 0.54, face = "bold"))
  }

  return(Plot)
}
