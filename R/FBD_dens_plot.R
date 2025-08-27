#' Density plots for each FBD parameter
#'
#' Produces a density or violin plot displaying the distribution of FBD
#' parameter samples by time bin.
#'
#' @details
#' Density plots are produced using
#' [ggplot2::stat_density()], and violin plots
#' are produced using [ggplot2::geom_violin()].
#' On violin plots, a horizontal line indicates the median (of the density),
#' and the black dot indicates the mean.
#'
#' @param posterior A data frame of posterior parameter estimates containing a
#' single "Time_bin" column and one column for each FBD parameter value. Such
#' data frame can be imported using [combine_log()] followed by
#' [FBD_reshape()].
#' @param parameter A string containing the name of an FBD parameter in the
#' data frame; abbreviations allowed.
#' @param type The type of plot; either `"density"` for a density plot or
#' `"violin"` for violin plots. Abbreviations allowed.
#' @param stack When `type = "density"`, whether to produce stacked
#' densities (`TRUE`) or overlapping densities (`FALSE`, the
#' default). Ignored otherwise.
#' @param color When `type = "violin"`, the color of the plotted
#' densities.
#'
#' @returns
#' A `ggplot` object, which can be modified using \pkg{ggplot2}
#' functions.
#'
#' @note
#' When setting `type = "violin"`, a warning may appear saying
#' something like `"In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) : collapsing to unique 'x' values"`. This warning can be ignored.
#'
#' @seealso
#' `vignette("fbd-params")` for the use of this function as part
#' of an analysis pipeline.
#'
#' [ggplot2::stat_density()],
#' [ggplot2::geom_violin()] for the underlying
#' functions to produce the plots.
#'
#' [combine_log()] for producing a single data frame of FBD parameter
#' posterior samples from multiple log files.
#'
#' [FBD_reshape()] for converting a single data frame of FBD
#' parameter estimates, such as those imported using [combine_log()],
#' from wide to long format.
#'
#' [FBD_summary()], [FBD_normality_plot()],
#' [FBD_tests1()], and [FBD_tests2()] for other functions
#' used to summarize and display the distributions of the parameters.
#'
#' @examples
#' data("posterior3p")
#'
#' posterior3p_long <- FBD_reshape(posterior3p)
#'
#' FBD_dens_plot(posterior3p_long, parameter = "net_speciation",
#'               type = "density", stack = FALSE)
#'
#' FBD_dens_plot(posterior3p_long, parameter = "net_speciation",
#'               type = "density", stack = TRUE)
#'
#' FBD_dens_plot(posterior3p_long, parameter = "net_speciation",
#'               type = "violin", color = "red")

#' @export
FBD_dens_plot <- function(posterior, parameter, type = "density", stack = FALSE, color = "red") {

  if (!is.data.frame(posterior) || !"Time_bin" %in% names(posterior) || is.null(attr(posterior, "variables"))) {
    stop("'posterior' must be a data frame of MCMC posterior samples of FBD parameters.", call. = FALSE)
  }

  parameters <- attr(posterior, "variables")

  type <- match.arg(type, c("density", "violin"))

  if (missing(parameter)) {
    stop(paste("'parameter' must be one of:", paste0(parameters, collapse = " ")), call. = FALSE)
  }
  parameter <- match.arg(parameter, parameters, several.ok = FALSE)

  if (attr(posterior, "log.type") == "MrBayes") {
    param.names <- setNames(gsub("_", " ", firstup(parameters), fixed = TRUE), parameters)
  }
  else {
    param.names <- setNames(beast2.names(parameters), parameters)
  }

  posterior <- posterior[c("Time_bin", parameter)]

  time.bins <- sort(unique(posterior$Time_bin))

  posterior$Time_bin <- factor(posterior$Time_bin, levels = time.bins)

  posterior_long <- reshape(posterior, direction = "long",
                            v.names = "vals",
                            varying = parameter,
                            timevar = "parameter",
                            times = parameter)
  posterior_long$parameter <- factor(posterior_long$parameter, levels = parameter,
                                     labels = param.names[parameter])

  if (type == "density") {
    p <- ggplot(data = posterior_long, aes(x = .data$vals , fill = .data$Time_bin)) +
      theme_bw() +
      stat_density(position = if (stack) "stack" else "identity",
                   alpha = if (stack) 1 else .6,
                   outline.type = "both", color = "black") +
      labs(x = "Value", y = "Density", fill = "Time bin", color = "Time bin",
           title = param.names[parameter]) +
      theme(plot.title = element_text(hjust = .5),
            legend.position = "bottom")
  }
  else if (type == "violin") {
    p <- ggplot(data = posterior_long, aes(x = .data$Time_bin, y = .data$vals)) +
      theme_bw() +
      geom_violin(color = color, fill = color,
                  alpha=0.5, draw_quantiles = 0.5, size=0.8, trim = FALSE) +
      stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
      guides(fill = "none", color = "none") +
      labs(x = "Time bin", y = "Value", title = param.names[parameter]) +
      theme(plot.title = element_text(hjust = .5))
  }

  p
}
