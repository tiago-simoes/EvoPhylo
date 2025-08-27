#' Inspect FBD parameter distributions visually
#'
#' Produces plots of the distributions of fossilized birth–death process (FBD)
#' parameters to facilitate the assessment of the assumptions of normality
#' within time bins and homogeneity of variance across time bins.
#'
#' @details
#' The plots produced include density plots for each parameter within each time
#' bin (residualized to have a mean of zero), scaled so that the top of the
#' density is at a value of one (in *black*). Superimposed onto these
#' densities are the densities of a normal distribution with the same mean and
#' variance (and scaled by the same amount) (in *red*). Deviations between
#' the normal density in *red* and the density of the parameters in
#' *black* indicate deviations from normality. The standard deviation of
#' each parameter is also displayed for each time bin to facilitate assessing
#' homogeneity of variance.
#'
#' @inheritParams FBD_dens_plot
#'
#' @returns
#' A `ggplot` object, which can be modified using \pkg{ggplot2}
#' functions.
#'
#' @seealso
#' `vignette("fbd-params")` for the use of this function as part
#' of an analysis pipeline.
#'
#' [combine_log()] for producing a single data set of parameter
#' posterior samples from individual parameter log files.
#'
#' [FBD_reshape()] for converting posterior parameter table from wide
#' to long format.
#'
#' [FBD_tests1()] for statistical tests of normality and homogeneity
#' of variance.
#'
#' [FBD_tests2()] for tests of differences in parameter means.
#'
#' @examples
#' data("posterior3p")
#'
#' posterior3p_long <- FBD_reshape(posterior3p)
#'
#' FBD_normality_plot(posterior3p_long)

#' @export
FBD_normality_plot <- function(posterior) {

  if (!is.data.frame(posterior) || !"Time_bin" %in% names(posterior) || is.null(attr(posterior, "variables"))) {
    stop("'posterior' must be a data frame of MCMC posterior samples of FBD parameters.", call. = FALSE)
  }

  parameters <- attr(posterior, "variables")

  if (attr(posterior, "log.type") == "MrBayes") {
    param.names <- setNames(gsub("_", " ", firstup(parameters), fixed = TRUE), parameters)
  }
  else {
    param.names <- setNames(beast2.names(parameters), parameters)
  }

  posterior <- posterior[c("Time_bin", parameters)]

  time.bins <- sort(unique(posterior$Time_bin))

  posterior$Time_bin <- factor(posterior$Time_bin, levels = time.bins,
                               labels = paste("Time bin", time.bins))

  posterior_long <- reshape(posterior, direction = "long",
                            v.names = "vals",
                            varying = parameters,
                            timevar = "parameter",
                            times = parameters)
  posterior_long$parameter <- factor(posterior_long$parameter, levels = parameters,
                                     labels = param.names[parameters])

  #Turn variables into residualized versions
  posterior_long[["resid_vals"]] <- with(posterior_long, vals - ave(vals, parameter, Time_bin))

  p <- ggplot(data = posterior_long, aes(x = .data$resid_vals)) +
    geom_density(aes(y = after_stat(.data$scaled))) +
    facet_grid(rows = vars(.data$Time_bin), cols = vars(.data$parameter), scales = "free")

  # Add normal densities
  ## Compute SDs for each density (all means = 0)
  aux_d <- aggregate(vals ~ parameter + Time_bin, data = posterior_long,
                     FUN = sd)
  names(aux_d)[names(aux_d) == "vals"] <- "sd"

  # ggplot layer_data uses PANEL for facets
  aux_d$PANEL <- factor(seq_len(nrow(aux_d)))

  ## Expand data range to ensure full normal density is displayed
  aux_d$x_high <- 3 * aux_d$sd
  aux_d$x_low <- -3 * aux_d$sd

  aux_d_long <- reshape(aux_d, direction = "long", timevar = "high_low",
                        times = c("high", "low"), v.names = "x",
                        varying = c("x_high", "x_low"))

  p <- p + geom_blank(data = aux_d_long, aes(x = .data$x))

  ## Extract density data from geom_ensity to correctly scale normal densities
  ggpbd <- merge(layer_data(p, 1), aux_d)

  ## Add normal densities
  ggpbd$norm_density <- dnorm(ggpbd$x, sd = ggpbd$sd)

  ## Add scaling factor that ensure all original densities have a height of 1,
  ## then apply that to normal densities
  maxs <- aggregate(density ~ PANEL, data = ggpbd, FUN = max)
  names(maxs)[names(maxs) == "density"] <- "dens_max"
  ggpbd <- merge(ggpbd, maxs)

  p <- p + geom_line(data = ggpbd, aes(x = .data$x, y = .data$norm_density/.data$dens_max), color = "red")

  ## Annotate with standard deviation

  ### Find which side is more empty, add annotation there
  aux_d$sd_loc_x <- NA_real_
  rel_pos <- .85 #higher numbers mean more to the right
  for (i in levels(aux_d$PANEL)) {

    ranges <- c(min(aux_d$x_low[aux_d$PANEL == i],
                    ggpbd$x[ggpbd$PANEL == i]),
                max(aux_d$x_high[aux_d$PANEL == i],
                    ggpbd$x[ggpbd$PANEL == i]))
    ranges_mid <- mean(ranges)

    if (max(ggpbd$scaled[ggpbd$PANEL == i & ggpbd$x < ranges_mid]) >=
        max(ggpbd$scaled[ggpbd$PANEL == i & ggpbd$x >= ranges_mid])) {
      aux_d$sd_loc_x[aux_d$PANEL == i] <- sum(c(1-rel_pos, rel_pos) * ranges)
    }
    else {
      aux_d$sd_loc_x[aux_d$PANEL == i] <- sum(c(rel_pos, 1-rel_pos) * ranges)
    }

  }

  p <- p + geom_label(data = aux_d, aes(x = .data$sd_loc_x, y = .85, label = paste0("italic(SD) == ", format(sd, digits = 2))),
                      size = 3.5, label.padding	= unit(.15, "lines"), parse = TRUE)

  p + theme_bw() + labs(x = "Residual", y = "Scaled density")
}
