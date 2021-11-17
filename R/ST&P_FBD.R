#Get summary (n, mean, sd, 5 number) of parameters values by time bin
FBD_summary <- function(AllRunsrMelted_MC, file = NULL, digits = 3) {
  time.bins <- sort(unique(AllRunsrMelted_MC$Time_bin))
  parameters <- c("net_speciation", "relative_extinction", "relative_fossilization")

  out <- expand.grid(parameter = parameters, Time_bin = time.bins,
                     stringsAsFactors = FALSE)

  summary.list <- lapply(seq_len(nrow(out)), function(i) {

    in_t <- which(AllRunsrMelted_MC[["Time_bin"]] == out[["Time_bin"]][i])

    oneSummary(AllRunsrMelted_MC[[out[["parameter"]][i]]][in_t], digits = digits)
  })

  out <- cbind(out, do.call("rbind", summary.list))

  out$parameter <- factor(out$parameter, levels = parameters)

  out <- out[with(out, order(parameter, Time_bin)),]

  rownames(out) <- NULL

  if (length(file) > 0) {
    write.csv(out, file = file)
    invisible(out)
  }
  else {
    return(out)
  }
}

#Plot density of one parameter by time bin; density or violin plots
FBD_dens_plot <- function(AllRunsrMelted_MC, parameter, type = "density", stack = FALSE, color = "red") {

  type <- match.arg(type, c("density", "violin"))

  params <- c("net_speciation", "relative_extinction",
              "relative_fossilization")
  param.names <- setNames(gsub("_", " ", firstup(params), fixed = TRUE),
                          params)

  if (missing(parameter)) {
    stop("'parameter' must be specified.", call. = FALSE)
  }
  parameter <- match.arg(parameter, params, several.ok = FALSE)

  AllRunsrMelted_MC <- AllRunsrMelted_MC[c("Time_bin", parameter)]

  time.bins <- sort(unique(AllRunsrMelted_MC$Time_bin))

  AllRunsrMelted_MC$Time_bin <- factor(AllRunsrMelted_MC$Time_bin, levels = time.bins)

  AllRunsrMelted_MC_long <- reshape(AllRunsrMelted_MC, direction = "long",
                                    v.names = "vals",
                                    varying = parameter,
                                    timevar = "parameter",
                                    times = parameter)
  AllRunsrMelted_MC_long$parameter <- factor(AllRunsrMelted_MC_long$parameter, levels = parameter,
                                             labels = param.names[parameter])

  if (type == "density") {
    p <- ggplot(AllRunsrMelted_MC_long, aes(x = vals , fill = Time_bin)) + theme_bw() +
      stat_density(position = if (stack) "stack" else "identity",
                   alpha = if (stack) 1 else .6,
                   outline.type = "both", color = "black") +
      labs(x = "Value", y = "Density", fill = "Time bin", color = "Time bin",
           title = param.names[parameter]) +
      theme(plot.title = element_text(hjust = .5),
            legend.position = "bottom")
  }
  else if (type == "violin") {
    p <- ggplot(AllRunsrMelted_MC_long, aes(x = Time_bin, y = vals)) + theme_bw() +
      geom_violin(color = color, fill = color,
                  alpha=0.5, draw_quantiles = 0.5, size=0.8, trim = FALSE) +
      stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
      guides(fill = "none", color = "none") +
      labs(x = "Time bin", y = "Value", title = param.names[parameter]) +
      theme(plot.title = element_text(hjust = .5))
  }

  p
}

#Reshape AllRuns from wide to long with Time_bins as time and parameters as varying
FBD_reshape <- function(AllRuns) {
  AllRuns_long <- reshape(AllRuns, direction = "long",
                          varying = list(names(AllRuns)[startsWith(names(AllRuns), "net_speciation_")],
                                         names(AllRuns)[startsWith(names(AllRuns), "relative_extinction_")],
                                         names(AllRuns)[startsWith(names(AllRuns), "relative_fossilization_")]),
                          v.names = c("net_speciation", "relative_extinction", "relative_fossilization"),
                          timevar = "Time_bin",
                          sep = "_")
  AllRuns_long[["Time_bin"]] <- factor(AllRuns_long[["Time_bin"]])
  AllRuns_long
}

#Test assumptions of normality and homoscedasticity for each parameter
#Normality tests within time bin, across time bin, and pooled within time bin
#Homscedasticity tests across time bin
FBD_tests1 <- function(AllRunsrMelted_MC, downsample = TRUE) {
  params <- c("net_speciation", "relative_extinction", "relative_fossilization")

  AllRunsrMelted_MC$Time_bin <- factor(AllRunsrMelted_MC$Time_bin)

  #Bartlett test for homogeneity of variance across Time_bins
  #Run for each param
  bartlett_df <- do.call("rbind", lapply(params, function(p) {
    d <- as.data.frame(matrix(nrow = 1, ncol = 3,
                              dimnames = list(NULL, c("parameter", "statistic", "p-value"))))
    d[["parameter"]] <- p
    bt <- bartlett.test(AllRunsrMelted_MC[[p]], AllRunsrMelted_MC$Time_bin)
    d[["statistic"]] <- bt$statistic
    d[["p-value"]] <- bt$p.value
    return(d)
  }))

  #Fligner-Killeen test for homogeneity of variance across Time_bins
  #Run for each param
  fligner_df <- do.call("rbind", lapply(params, function(p) {
    d <- as.data.frame(matrix(nrow = 1, ncol = 3,
                              dimnames = list(NULL, c("parameter", "statistic", "p-value"))))
    d[["parameter"]] <- p
    bt <- fligner.test(AllRunsrMelted_MC[[p]], AllRunsrMelted_MC$Time_bin)
    d[["statistic"]] <- bt$statistic
    d[["p-value"]] <- bt$p.value
    return(d)
  }))

  #Shapiro-Wilk test for normality
  #Run within each Time_bin group, overall, and on residuals for each param
  #Need to downsample for SW test
  max.n <- floor(5000/nlevels(AllRunsrMelted_MC$Time_bin))

  if (downsample) {
    keep <- unlist(lapply(levels(AllRunsrMelted_MC$Time_bin),
                          function(i) {
                            if (sum(AllRunsrMelted_MC$Time_bin == i) > max.n) {
                              which(AllRunsrMelted_MC$Time_bin == i)[round(seq(1, sum(AllRunsrMelted_MC$Time_bin == i), length.out = max.n))]
                            }
                            else which(AllRunsrMelted_MC$Time_bin == i)
                          }))
    AllRunsrMelted_MC <- AllRunsrMelted_MC[keep,,drop=FALSE]
    run.sw.test <- rep(TRUE, nlevels(AllRunsrMelted_MC$Time_bin))
  }
  else {
    t <- table(AllRunsrMelted_MC$Time_bin)
    run.sw.test <- t[levels(AllRunsrMelted_MC$Time_bin)] <= max.n
    if (!any(run.sw.test)) warning("Shapiro-Wilk normality tests requiring downsampling and will not be run. Set downsample = TRUE to run these tests.", call. = FALSE)
  }

  shapiro_df <- setNames(lapply(params, function(p) {
    d <- as.data.frame(matrix(nrow = nlevels(AllRunsrMelted_MC$Time_bin) + 2, ncol = 3,
                              dimnames = list(c(paste("Time bin", levels(AllRunsrMelted_MC$Time_bin)), "Overall", "Residuals"),
                                              c("parameter", "statistic", "p-value"))))
    d[["parameter"]] <- p

    for (i in seq_len(nlevels(AllRunsrMelted_MC$Time_bin))[run.sw.test]) {
      sw <- shapiro.test(AllRunsrMelted_MC[[p]][AllRunsrMelted_MC$Time_bin == levels(AllRunsrMelted_MC$Time_bin)[i]])
      d[["statistic"]][i] <- sw[["statistic"]]
      d[["p-value"]][i] <- sw[["p.value"]]
    }

    if (all(run.sw.test)) {
      #Overall
      sw <- shapiro.test(AllRunsrMelted_MC[[p]])
      d[["statistic"]][nlevels(AllRunsrMelted_MC$Time_bin) + 1] <- sw[["statistic"]]
      d[["p-value"]][nlevels(AllRunsrMelted_MC$Time_bin) + 1] <- sw[["p.value"]]

      #Residuals
      means <- ave(AllRunsrMelted_MC[[p]], AllRunsrMelted_MC$Time_bin)
      sw <- shapiro.test(AllRunsrMelted_MC[[p]] - means)
      d[["statistic"]][nlevels(AllRunsrMelted_MC$Time_bin) + 2] <- sw[["statistic"]]
      d[["p-value"]][nlevels(AllRunsrMelted_MC$Time_bin) + 2] <- sw[["p.value"]]
    }

    return(d)
  }), params)

  list(shapiro = shapiro_df,
       bartlett = bartlett_df,
       fligner = fligner_df)
}

#Test differences in location for each parameter between time bins
#or between analysis types for each time bin
FBD_tests2 <- function(AllRunsrMelted_MC) {

  params <- c("net_speciation", "relative_extinction", "relative_fossilization")

  AllRunsrMelted_MC$Time_bin <- factor(AllRunsrMelted_MC$Time_bin)

  df_dimnames <- list(NULL, c("parameter", "Time_bin1", "Time_bin2", "n1", "n2", "p-value", "p-value adj"))
  n_tests <- choose(nlevels(AllRunsrMelted_MC[["Time_bin"]]), 2)

  #T-tests between pairs of Time_bins for each param
  t_test_list <- setNames(lapply(params, function(p) {
      t <- pairwise.t.test(AllRunsrMelted_MC[[p]],
                           AllRunsrMelted_MC[["Time_bin"]],
                           p.adjust.method = "none")
      d <- as.data.frame(matrix(nrow = n_tests, ncol = length(df_dimnames[[2]]),
                                dimnames = df_dimnames))
      d$parameter <- p

      d[paste0("Time_bin", 1:2)] <- as.data.frame(t(combn(levels(AllRunsrMelted_MC[["Time_bin"]]), 2)))

      for (i in seq_len(nrow(d))) {
        d[["n1"]][i] <- sum(AllRunsrMelted_MC[["Time_bin"]] == d[["Time_bin1"]][i])
        d[["n2"]][i] <- sum(AllRunsrMelted_MC[["Time_bin"]] == d[["Time_bin2"]][i])
        d[["p-value"]][i] <- t$p.value[d[["Time_bin2"]][i], d[["Time_bin1"]][i]]
      }

    d[["p-value adj"]] <- p.adjust(d[["p-value"]], "fdr")
    return(d)
  }), params)

  #Mann-Whitney U tests between pairs of Time_bins for each param
  mwu_test_list <- setNames(lapply(params, function(p) {

      t <- pairwise.wilcox.test(AllRunsrMelted_MC[[p]],
                                AllRunsrMelted_MC[["Time_bin"]],
                                p.adjust.method = "none")
      d <- as.data.frame(matrix(nrow = n_tests, ncol = length(df_dimnames[[2]]),
                                dimnames = df_dimnames))
      d$parameter <- p

      d[paste0("Time_bin", 1:2)] <- as.data.frame(t(combn(levels(AllRunsrMelted_MC[["Time_bin"]]), 2)))

      for (i in seq_len(nrow(d))) {
        d[["n1"]][i] <- sum(AllRunsrMelted_MC[["Time_bin"]] == d[["Time_bin1"]][i])
        d[["n2"]][i] <- sum(AllRunsrMelted_MC[["Time_bin"]] == d[["Time_bin2"]][i])
        d[["p-value"]][i] <- t$p.value[d[["Time_bin2"]][i], d[["Time_bin1"]][i]]
      }


    d[["p-value adj"]] <- p.adjust(d[["p-value"]], "fdr")
    return(df)
  }), params)

  list(t_tests = t_test_list,
       mwu_tests = mwu_test_list)
}

#Visualize deviations from normality for each parameter by time bin using
#density plots
FBD_normality_plot <- function(AllRunsrMelted_MC) {
  params <- c("net_speciation", "relative_extinction", "relative_fossilization")

  param.names <- setNames(gsub("_", " ", firstup(params), fixed = TRUE),
                          params)

  AllRunsrMelted_MC <- AllRunsrMelted_MC[c("Time_bin", params)]

  time.bins <- sort(unique(AllRunsrMelted_MC$Time_bin))

  AllRunsrMelted_MC$Time_bin <- factor(AllRunsrMelted_MC$Time_bin, levels = time.bins,
                                       labels = paste("Time bin", time.bins))

  AllRunsrMelted_MC_long <- reshape(AllRunsrMelted_MC, direction = "long",
                                    v.names = "vals",
                                    varying = params,
                                    timevar = "parameter",
                                    times = params)
  AllRunsrMelted_MC_long$parameter <- factor(AllRunsrMelted_MC_long$parameter, levels = params,
                                             labels = param.names[params])

  #Turn variables into residualized versions
  AllRunsrMelted_MC_long[["resid_vals"]] <- with(AllRunsrMelted_MC_long, vals - ave(vals, parameter, Time_bin))

  p <- ggplot(AllRunsrMelted_MC_long, aes(x = resid_vals)) +
    geom_density(aes(y=after_stat(scaled))) +
    facet_grid(rows = vars(Time_bin), cols = vars(parameter), scales = "free")

  # Add normal densities
  ## Compute SDs for each density (all means = 0)
  aux_d <- aggregate(vals ~ parameter + Time_bin, data = AllRunsrMelted_MC_long,
                     FUN = sd)
  names(aux_d)[names(aux_d) == "vals"] <- "sd"

  # ggplot layer_data uses PANEL for facets
  aux_d$PANEL <- factor(seq_len(nrow(aux_d)))

  ## Expand data range to ensure full normal density is displayed
  aux_d$x_high <- 3*aux_d$sd
  aux_d$x_low <- -3*aux_d$sd

  aux_d_long <- reshape(aux_d, direction = "long", timevar = "high_low",
                        times = c("high", "low"), v.names = "x",
                        varying = c("x_high", "x_low"))

  p <- p + geom_blank(data = aux_d_long, aes(x = x))

  ## Extract density data from geom_ensity to correctly scale normal densities
  ggpbd <- merge(layer_data(p, 1), aux_d)

  ## Add normal densities
  ggpbd$norm_density <- dnorm(ggpbd$x, sd = ggpbd$sd)

  ## Add scaling factor that ensure all original densities have a height of 1,
  ## then apply that to normal densities
  maxs <- aggregate(density ~ PANEL, data = ggpbd, FUN = max)
  names(maxs)[names(maxs) == "density"] <- "dens_max"
  ggpbd <- merge(ggpbd, maxs)

  p <- p + geom_line(data = ggpbd, aes(x = x, y = norm_density/dens_max), color = "red")

  ## Annotate with standard deviation

  ### Find which side is more empty, add annotation there
  aux_d$sd_loc_x <- NA_real_
  for (i in levels(aux_d$PANEL)) {

    ranges <- c(min(aux_d$x_low[aux_d$PANEL == i], ggpbd$x[ggpbd$PANEL == i]), max(aux_d$x_high[aux_d$PANEL == i], ggpbd$x[ggpbd$PANEL == i]))
    ranges_mid <- mean(ranges)

    if (max(ggpbd$scaled[ggpbd$PANEL == i & ggpbd$x < ranges_mid]) >= max(ggpbd$scaled[ggpbd$PANEL == i & ggpbd$x >= ranges_mid])) {
      aux_d$sd_loc_x[aux_d$PANEL == i] <- sum(c(.25, .75) * ranges)
    }
    else {
      aux_d$sd_loc_x[aux_d$PANEL == i] <- sum(c(.75, .25) * ranges)
    }

  }

  p <- p + geom_label(data = aux_d, aes(x = sd_loc_x, y = .8, label = paste0("italic(SD) == ", format(sd, digits = 2))), parse = TRUE)


  p <- p + theme_bw() +
    labs(x = "Residual", y = "Scaled density")
  p
}

