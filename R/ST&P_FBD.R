#Import and combine log (.p) files from Mr. Bayes, optionally downsampling; produces posterior1p-style object
combine_log <- function(path = ".", burnin = .25, downsample = 1e4) {
  #Get FBD parameter estimates from collection of log files (.p)
  if (!is.character(path)) files <- NULL
  else if (all(utils::file_test(path, op = "-f"))) {
    files <- path
  }
  else if (length(path) == 1 && utils::file_test(path, op = "-d")) {
    files <- paste0(path, "/", list.files(path, pattern = '\\.p$'))
  }
  else files <- NULL
  
  if (length(files) == 0) {
    stop("The value supplied to 'path' must be a character vector containing the names of parameter log (.p) files or of a folder containg such files.", call. = FALSE)
  }
  
  if (length(burnin) != 1 || !is.numeric(burnin) || burnin < 0) {
    stop("'burnin' must be a single number corresponding to the number or percentage of rows to drop from the beginning of each log file.",
         call. = FALSE)
  }
  
  if (length(downsample) != 1 || !is.numeric(downsample) || downsample < 0) {
    stop("'downsample' must be a single number corresponding to the number or percentage of rows to remain after downsampling from each log file.",
         call. = FALSE)
  }
  
  L <- lapply(files, function(x) {
    rtest <- read.table(x, skip = 1, header = TRUE, nrows = 3)
    if (!all(c("Gen", "LnL", "LnPr") %in% names(rtest))) return(NULL)
    
    r <- read.table(x, skip = 1, header = TRUE)
    if (burnin > 0) {
      if (burnin < 1) {
        b <- seq_len(round(burnin * NROW(r)))
      }
      else {
        b <- seq_len(round(min(burnin, NROW(r))))
      }
      r <- r[-b,,drop = FALSE]
    }
    r
  })
  
  if (all(lengths(L) == 0)) {
    stop("No parameter log files were found in the supplied path.", call. = FALSE)
  }
  
  L[lengths(L) == 0] <- NULL
  
  if (!all(vapply(L, function(i) identical(names(i), names(L[[1]])), logical(1L)))) {
    stop("All parameter log files must have the same column names.", call. = FALSE)
  }
  
  samples <- do.call("rbind", L)
  
  if (downsample > 0) {
    if (downsample < 1) {
      d <- as.integer(round(seq(1, NROW(samples), length.out = NROW(samples)*downsample)))
    }
    else {
      d <- as.integer(round(seq(1, NROW(samples), length.out = min(NROW(samples), downsample))))
    }
    samples <- samples[d,,drop = FALSE]
  }
  
  rownames(samples) <- NULL
  
  return(samples)
}

#Reshape AllRuns from wide to long with Time_bins as time and parameters as varying
#' @param variables Names of FBD rate variables in the log. If NULL (default), will attempt to auto-detect the names and log type.
#' @param log.type Name of the software which produced the log (currently supported: MrBayes or BEAST2).
FBD_reshape <- function(posterior, variables = NULL, log.type = c("MrBayes", "BEAST2")) {
  if (!is.data.frame(posterior)) {
    stop("'posterior' must be a data frame.", call. = FALSE)
  }
  if(!is.null(variables)) {
    exist = sapply(names, function(nm) {
      any(startsWith(names(posterior), nm))
    })
    if(any(!exist)) stop("Specified variables not found in posterior")
    
    if(length(log.type) > 1 || !log.type %in% c("MrBayes", "BEAST2")) {
      stop("Log type must be one of 'MrBayes' or 'BEAST2'")
    }
  }
  else {
    autodetect = detect_posterior(posterior)
    variables = autodetect$variables
    log.type = autodetect$log.type
  }
  
  varying = lapply(variables, function(v) {
    names(posterior)[startsWith(names(posterior), v)]
  })
  idname = if(log.type == "MrBayes") "Gen" else "Sample"
  
  posterior_long <- reshape(posterior, direction = "long",
                            varying = varying,
                            v.names = variables,
                            timevar = "Time_bin",
                            sep = if(log.type == "MrBayes") "_" else ".",
                            idvar = "Gen", ids = posterior[[idname]])
  
  posterior_long[["Time_bin"]] <- factor(posterior_long[["Time_bin"]])
  rownames(posterior_long) <- NULL
  attr(posterior_long, "reshapeLong") <- NULL
  
  attr(posterior_long, "log.type") = log.type
  attr(posterior_long, "variables") = variables
  
  posterior_long
}

#Get summary (n, mean, sd, 5 number) of parameters values by time bin
FBD_summary <- function(posterior, file = NULL, digits = 3) {
  
  if (!is.data.frame(posterior) || !"Time_bin" %in% names(posterior) || is.null(attr(posterior, "variables"))) {
    stop("'posterior' must be a data frame of MCMC posterior samples of FBD parameters.", call. = FALSE)
  }
  
  parameters = attr(posterior, "variables")
  
  time.bins <- sort(unique(posterior$Time_bin))
  
  out <- expand.grid(parameter = parameters, Time_bin = time.bins,
                     stringsAsFactors = FALSE)
  
  summary.list <- lapply(seq_len(nrow(out)), function(i) {
    
    in_t <- which(posterior[["Time_bin"]] == out[["Time_bin"]][i])
    
    oneSummary(posterior[[out[["parameter"]][i]]][in_t], digits = digits)
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

#New summary function
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

#Plot density of one parameter by time bin; density or violin plots
FBD_dens_plot <- function(posterior, parameter, type = "density", stack = FALSE, color = "red") {
  
  if (!is.data.frame(posterior) || !"Time_bin" %in% names(posterior) || is.null(attr(posterior, "variables"))) {
    stop("'posterior' must be a data frame of MCMC posterior samples of FBD parameters.", call. = FALSE)
  }
  
  parameters = attr(posterior, "variables")
  
  type <- match.arg(type, c("density", "violin"))
  
  if (missing(parameter)) {
    stop(paste("'parameter' must be one of:", paste0(parameters, collapse = " ")), call. = FALSE)
  }
  parameter <- match.arg(parameter, parameters, several.ok = FALSE)
  
  if(attr(posterior, "log.type") == "MrBayes") {
    param.names <- setNames(gsub("_", " ", firstup(parameters), fixed = TRUE), parameters)
  }
  else param.names = setNames(beast2.names(parameters), parameters)
  
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

#Test assumptions of normality and homoscedasticity for each parameter
#Normality tests within time bin, across time bin, and pooled within time bin
#Homscedasticity tests across time bin
FBD_tests1 <- function(posterior, downsample = TRUE) {
  if (!is.data.frame(posterior) || !"Time_bin" %in% names(posterior) || is.null(attr(posterior, "variables"))) {
    stop("'posterior' must be a data frame of MCMC posterior samples of FBD parameters.", call. = FALSE)
  }
  
  parameters = attr(posterior, "variables")
  
  posterior$Time_bin <- factor(posterior$Time_bin)
  
  #Bartlett test for homogeneity of variance across Time_bins
  #Run for each param
  d <- as.data.frame(matrix(nrow = 1, ncol = 3,
                            dimnames = list(NULL, c("parameter", "statistic", "p-value"))))
  
  bartlett_df <- do.call("rbind", lapply(parameters, function(p) {
    d[["parameter"]] <- p
    bt <- bartlett.test(posterior[[p]], posterior$Time_bin)
    d[["statistic"]] <- bt$statistic
    d[["p-value"]] <- bt$p.value
    return(d)
  }))
  
  #Fligner-Killeen test for homogeneity of variance across Time_bins
  #Run for each param
  fligner_df <- do.call("rbind", lapply(parameters, function(p) {
    d[["parameter"]] <- p
    bt <- fligner.test(posterior[[p]], posterior$Time_bin)
    d[["statistic"]] <- bt$statistic
    d[["p-value"]] <- bt$p.value
    return(d)
  }))
  
  #Shapiro-Wilk test for normality
  #Run within each Time_bin group, overall, and on residuals for each param
  #Need to downsample for SW test
  max.n <- floor(5000/nlevels(posterior$Time_bin))
  
  if (downsample) {
    keep <- unlist(lapply(levels(posterior$Time_bin),
                          function(i) {
                            if (sum(posterior$Time_bin == i) > max.n) {
                              which(posterior$Time_bin == i)[round(seq(1, sum(posterior$Time_bin == i), length.out = max.n))]
                            }
                            else which(posterior$Time_bin == i)
                          }))
    posterior <- posterior[keep,,drop=FALSE]
    run.sw.test <- rep(TRUE, nlevels(posterior$Time_bin))
  }
  else {
    t <- table(posterior$Time_bin)
    run.sw.test <- t[levels(posterior$Time_bin)] <= max.n
    if (!any(run.sw.test)) warning("Shapiro-Wilk normality tests require downsampling and will not be run. Set downsample = TRUE to run these tests.", call. = FALSE)
  }
  
  d <- as.data.frame(matrix(nrow = nlevels(posterior$Time_bin) + 2, ncol = 3,
                            dimnames = list(c(paste("Time bin", levels(posterior$Time_bin)), "Overall", "Residuals"),
                                            c("parameter", "statistic", "p-value"))))
  shapiro_df <- setNames(lapply(parameters, function(p) {
    d[["parameter"]] <- p
    
    for (i in seq_len(nlevels(posterior$Time_bin))[run.sw.test]) {
      sw <- shapiro.test(posterior[[p]][posterior$Time_bin == levels(posterior$Time_bin)[i]])
      d[["statistic"]][i] <- sw[["statistic"]]
      d[["p-value"]][i] <- sw[["p.value"]]
    }
    
    if (all(run.sw.test)) {
      #Overall
      sw <- shapiro.test(posterior[[p]])
      d[["statistic"]][nlevels(posterior$Time_bin) + 1] <- sw[["statistic"]]
      d[["p-value"]][nlevels(posterior$Time_bin) + 1] <- sw[["p.value"]]
      
      #Residuals
      means <- ave(posterior[[p]], posterior$Time_bin)
      sw <- shapiro.test(posterior[[p]] - means)
      d[["statistic"]][nlevels(posterior$Time_bin) + 2] <- sw[["statistic"]]
      d[["p-value"]][nlevels(posterior$Time_bin) + 2] <- sw[["p.value"]]
    }
    
    return(d)
  }), parameters)
  
  list(shapiro = shapiro_df,
       bartlett = bartlett_df,
       fligner = fligner_df)
}

#Test differences in location for each parameter between time bins
FBD_tests2 <- function(posterior, p.adjust.method = "fdr") {
  
  if (!is.data.frame(posterior) || !"Time_bin" %in% names(posterior) || is.null(attr(posterior, "variables"))) {
    stop("'posterior' must be a data frame of MCMC posterior samples of FBD parameters.", call. = FALSE)
  }
  
  parameters = attr(posterior, "variables")
  
  posterior$Time_bin <- factor(posterior$Time_bin)
  
  df_dimnames <- list(NULL, c("parameter", "Time_bin1", "Time_bin2", "n1", "n2", "p-value", "p-value adj"))
  n_tests <- choose(nlevels(posterior[["Time_bin"]]), 2)
  
  d <- as.data.frame(matrix(nrow = n_tests, ncol = length(df_dimnames[[2]]),
                            dimnames = df_dimnames))
  
  #T-tests between pairs of Time_bins for each param
  t_test_list <- setNames(lapply(parameters, function(p) {
    t <- pairwise.t.test(posterior[[p]],
                         posterior[["Time_bin"]],
                         p.adjust.method = "none")
    
    d$parameter <- p
    
    d[paste0("Time_bin", 1:2)] <- as.data.frame(t(combn(levels(posterior[["Time_bin"]]), 2)))
    
    for (i in seq_len(nrow(d))) {
      d[["n1"]][i] <- sum(posterior[["Time_bin"]] == d[["Time_bin1"]][i])
      d[["n2"]][i] <- sum(posterior[["Time_bin"]] == d[["Time_bin2"]][i])
      d[["p-value"]][i] <- t$p.value[d[["Time_bin2"]][i], d[["Time_bin1"]][i]]
    }
    
    d[["p-value adj"]] <- p.adjust(d[["p-value"]], p.adjust.method)
    return(d)
  }), parameters)
  
  #Mann-Whitney U tests between pairs of Time_bins for each param
  mwu_test_list <- setNames(lapply(parameters, function(p) {
    
    t <- pairwise.wilcox.test(posterior[[p]],
                              posterior[["Time_bin"]],
                              p.adjust.method = "none")
    
    d$parameter <- p
    
    d[paste0("Time_bin", 1:2)] <- as.data.frame(t(combn(levels(posterior[["Time_bin"]]), 2)))
    
    for (i in seq_len(nrow(d))) {
      d[["n1"]][i] <- sum(posterior[["Time_bin"]] == d[["Time_bin1"]][i])
      d[["n2"]][i] <- sum(posterior[["Time_bin"]] == d[["Time_bin2"]][i])
      d[["p-value"]][i] <- t$p.value[d[["Time_bin2"]][i], d[["Time_bin1"]][i]]
    }
    
    
    d[["p-value adj"]] <- p.adjust(d[["p-value"]], p.adjust.method)
    return(d)
  }), parameters)
  
  list(t_tests = t_test_list,
       mwu_tests = mwu_test_list)
}

#Visualize deviations from normality for each parameter by time bin using
#density plots
FBD_normality_plot <- function(posterior) {
  
  if (!is.data.frame(posterior) || !"Time_bin" %in% names(posterior) || is.null(attr(posterior, "variables"))) {
    stop("'posterior' must be a data frame of MCMC posterior samples of FBD parameters.", call. = FALSE)
  }
  
  parameters = attr(posterior, "variables")
  
  if(attr(posterior, "log.type") == "MrBayes") {
    param.names <- setNames(gsub("_", " ", firstup(parameters), fixed = TRUE), parameters)
  }
  else param.names = setNames(beast2.names(parameters), parameters)
  
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
  aux_d$x_high <- 3*aux_d$sd
  aux_d$x_low <- -3*aux_d$sd
  
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

