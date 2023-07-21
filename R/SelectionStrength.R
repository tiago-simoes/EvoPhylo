#Reshape rate_table with a column for each rate to long
clock_reshape <- function(rate_table) {
  rate_table_long <- reshape(rate_table, direction = "long",
                             varying = which(startsWith(names(rate_table), "rates")),
                             v.names = "rates",
                             timevar = "clock",
                             idvar = "nodes",
                             sep = "_")
  rate_table_long[["clock"]] <- factor(rate_table_long[["clock"]])
  rownames(rate_table_long) <- NULL
  attr(rate_table_long, "reshapeLong") <- NULL

  rate_table_long
}

# Check background rates distribution and if they need transformation
plot_back_rates <- function(type = c("MrBayes", "BEAST2"),
                           posterior,
                           clock = 1,
                           trans = c("none", "log", "log10"),
                           size = 12, quantile = 0.95) {

  if(!type %in% c("MrBayes", "BEAST2")) stop("Bad type call")

  if (type == "BEAST2" && (missing(posterior) || !is.data.frame(posterior))) {
      stop("'posterior' must be a data frame.", call. = FALSE)
    }

  if(type == "MrBayes") {
      if (missing(posterior) || !is.data.frame(posterior)) {
        stop("'posterior' must be a data frame.", call. = FALSE)
      }
      if (hasName(posterior, "clockrate.all.")) {
        names(posterior)[which(names(posterior) == "clockrate.all.")] <- "clockrate"
      }
      if (!hasName(posterior, "clockrate")) {
        stop("A 'clockrate' column must be present in 'posterior'.", call. = FALSE)
      }}

   rates.post = ..density.. = clockrate = Rel.Back.Rate = NULL

  if(type == "BEAST2") {
      #get BEAST2 relative background clock rate (for the desired clock partition) and data transform
      posterior.clockrate<-get_clockrate_posterior(posterior)                     #get rate table from posterior sample using 'get_clockrate_posterior' helper
      posterior.clockrate.long<-posterior_clockrate_reshape(posterior.clockrate)  #convert posterior mean rates table to long using 'posterior_clockrate_reshape' helper
      posterior.final<-posterior.clockrate.long[posterior.clockrate.long$clock == clock,]                       #keep values for desired clock

      #check original absolute data distribution
      P1<-ggplot2::ggplot(posterior.final, aes(x=rates.post)) +
        geom_histogram(aes(y=..density..), bins = 50, colour="black", fill="white")+
        geom_density(alpha=.2, fill="cyan")+
        geom_vline(aes(xintercept=mean(rates.post)),color="red", linetype="dashed", size=1)+
        labs(x = "Absolute rates", y = "Density") +
        xlim(c(0,quantile(posterior.final$rates.post, quantile)))+
        theme_bw()+
        #labs(title = paste0("Clock partition ", clock), call. = FALSE) +
        theme(plot.title = element_text(size = size, face = "bold", hjust = 0.5))

      if (trans == "none") {
      #get relative background clock rate
      posterior.rel.clockrate <- posterior.final$rates.post/mean(posterior.final$rates.post)
      }
      else if (trans == "log10"){
      # log10 transform data
      posterior.final$clockrate.log<-log10(posterior.final$rates.post)
      }
      else {
      # ln transform data
      posterior.final$clockrate.log<-log(posterior.final$rates.post)
      }
      #get relative background clock rate
      posterior.rel.clockrate <- as.data.frame(posterior.final$clockrate.log/mean(posterior.final$clockrate.log))
      names(posterior.rel.clockrate) <- "Rel.Back.Rate"
      }

   else {
      #get Mr. Bayes relative background clock rate (shared among all partitions)
      ## check original abosulte background rate distribution
      posterior.final<-posterior
      P1<-ggplot2::ggplot(posterior.final, aes(x=clockrate)) +
        geom_histogram(aes(y=..density..), bins = 50, colour="black", fill="white")+
        geom_density(alpha=.2, fill="cyan")+
        geom_vline(aes(xintercept=mean(clockrate)),color="red", linetype="dashed", size=1)+
        labs(x = "Absolute rates", y = "Density") +
        xlim(c(0,quantile(posterior.final$clockrate, quantile)))+
        theme_bw() +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

      if (trans == "none") {
      posterior.rel.clockrate <- posterior.final$clockrate/mean(posterior.final$clockrate)
      }
      else if (trans == "log10"){
      # log10 transform data
      posterior.final$clockrate.log<-log10(posterior.final$clockrate)
      }
      else {
      # ln transform data
      posterior.final$clockrate.log<-log(posterior.final$clockrate)
      }
      #get relative background clock rate
      posterior.rel.clockrate <- as.data.frame(posterior.final$clockrate.log/mean(posterior.final$clockrate.log))
      names(posterior.rel.clockrate) <- "Rel.Back.Rate"
      }

  #Plot final rel background rate dist
  if(trans == "none"){
      ## get final relative background rate distribution
    P3<- ggplot2::ggplot(posterior.rel.clockrate, aes(x=Rel.Back.Rate)) +
        geom_histogram(aes(y=..density..), bins = 50, colour="black", fill="white")+
        geom_density(alpha=.2, fill="cyan")+
        geom_vline(aes(xintercept=mean(Rel.Back.Rate)),color="red", linetype="dashed", size=1)+
        #xlim(c(0,quantile(posterior.rel.clockrate$Rel.Back.Rate, quantile)))+
        labs(x = "Relative rates", y = "Density") +
        theme_bw() +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

    if(type == "BEAST2") {
    All_P<- P1 + P3
    Plot <- All_P + patchwork::plot_annotation(
                         title = paste0("Background rate distribution for clock partition ", clock)) &
                         theme(plot.title = element_text(size = 12, hjust = 0.54, face = "bold"))
    }
    else {
    #Type = "Mr.Bayes" (shared background rates)
    All_P<- P1 + P3
    Plot <- All_P + patchwork::plot_annotation(
                         title = 'Background rate distribution')&
                         theme(plot.title = element_text(size = 12, hjust = 0.54, face = "bold"))
    }
    }
  else{
    # plot transformed background rate values
    ## transformed absolute rates
    P2<-ggplot2::ggplot(posterior.final, aes(x=clockrate.log)) +
        geom_histogram(aes(y=..density..), bins = 50, colour="black", fill="white")+
        geom_density(alpha=.2, fill="cyan")+
        geom_vline(aes(xintercept=mean(clockrate.log)),color="red", linetype="dashed", size=1)+
        labs(x = "Absolute rates (transformed)", y = "Density") +
        theme_bw() +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

    ## transformed relative rates
    P3<- ggplot2::ggplot(posterior.rel.clockrate, aes(x=Rel.Back.Rate)) +
        geom_histogram(aes(y=..density..), bins = 50, colour="black", fill="white")+
        geom_density(alpha=.2, fill="cyan")+
        geom_vline(aes(xintercept=mean(Rel.Back.Rate)),color="red", linetype="dashed", size=1)+
        labs(x = "Relative rates (transformed)", y = "Density") +
        theme_bw() +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

    if(type == "BEAST2") {
    All_P <- P1 + P2 + P3
    Plot <- All_P + patchwork::plot_annotation(
                         title = paste0("Background rate distribution for clock partition ", clock)) &
                         theme(plot.title = element_text(size = 12, hjust = 0.54, face = "bold"))
    }
    else {
    #Type = "Mr.Bayes" (shared background rates)
    All_P <- P1 + P2 + P3
    Plot <- All_P + patchwork::plot_annotation(
                         title = 'Background rate distribution')&
                         theme(plot.title = element_text(size = 12, hjust = 0.54, face = "bold"))
    }
  }
      return(Plot)
  }

#t-tests with MrBayes output
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
  post.ts <- abs(post.mean - rate_table_long$abs_rate)/post.se
  pvals <- 2*pt(post.ts, df = post.df, lower.tail = FALSE)

  out <- data.frame(rate_table_long$clade,
                    rate_table_long$nodes,
                    rate_table_long$clock,
                    rate_table_long$rate,
                    rate_table_long$abs_rate,
                    pvals)
  names(out) <- c("clade", "nodes", "clock", "relative.rate.mean", "absolute.rate.mean", "p.value")
  return(out)
}


#t-tests with BEAST2 output
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
  clade = nodes = clock = rates = mean.rates.post <- NULL

  posterior.clockrate<-get_clockrate_posterior(posterior)       #get rate table from posterior sample using 'get_clockrate_posterior' helper

  post.df <- nrow(posterior.clockrate) - 1

  #Get mean of posterior means
  post.mean <- as.data.frame(lapply(posterior.clockrate, mean))
  rownames(post.mean)<- "rates"
  post.mean.long <- post_mean_reshape(post.mean)                 #convert posterior mean rates table to long using 'post_mean_reshape' helper

  #Get standard error of posterior means
  post.se <-as.data.frame(lapply(posterior.clockrate, function(x){
    sd(x)/sqrt(length(x))
  }))
  rownames(post.se)<- "se"
  post.se.long <- post_se_reshape(post.se)                       #convert post.se table to long using 'post_se_reshape' helper

  #Get branch rates in long format
  rate.table.long <- clock_reshape(rate_table)

  #Combine all tables by background clock partition
  comb_rates<- Reduce(function(x, y) merge(x, y, all = FALSE, by = "clock"),
                      list(post.mean.long, post.se.long, rate.table.long))
  comb_rates<-subset(comb_rates, select = c(clade, nodes, clock, rates, mean.rates.post, post.se))

  #Calculate relative branch rates
  comb_rates$rel_rate <- comb_rates$rates/comb_rates$mean.rates.post

  #Calculate pairwise t-tests (post.ts) and p-values (pvals)
  post.ts <- abs(comb_rates$mean.rates.post - comb_rates$rate)/comb_rates$post.se
  pvals <- 2*pt(post.ts, df = post.df, lower.tail = FALSE)

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



#Plot tree with colored thresholds

plot_treerates_sgn <- function(type = c("MrBayes", "BEAST2"),
                                 tree, posterior,
                                 trans = c("none", "log", "log10"),
                                 summary = "mean", drop.dummyextant = TRUE,
                                 clock = 1, threshold = c("1 SD", "2 SD"),
                                 low = "blue", mid = "gray90", high = "red",
                                 branch_size = 2, tip_size = 2,
                                 xlim = NULL, nbreaks = 10, geo_size = list(2, 3),
                                 geo_skip = c("Quaternary", "Holocene", "Late Pleistocene")){


  if(!type %in% c("MrBayes", "BEAST2")) stop("Bad type call")

  #Drop extant "dummy" tip
  if (type == "MrBayes" && drop.dummyextant) {
    tree <- treeio::drop.tip(tree, "Dummyextant")
  }

  #Process threshold
  if (length(threshold) > 0) {
    if (type == "BEAST2" && (missing(posterior) || !is.data.frame(posterior))) {
      stop("'posterior' must be a data frame.", call. = FALSE)
    }

    if(type == "MrBayes") {
      if (missing(posterior) || !is.data.frame(posterior)) {
        stop("'posterior' must be a data frame.", call. = FALSE)
      }
      if (hasName(posterior, "clockrate.all.")) {
        names(posterior)[which(names(posterior) == "clockrate.all.")] <- "clockrate"
      }
      if (!hasName(posterior, "clockrate")) {
        stop("A 'clockrate' column must be present in 'posterior'.", call. = FALSE)
      }
    }

    if (!is.character(threshold))
      stop("'threshold' must be a character vector.", call. = FALSE)
    thresh_conf <- endsWith(threshold, "%")
    thresh_sd <- endsWith(tolower(threshold), "sd")

    if (any(!thresh_conf & !thresh_sd)) {
      stop("All entries in 'threshold' must end in '%' for confidence intervals or 'SD' for standard deviations.", call. = FALSE)
    }

    thresh_vals <- character(length(threshold))

    if (any(thresh_conf)) thresh_vals[thresh_conf] <- substring(threshold[thresh_conf], 1, nchar(threshold[thresh_conf]) - 1)
    if (any(thresh_sd)) thresh_vals[thresh_sd] <- substring(threshold[thresh_sd], 1, nchar(threshold[thresh_sd]) - 2)

    if (anyNA(suppressWarnings(as.numeric(thresh_vals)))) {
      stop("All entries in 'threshold' must be confidence levels or the number of standard deviations to use as the thresholds.", call. = FALSE)
    }
    thresh_vals <- as.numeric(thresh_vals)

    if(type == "BEAST2") {
      #get BEAST2 relative background clock rate (for the desired clock partition) and data transform
      posterior.clockrate<-get_clockrate_posterior(posterior)                     #get rate table from posterior sample using 'get_clockrate_posterior' helper
      posterior.clockrate.long<-posterior_clockrate_reshape(posterior.clockrate)  #convert posterior mean rates table to long using 'posterior_clockrate_reshape' helper
      posterior.final<-posterior.clockrate.long[posterior.clockrate.long$clock == clock,]                       #keep values for desired clock

      if (trans == "none") {
      #get relative background clock rate
      posterior.rel.clockrate <- posterior.final$rates.post/mean(posterior.final$rates.post)
      }
      else if (trans == "log10"){
      posterior.final$rates.post.log<-log10(posterior.final$rates.post)
      #get relative background clock rate
      posterior.rel.clockrate <- posterior.final$rates.post.log/mean(posterior.final$rates.post.log)
      }
      else {
      # ln transform data
      posterior.final$rates.post.log<-log(posterior.final$rates.post)
      #get relative background clock rate
      posterior.rel.clockrate <- posterior.final$rates.post.log/mean(posterior.final$rates.post.log)
      }}

     else {
      #get Mr. Bayes relative background clock rate (shared among all partitions)
       if (trans == "none") {
      posterior.rel.clockrate <- posterior$clockrate/mean(posterior$clockrate)
      }
      else if (trans == "log10"){
      posterior$clockrate.log<-log10(posterior$clockrate)
      #get relative background clock rate
      posterior.rel.clockrate <- posterior$clockrate.log/mean(posterior$clockrate.log)
      }
      else {
      # ln transform data
      posterior$clockrate.log<-log(posterior$clockrate)
      #get relative background clock rate
      posterior.rel.clockrate <- posterior$clockrate.log/mean(posterior$clockrate.log)
      }}

    mean.posterior.rel.clockrate <- 1

    breaks <- numeric(2*length(threshold))
    labels <- character(2*length(threshold))

    # Get threshold sd/conf values for each clock partition
    if (any(thresh_sd)) {
      breaks[c(thresh_sd, thresh_sd)] <- mean.posterior.rel.clockrate + c(-thresh_vals[thresh_sd], thresh_vals[thresh_sd]) * sd(posterior.rel.clockrate)
      labels[c(thresh_sd, thresh_sd)] <- c(sprintf("-%s SD", round(thresh_vals[thresh_sd], 2)),
                                           sprintf("+%s SD", round(thresh_vals[thresh_sd], 2)))
    }
    if (any(thresh_conf)) {
      n <- length(posterior.rel.clockrate)
      tcrits <- qt(.5*(1 + thresh_vals[thresh_conf]/100), n - 1)
      breaks[c(thresh_conf, thresh_conf)] <- mean.posterior.rel.clockrate + c(-tcrits, tcrits) * sd(posterior.rel.clockrate)/sqrt(n)
      labels[c(thresh_conf, thresh_conf)] <- c(sprintf("Lower %s%%CI", round(thresh_vals[thresh_conf],2)),
                                               sprintf("Lower %s%%CI", round(thresh_vals[thresh_conf],2)))
    }

    break_order <- order(breaks)
    breaks <- breaks[break_order]
    labels <- labels[break_order]
  }
  else {
    breaks <- numeric(0)
    labels <- character(0)
  }

  if(type == "MrBayes") {
    #Getting multiple clock rates
    p <- unglue::unglue_data(names(tree@data), "rate<model>Brlens<clock>_<summary>",
                             open = "<", close = ">")
    rownames(p) <- names(tree@data)

    p <- p[rowSums(is.na(p)) < ncol(p),,drop=FALSE]

    p$clock <- gsub("\\{|\\}", "", p$clock)

    summary <- match.arg(summary, c("mean", "median"))

    #Process lens
    if (all(p[["clock"]] == "") || all(p[["clock"]] == "all")) {
      if (!isTRUE(clock == 1)) {
        warning("Only one clock is available; ignoring 'clock' value", call. = FALSE)
      }
      if (all(p[["clock"]] == "all")) {
        warning("All data partitions were analyzed as a single clock ('all'); ignoring value in argument 'clock' ", call. = FALSE)
      }
      p[["clock"]] <- 1L
      clock <- 1L
    }
    else if (!is.numeric(clock) || !clock %in% as.numeric(p[["clock"]])) {
      stop(paste0("All 'clock' values must be numeric, but the following values were found in the tree file: ", paste(unique(p[["clock"]]), collapse = ", ")), call. = FALSE)
    }

    rate_var <- rownames(p)[as.numeric(p[["clock"]]) == clock & p[["summary"]] == summary]

    tree@data$`prob+-sd` <- as.factor(tree@data$`prob+-sd`)

    char_col <- vapply(tree@data, is.character, logical(1L))
    tree@data[char_col] <- lapply(tree@data[char_col], as.numeric)

    offset <- min(tree@data$age_median)
  }
  else {
    offset <- min(tree@data$height_median)
  }
  if (type == "MrBayes" && is.null(xlim)) {
    x1 <- -round(max(tree@data$age_median) + 15, -1)
    x2 <- -round(min(tree@data$age_median) - 15, -1)
  }
   else if (type == "BEAST2" && is.null(xlim)) {
    x1 <- -round(max(tree@data$height_median) + 15, -1)
    x2 <- -round(min(tree@data$height_median) - 15, -1)
  }
  else {
    x1 <- round(xlim[1], -1)
    x2 <- round(xlim[2], -1)
  }


  #Create integer version of rate variable split up by breaks
  #MrBayes already uses relative branch rates (normalized)
  if(type == "BEAST2") {
    #get relative branch rates (normalize) and split up by breaks
    tree@data$rel.rate <- tree@data$rate/mean(posterior.final$rates.post)
    tree@data$clockfac <- as.numeric(cut(tree@data$rel.rate, breaks = c(-Inf, breaks, Inf)))
   }
  else {
    tree@data$clockfac <- as.numeric(cut(tree@data[[rate_var]], breaks = c(-Inf, breaks, Inf)))
    }

  #Make tree plot
  selection_plot <- ggtree::ggtree(tree, layout = "rectangular", ladderize = TRUE, right = TRUE,
                                   position = position_nudge(x = -offset),
                                   size = branch_size,
                                   mapping = aes(color = .data$clockfac)) +
    ggtree::geom_tiplab(size = tip_size, linesize = 0.01, fontface = "italic",
                        color = "black", offset = -offset + .5) +
    scale_color_steps2("Background Rate\nThreshold",
                       low = low, high = high,
                       mid = mid,
                       midpoint = mean(c(1, length(breaks) + 1)),
                       breaks = seq_along(breaks) + .5,
                       labels = labels,
                       limits = c(1.5 - 1e-8, length(breaks) + .5 + 1e-8)) +
    deeptime::coord_geo(xlim = c(x1, x2), ylim = c(-1, treeio::Ntip(tree) + 2), expand = FALSE,
                        dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
                        skip = geo_skip,
                        pos = list("bottom", "bottom"), alpha = 1, height = unit(1, "line"),
                        rot = 0, size = geo_size, neg = TRUE) +
    scale_x_continuous(n.breaks = nbreaks, labels = abs) +
    ggtree::theme_tree2() +
    labs(title = sprintf("Selection for partition %s", clock), call. = FALSE) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = c(.05, .25),
          legend.title = element_text(size = 8, face = "bold"),
          legend.key.size = unit(0.5,'cm'),
          legend.text = element_text(size=7))

  selection_plot <- ggtree::revts(selection_plot)
  return(selection_plot)
}
