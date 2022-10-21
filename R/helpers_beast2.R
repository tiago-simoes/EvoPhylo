# attempts to detect type of posterior and corresponding FBD variable names
detect_posterior <- function(df) {
  potential_names <- list(

    MrBayes = c("net_speciation", "relative_extinction", "relative_fossilization"),
    SA = c("diversificationRateFBD", "turnoverFBD", "samplingProportionFBD"),
    BDSky = c("birthRate", "deathRate", "samplingRate"),
    BDSky_alt = c("netDiversification", "turnOver", "samplingProportion"))

  found <- FALSE
  for(xn in seq_len(length(potential_names))) {

    exist <- sapply(potential_names[[xn]], function(nm) {
      any(startsWith(names(df), nm))
    })
    if(all(exist)) {
      if(found) stop("Multiple sets of variables found, please specify variables and log.type directly")
      found <- TRUE
      fd_names <- potential_names[[xn]]
      fd_typ <- if(xn == 1) "MrBayes" else "BEAST2"
    }
  }

  if(found) return(list(variables = fd_names, log.type = fd_typ))
  stop("No sets of variables found, please specify variables and log.type directly")
}

#attempts to make plot titles from accepted variable names - defaults to variable name if no match found
beast2.names <- function(variables) {

  vs <- c("diversificationRateFBD", "turnoverFBD", "samplingProportionFBD", "birthRate", "deathRate", "samplingRate",
         "netDiversification", "turnOver", "samplingProportion")
  vnames <- c("Diversification rate", "Turnover", "Sampling proportion", "Birth rate", "Death rate", "Sampling rate",
             "Net diversification", "Turnover", "Sampling proportion")
  vidxs <- match(variables, vs)
  param.names <- sapply(seq_len(length(vidxs)), function(v) {
    if(!is.na(vidxs[v])) vnames[vidxs[v]]
    else variables[v]
  })
  param.names
}


#####################################
#Helper functions: rates & selection#
#####################################

# collect partition specific background clock rate information from BEAST2 log files
get_clockrate_posterior <- function(posterior){

  p <- unglue::unglue_data(names(posterior), "rate.<filename>part<clock>.mean",
                           open = "<", close = ">")

  rownames(p) <- names(posterior)

  p <- p[rowSums(is.na(p)) < ncol(p),,drop=FALSE]

  p$clock <- gsub("\\{|\\}", "", p$clock)

  rates <- rownames(p)

  rate_table_posterior <- setNames(data.frame(posterior[rates]),
                        c(paste0("rates", p[rates, "clock"])))

  return(rate_table_posterior)
}


# reshape background clock rate tables
post_mean_reshape <- function(post.mean) {
  post.mean.long <- reshape(post.mean, direction = "long",
                             varying = which(startsWith(names(post.mean), "rates")),
                             v.names = "mean.rates.post",
                             timevar = "clock",
                             sep = "_")
  post.mean.long[["clock"]] <- factor(post.mean.long[["clock"]])
  rownames(post.mean.long) <- NULL
  attr(post.mean.long, "reshapeLong") <- NULL

  return(post.mean.long)
}


# reshape tables with standard error
post_se_reshape <- function(post.mean) {
  post.se.long <- reshape(post.mean, direction = "long",
                             varying = which(startsWith(names(post.mean), "rates")),
                             v.names = "post.se",
                             timevar = "clock",
                             sep = "_")
  post.se.long[["clock"]] <- factor(post.se.long[["clock"]])
  rownames(post.se.long) <- NULL
  attr(post.se.long, "reshapeLong") <- NULL

  return(post.se.long)
}
