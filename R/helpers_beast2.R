# attempts to detect type of posterior and corresponding FBD variable names
detect_posterior <- function(df) {
  potential_names <- list(
    MrBayes = c("net_speciation", "relative_extinction", "relative_fossilization"),
    SA = c("diversificationRateFBD", "turnoverFBD", "samplingProportionFBD"),
    BDSky = c("birthRate", "deathRate", "samplingRate"), 
    BDSky_alt = c("netDiversification", "turnOver", "samplingProportion"))
  
  found <- FALSE
  for(xn in 1:length(potential_names)) {
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
  param.names <- sapply(1:length(vidxs), function(v) {
    if(!is.na(vidxs[v])) vnames[vidxs[v]]
    else variables[v]
  })
  param.names
}