detect_posterior = function(df) {
  potential_names = list(
    MrBayes = c("net_speciation", "relative_extinction", "relative_fossilization"),
    SA = c("diversificationRateFBD", "turnoverFBD", "samplingProportionFBD"),
    BDSky = c("birthRate", "deathRate", "samplingRate"), 
    BDSky_alt = c("netDiversification", "turnOver", "samplingProportion"))
  
  found = F
  for(xn in 1:length(potential_names)) {
    exist = sapply(potential_names[[xn]], function(nm) {
      any(startsWith(names(df), nm))
    })
    if(all(exist)) {
      if(found) stop("Multiple sets of variables found, please specify log.type and variables directly")
      found = T
      fd_names = potential_names[[xn]]
      fd_typ = if(xn == 1) "MrBayes" else "BEAST2"
    }
  }

  if(found) return(list(variables = fd_names, log.type = fd_typ))
  stop("No sets of variables found, please specify log.type and variables directly")
}