get_gower_dist <- function(file, numeric = FALSE) {
  if (is.matrix(file) || is.data.frame(file)) {
    Data_M <- as.data.frame(as.matrix(file))
  }
  else if (is.character(file)) {
    Data_M <- as.data.frame(as.matrix(read.csv(file, colClasses = "character")))
  }
  else stop("'file' must be a file path to a csv file or a data frame.", call. = FALSE)
  
  if (numeric) Data_M[] <- lapply(Data_M, as.numeric)
  
  Dmatrix <- StatMatch::gower.dist(Data_M, KR.corr = TRUE)
  
  return(Dmatrix)
}