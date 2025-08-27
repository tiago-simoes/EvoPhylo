#' Write character partitions as separate Nexus files (for use in BEAUti)
#'
#' @param x character data matrix as Nexus file (.nex) or data frame (with taxa as rows and characters as columns) read directly from local directory
#' @param cluster_df cluster partitions as returned by [make_clusters()]
#' @param file path to save the alignments. If `file = "example.nex"`, alignments will be saved to files `"example_part1.nex"`, `"example_part2.nex"`, etc.
#'
#' @returns No return value.
#'
#' @examples
#' # Load example phylogenetic data matrix
#' data("characters")
#'
#' # Create distance matrix
#' Dmatrix <- get_gower_dist(characters)
#'
#' # Find optimal partitioning scheme using PAM under k=3 partitions
#' cluster_df <- make_clusters(Dmatrix, k = 3)
#'
#' # Write to Nexus files
#' \dontrun{write_partitioned_alignments(characters, cluster_df, "example.nex")}

#' @export
write_partitioned_alignments <- function(x, cluster_df, file) {
  if (is.matrix(x) || is.data.frame(x) ||
      (is.list(x) && all(lengths(x) == length(x[[1L]])))) {
    data <- as.data.frame(as.matrix(x))
  }
  else if (length(x) == 1L && is.character(x) && endsWith(x, ".nex")) {
    data <- as.data.frame(ape::read.nexus.data(x))
  }
  else stop("'x' must be a data frame, matrix, or file path to a .nex file.", call. = FALSE)

  nk <- length(levels(cluster_df$cluster))
  file <- tools::file_path_sans_ext(file)

  for(ii in seq_len(nk)) {
    charset <- cluster_df$character_number[cluster_df$cluster == ii]
    aln <- lapply(as.data.frame(data), function(char) char[charset])

    fn <- paste0(file, "_part", ii, ".nex")
    ape::write.nexus.data(aln, file = fn, format = "standard", interleaved = FALSE)

  }
}
