#' Export character partitions to a Nexus file
#'
#' Creates and exports a Nexus file with a list of characters and their
#' respective partitions as inferred by [make_clusters()]. The contents can be copied and pasted directly into a Mr. Bayes
#' commands block for a partitioned clock Bayesian inference analysis.
#'
#' @param cluster_df A `cluster_df` object; the output of a call to
#' [make_clusters()].
#' @param file The path of the text file to be created containing the
#' partitioning information in Nexus format. If `NULL` (the default), no
#' file will be written and the output will be returned as a string. If
#' `""`, the text will be printed to the console. Passed directly to the
#' `file` argument of [cat()].
#'
#' @returns
#' The text as a string, returned invisibly if `file` is not
#' `NULL`. Use [cat()] on the resulting output to format it
#' correctly (i.e., to turn `"\n"` into line breaks).
#'
#' @seealso
#' `vignette("char-part")` for the use of this function as part
#' of an analysis pipeline.
#'
#' [make_clusters()]
#'
#' @examples
#' # Load example phylogenetic data matrix
#' data("characters")
#'
#' # Create distance matrix
#' Dmatrix <- get_gower_dist(characters)
#'
#' # Find optimal partitioning scheme using PAM under k=3
#' # partitions
#' cluster_df <- make_clusters(Dmatrix, k = 3)
#'
#' # Write to Nexus file and export to .txt file:
#' file <- tempfile(fileext = ".txt")
#'
#' # You would set, e.g.,
#' # file <- "path/to/file.txt"
#'
#' cluster_to_nexus(cluster_df, file = file)
#'

#' @export
cluster_to_nexus <- function(cluster_df, file = NULL) {

  cluster <- levels(cluster_df$cluster)
  vector <- vapply(levels(cluster_df$cluster), function(i) {
    paste(cluster_df$character_number[cluster_df$cluster == i], collapse = " ")
  }, character(1L))

  out <- paste0("#NEXUS\n\n",
                paste(paste0("charset morph_p", cluster, " = ", vector, ";"), collapse = "\n"),
                "\nset partition=all;\n")

  if (is.null(file)) {
    return(out)
  }

  cat(out, file = file)
  return(invisible(out))
}
