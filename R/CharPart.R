#All functions documented with examples

# R/globals.R or top of any script file
utils::globalVariables(c("value", "%>%"))

#Load data and compute Gower distance matrix
get_gower_dist <- function(x, numeric = FALSE) {
  if (is.matrix(x) || is.data.frame(x) ||
      (is.list(x) && all(lengths(x) == length(x[[1]])))) {
    data <- as.data.frame(as.matrix(x))
  }
  else if (length(x) == 1 && is.character(x) && endsWith(x, ".nex")) {
    data <- as.data.frame(ape::read.nexus.data(x))
  }
  else stop("'x' must be a data frame, matrix, or file path to a .nex file.", call. = FALSE)

  if (numeric) data[] <- safe_as.numeric(data)

  return(gowerdist(data))
}


#--------------------------------------------------------------------------------------
#Compute silhouette widths for different numbers of clusters; optionally plot
get_sil_widths <- function(dist_mat, max.k = 10) {

  sil_width <- vapply(2:max.k, function(i) {
    cluster::pam(dist_mat, diss = TRUE, k = i)$silinfo$avg.width
  }, numeric(1L))

  sil_width_df <- data.frame(k = 2:max.k, sil_width = sil_width)
  class(sil_width_df) <- c("sil_width_df", "data.frame")

  return(sil_width_df)
}

#--------------------------------------------------------------------------------------
#Plot silhouette widths
plot.sil_width_df <- function(x, ...) {
  ggplot(data = x) +
    geom_line(aes(x = .data$k, y = .data$sil_width), ...) +
    labs(x = "Number of clusters", y = "Silhouette Width") +
    theme_bw()
}


#--------------------------------------------------------------------------------------
#Perform clustering,
make_clusters <- function(dist_mat, k, tsne = FALSE, tsne_dim = 2,
                          tsne_theta = 0, ...) {
  pam.fit <- cluster::pam(dist_mat, diss = TRUE, k = k)

  df.Clusters <- data.frame(character_number = seq_along(pam.fit$clustering),
                            cluster = factor(pam.fit$clustering, levels = seq_len(k)))

  if (tsne) {
    if (!requireNamespace("Rtsne", quietly = TRUE)) {
      stop("The {Rtsne} package must be installed to use tSNE.\nInstall it using install.packages(\"Rtsne\").", call. = FALSE)
    }
    args <- list(...)
    args[["X"]] <- dist_mat
    args[["dims"]] <- tsne_dim
    args[["theta"]] <- tsne_theta
    args[["is_distance"]] <- TRUE

    tsne <- do.call(Rtsne::Rtsne, args)

    for (i in seq_len(tsne_dim)) {
      df.Clusters[[paste0("tSNE_Dim", i)]] <- tsne$Y[,i]
    }
  }

  class(df.Clusters) <- c("cluster_df", "data.frame")

  attr(df.Clusters, "pam.fit") <- pam.fit

  return(df.Clusters)
}


#--------------------------------------------------------------------------------------
plot.cluster_df <- function(x, seed = NA, nrow = 1, ...) {

  if (sum(startsWith(names(x), "tSNE_Dim")) > 1) {

    names(x) <- gsub("tSNE_Dim", "tSNE Dimension ", names(x), fixed = TRUE)

    dim_combs <- combn(names(x)[startsWith(names(x), "tSNE Dimension")], 2, simplify = FALSE)

    plots <- lapply(seq_along(dim_combs), function(i) {
      ggplot(x, aes(x = .data[[dim_combs[[i]][1]]],
                    y = .data[[dim_combs[[i]][2]]])) +
        geom_point(aes(color = .data$cluster)) +
          ggrepel::geom_text_repel(aes(label = .data$character_number, color = .data$cluster),
                                   show.legend = FALSE, seed = seed, ...) +
        theme_bw() + labs(color = "Cluster")
    })

    p <- patchwork::wrap_plots(plots, nrow = nrow, guides = "collect")
  }
  else {
    pos <- position_jitter(seed = seed, width = .3)

    p <- ggplot(x, aes(x = .data$cluster, y = 1)) +
      geom_point(aes(color = .data$cluster), position = pos) +
      ggrepel::geom_text_repel(aes(label = .data$character_number, color = .data$cluster),
                      show.legend = FALSE, position = pos, seed = seed, ...) +
      theme_bw() + theme(axis.title.y = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks.y = element_blank(),
                         panel.grid = element_blank()) +
      labs(color = "Cluster")
  }
  p
}

#--------------------------------------------------------------------------------------
#Write cluster output to a Nexus file
cluster_to_nexus <- function(cluster_df, file = NULL) {

  cluster <- levels(cluster_df$cluster)
  vector <- vapply(levels(cluster_df$cluster), function(i) {
    paste(cluster_df$character_number[cluster_df$cluster==i], collapse = " ")
  }, character(1L))

  out <- paste0("#NEXUS\n\n",
                paste(paste0("charset morph_p", cluster," = ", vector, ";"), collapse = "\n"),
                "\nset partition=all;\n")

  if (is.null(file)) {
    return(out)
  }
  else {
    cat(out, file = file)
    return(invisible(out))
  }
}

#-------------------------------------------------------------------------------------------------


#' Write character partitions as separate Nexus files (for use in BEAUti)
#'
#' @param x character data matrix as Nexus file (.nex) or data frame (with taxa as rows and characters as columns) read directly from local directory
#' @param cluster_df cluster partitions as outputted by \code{make.clusters}
#' @param file path to save the alignments. If \code{file = "example.nex"}, alignments will be saved to files \code{"example_part1.nex"}, \code{"example_part2.nex"}, etc.
#'
#' @seealso \code{\link{write_partitioned_alignments2}} expanded function for exporting both  morphological and molecular partitions.
#' @return No return value. This function is called for its side effect of writing alignment files to disk.
#' @export
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
#'
write_partitioned_alignments <- function(x, cluster_df, file) {
    if (is.matrix(x) || is.data.frame(x) ||
      (is.list(x) && all(lengths(x) == length(x[[1]])))) {
    data <- as.data.frame(as.matrix(x))
  }
  else if (length(x) == 1 && is.character(x) && endsWith(x, ".nex")) {
    data <- as.data.frame(ape::read.nexus.data(x))
  }
  else stop("'x' must be a data frame, matrix, or file path to a .nex file.", call. = FALSE)

  nk <- length(levels(cluster_df$cluster))
  file <- tools::file_path_sans_ext(file)

  for(ii in 1:nk) {
    charset <- cluster_df$character_number[cluster_df$cluster==ii]
    aln <- lapply(as.data.frame(data), function(char) char[charset])

    fn <- paste0(file, "_part", ii, ".nex")
    ape::write.nexus.data(aln, file = fn, format = "standard", interleaved = FALSE)

  }
}


#-------------------------------------------------------------------------------------------------
#' Write alignment partitions as separate alignment files for various data types
#'
#' @param x concatenated alignment file in Nexus or Phyllip format read directly from local directory
#' @param cluster_df cluster partitions as outputted by \code{make.clusters}
#' @param partition_file name of text file with user provided partitions, with names and start&end positions
#' @param in_format Format of the input alignment file. One of "phylip", "interleaved", "sequential", "clustal", "fasta", or "nexus", or any unambiguous abbreviation. Passed to \code{phangorn::read.phyDat}.
#' @param in_type Type of input sequences. One of "DNA", "AA", "CODON" or "USER". Passed to \code{phangorn::read.phyDat}.
#' @param out_file Path to save the alignments. If \code{out_file = "example.nex"}, files will be saved as \code{"example_part1.nex"}, \code{"example_part2.nex"}, etc.
#' @param out_type Output format type. One of "dna" (default), "protein", "standard", or "continuous".
#'
#' @seealso \code{\link{write_partitioned_alignments}} for the older version supporting morphological data only.
#' @return No return value. This function is called for its side effect of writing alignment files to disk.
#' @export
#' @importFrom magrittr %>%
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
#' # Write morphological partitions into multiple Nexus files
#' \dontrun{write_partitioned_alignments2(x = "characters.nex",
#'                              cluster_df = cluster_df,
#'                              out_file = "test", out_type = "standard")}
#'
#' # Write to molecular partitions into multiple Phyllip files
#' \dontrun{write_partitioned_alignments2(x = "alignments.phy",
#'                              partition_file = "sorted_partitions_50genes.txt",
#'                              in_format = "phylip", in_type = "dna",
#'                              out_file = "test", out_type = "dna")}
#'
#' @md
### Function

write_partitioned_alignments2 <- function(x, cluster_df, partition_file,
                                          in_format = NULL, in_type = NULL,
                                          out_file, out_type = "standard") {

#Get data matrix as data frame
  if (is.matrix(x) || is.data.frame(x) ||
      (is.list(x) && all(lengths(x) == length(x[[1]])))) {
    data <- as.data.frame(as.matrix(x))
  }
  else if (length(x) == 1 && is.character(x) && endsWith(x, ".nex")) {
    data <- as.data.frame(ape::read.nexus.data(x))
  }
  else if (length(x) == 1 && is.character(x) && (endsWith(x, ".phy") || endsWith(x, ".fa"))) {
    data <- as.data.frame(phangorn::read.phyDat(x, format = in_format, type = in_type))
  }
  else stop("'x' must be a data frame, matrix, or file path to a Nexus(.nex), Phylip (.phy), or Fasta (.fa) file", call. = FALSE)

#Get partition names and start&end positions
#From a user provided partition file

 if(length(partition_file) == 1 && is.character(partition_file)) {

    partitions <- read.table(partition_file, sep = ' ')
    loc_names <- as.character(unlist(tibble::enframe(partitions[,(which(partitions[1,] == '=') - 1)], name = NULL), use.names = F))
    partitions <- tibble::enframe(partitions[,ncol(partitions)], name = NULL)
    partitions <- partitions %>% tidyr::separate(value, into = c('Start', 'End'), sep = '-') %>%
                      dplyr::mutate_if(is.character, as.numeric)

#Write data for each partition as a separate file
  nk <- length(loc_names)
  out_file <- tools::file_path_sans_ext(out_file)

   for(j in 1:nk) {
    charset <- c()
    charset <- c(charset, partitions$Start[j]:partitions$End[j])
    aln <- lapply(as.data.frame(data), function(char) char[charset])

    fn <- paste0(out_file, "_part", j, "_locus", loc_names[j],".nex")
    ape::write.nexus.data(aln, file = fn, format = out_type, interleaved = FALSE)
    }}

#From a cluster file generated with 'make_clusters'
 else if (is.data.frame(cluster_df))  {

  nk <- length(levels(cluster_df$cluster))
  out_file <- tools::file_path_sans_ext(out_file)

  for(j in 1:nk) {
    charset <- cluster_df$character_number[cluster_df$cluster==j]
    aln <- lapply(as.data.frame(data), function(char) char[charset])

    fn <- paste0(out_file, "_part", j, ".nex")
    ape::write.nexus.data(aln, file = fn, format = out_type, interleaved = FALSE)
  }}

  else stop("either 'cluster_df' must be a data frame or 'partition_file' must be a partition file", call. = FALSE)

}

