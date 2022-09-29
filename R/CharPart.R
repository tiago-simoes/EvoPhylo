#All functions documented with examples

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

#Compute silhouette widths for different numbers of clusters; optionally plot
get_sil_widths <- function(dist_mat, max.k = 10) {

  sil_width <- vapply(2:max.k, function(i) {
    cluster::pam(dist_mat, diss = TRUE, k = i)$silinfo$avg.width
  }, numeric(1L))

  sil_width_df <- data.frame(k = 2:max.k, sil_width = sil_width)
  class(sil_width_df) <- c("sil_width_df", "data.frame")

  return(sil_width_df)
}

#Plot silhouette widths
plot.sil_width_df <- function(x, ...) {
  ggplot(data = x) +
    geom_line(aes(x = .data$k, y = .data$sil_width), ...) +
    labs(x = "Number of clusters", y = "Silhouette Width") +
    theme_bw()
}

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


#' Write character partitions as separate Nexus files (for use in BEAUti)
#'
#' @param characters character data matrix as matrix or data frame (with taxa as columns and characters as rows)
#' @param cluster_df cluster partitions as outputted by \code{make.clusters}
#' @param file path to save the alignments. If \code{file = "example.nex"}, alignments will be saved to files \code{"example_part1.nex"}, \code{"example_part2.nex"}, etc.
#'
#' @return no return value
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
write_partitioned_alignments <- function(characters, cluster_df, file) {
  nk <- length(levels(cluster_df$cluster))
  file <- tools::file_path_sans_ext(file)

  for(ii in 1:nk) {
    charset <- cluster_df$character_number[cluster_df$cluster==ii]
    aln <- lapply(as.data.frame(characters), function(char) char[charset])
    
    fn <- paste0(file, "_part", ii, ".nex")
    ape::write.nexus.data(aln, file = fn, format = "standard", interleaved = FALSE)

  }
}
