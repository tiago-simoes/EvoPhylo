#' Estimate and plot character partitions
#'
#' Determines cluster (partition) membership for phylogenetic morphological
#' characters from the supplied Gower distance matrix and requested number of
#' clusters using partitioning around medoids (PAM, or K-medoids). For further
#' and independently testing the quality of the chosen partitioning scheme,
#' users may also poduce graphic clustering (tSNEs), coloring data points
#' according to PAM clusters, to verify PAM clustering results.
#'
#' @details
#' `make_clusters()` calls [cluster::pam()] on the
#' supplied Gower distance matrix with the specified number of clusters to
#' determine cluster membership for each character. PAM is analogous to
#' K-means, but it has its clusters centered around medoids instead of centered
#' around centroids, which are less prone to the impact from outliers and
#' heterogeneous cluster sizes. PAM also has the advantage over k-means of
#' utilizing Gower distance matrices instead of Euclidean distance matrices
#' only.
#'
#' When `tsne = TRUE`, a Barnes-Hut t-distributed stochastic neighbor
#' embedding is used to compute a multi-dimensional embedding of the distance
#' matrix, coloring data points according to the PAM-defined clusters, as
#' estimated by the function `make_clusters`. This graphic clustering
#' allows users to independently test the quality of the chosen partitioning
#' scheme from PAM, and can help in visualizing the resulting clusters.
#' `Rtsne::Rtsne()` is used to do this. The resulting
#' dimensions will be included in the output; see Value below.
#'
#' `plot()` plots all morphological characters in a scatterplot with
#' points colored based on cluster membership. When `tsne = TRUE` in the
#' call to `make_clusters()`, the x- and y-axes will correspond to
#' requested tSNE dimensions. With more than 2 dimensions, several plots will
#' be produced, one for each pair of tSNE dimensions. These are displayed
#' together using [patchwork::plot_layout()].
#' When `tsne = FALSE`, the points will be arrange horizontally by cluster
#' membership and randomly placed vertically.
#'
#' @param dist_mat A Gower distance matrix, the output of a call to
#' [get_gower_dist()].
#' @param k The desired number of clusters (or character partitions), the
#' output from [get_sil_widths()].
#' @param tsne Whether to perform Barnes-Hut t-distributed stochastic neighbor
#' embedding (tSNE) to produce a multi-dimensional representation of the
#' distance matrix using `Rtsne::Rtsne()`. The number of
#' dimensions is controlled by the `tsne_dim` argument. See Details.
#' Default is `FALSE`.
#' @param tsne_dim When `tsne = TRUE`, the number of dimensions for the
#' tSNE multidimensional scaling plots. This is passed to the `dims`
#' argument of `Rtsne::Rtsne()`. Default is 2.
#' @param tsne_theta When `tsne = TRUE`, a parameter controlling the
#' speed/accuracy trade-off (increase for faster but less accurate results).
#' This is passed to the `theta` argument of
#' `Rtsne::Rtsne()`. Default is 0 for exact tSNE.
#' @param \dots For `make_clusters()`, other arguments passed to
#' `Rtsne::Rtsne()` when `tsne = TRUE`.
#'
#' For `plot()`, when plotting a `cluster_df` object, other arguments
#' passed to [ggrepel::geom_label_repel()]
#' to control display of the observation labels.
#'
#' @param x For `plot()`, a `cluster_df` object; the output of a call
#' to `make_clusters()`.
#' @param seed For `plot()`, the seed used to control the placement of the
#' labels and the jittering of the points. Jittering only occurs when
#' `tsne = FALSE` in the call to `make_clusters()`. Using a
#' non-`NA` seed ensure replicability across uses.
#' @param nrow For `plot()`, when `tsne = TRUE` in the call to
#' `make_clusters()` and `tsne_dim` is greater than 2, the number of
#' rows used to display the resulting 2-dimensional plots. Default is 1 for
#' side-by-side plots.
#'
#' @returns
#' A data frame, inheriting from class `"cluster_df"`, with a row
#' for each character with its number (`character_number`) and cluster
#' membership (`cluster`). When `tsne = TRUE`, additional columns
#' will be included, one for each requested tSNE dimension, labeled
#' `tSNE_Dim1`, `tSNE_Dim2`, etc., containing the values on the
#' dimensions computed using `Rtsne()`.
#'
#' The `pam` fit resulting from `cluster::pam` is returned in the
#' `"pam.fit"` attribute of the outut object.
#'
#' @note
#' When using `plot()` on a `cluster_df` object, warnings may
#' appear from `ggrepel` saying something along the lines of "unlabeled
#' data points (too many overlaps). Consider increasing max.overlaps". See [ggrepel::geom_label_repel()] for
#' details; the `max.overlaps` argument can be supplied to `plot()`
#' to increase the maximum number of element overlap in the plot.
#' Alternatively, users can increase the size of the plot when exporting it, as
#' it will increase the plot area and reduce the number of elements overlap.
#' This warning can generally be ignored, though.
#'
#' @seealso
#' `vignette("char-part")` for the use of this function as part
#' of an analysis pipeline.
#'
#' [get_gower_dist()], [get_sil_widths()],
#' [cluster_to_nexus()]
#'
#' [cluster::pam()], `Rtsne::Rtsne()`
#'
#' @examples
#' data("characters")
#'
#' # Reading example file as categorical data
#' Dmatrix <- get_gower_dist(characters)
#'
#' sil_widths <- get_sil_widths(Dmatrix, max.k = 7)
#'
#' sil_widths
#' # 3 clusters yields the highest silhouette width
#'
#' # Create clusters with PAM under k=3 partitions
#' cluster_df <- make_clusters(Dmatrix, k = 3)
#'
#' # Simple plot of clusters
#' plot(cluster_df, seed = 12345)
#'
#' # Create clusters with PAM under k=3 partitions and perform
#' # tSNE (3 dimensions; default is 2)
#' cluster_df_tsne <- make_clusters(Dmatrix, k = 3, tsne = TRUE,
#'                                  tsne_dim = 2)
#'
#' # Plot clusters, plots divided into 2 rows, and increasing
#' # overlap of text labels (default = 10)
#' plot(cluster_df_tsne, nrow = 2, max.overlaps = 20)
#'

#' @export
make_clusters <- function(dist_mat, k, tsne = FALSE, tsne_dim = 2,
                          tsne_theta = 0, ...) {
  pam.fit <- cluster::pam(dist_mat, diss = TRUE, k = k)

  df.Clusters <- data.frame(character_number = seq_along(pam.fit$clustering),
                            cluster = factor(pam.fit$clustering, levels = seq_len(k)))

  if (tsne) {
    rlang::check_installed("Rtsne")

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

#' @exportS3Method plot cluster_df
#' @rdname make_clusters
plot.cluster_df <- function(x, seed = NA, nrow = 1, ...) {

  if (sum(startsWith(names(x), "tSNE_Dim")) > 1) {

    names(x) <- gsub("tSNE_Dim", "tSNE Dimension ", names(x), fixed = TRUE)

    dim_combs <- combn(names(x)[startsWith(names(x), "tSNE Dimension")], 2L, simplify = FALSE)

    plots <- lapply(seq_along(dim_combs), function(i) {
      ggplot(x, aes(x = .data[[dim_combs[[i]][1L]]],
                    y = .data[[dim_combs[[i]][2L]]])) +
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
