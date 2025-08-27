#' Calculate silhouette widths index for various numbers of partitions
#'
#' Computes silhouette widths index for several possible numbers of
#' clusters(partitions) `k`, which determines how well an object falls
#' within their cluster compared to other clusters. The best number of clusters
#' `k` is the one with the highest silhouette width.
#'
#' @details
#' `get_sil_widths()` calls [cluster::pam()] on the
#' supplied Gower distance matrix with each number of clusters (partitions) up
#' to `max.k` and stores the average silhouette widths across the
#' clustered characters. When `plot = TRUE`, a plot of the sillhouette
#' widths against the number of clusters is produced, though this can also be
#' produced separately on the resulting data frame using
#' `plot.sil_width_df()`. The number of clusters with the greatest
#' silhouette width should be selected for use in the final clustering
#' specification.
#'
#' @param dist_mat A Gower distance matrix, the output of a call to [get_gower_dist()].
#' @param max.k The maximum number of clusters(partitions) to search across.
#' @param x A `sil_width_df` object; the output of a call to `get_sil_widths()`.
#' @param \dots Further arguments passed to [ggplot2::geom_path()] to control the
#' appearance of the plot.
#'
#' @returns
#' For `get_sil_widths()`, it produces a data frame, inheriting
#' from class `"sil_width_df"`, with two columns: `k` is the number
#' of clusters, and `sil_width` is the silhouette widths for each number
#' of clusters. If `plot = TRUE`, the output is returned invisibly.
#'
#' For `plot()` on a `get_sil_widths()` object, it produces a
#' `ggplot` object that can be manipulated using \pkg{ggplot2} syntax
#' (e.g., to change the `theme` or labels).
#'
#' @seealso
#' `vignette("char-part")` for the use of this function as part
#' of an analysis pipeline.
#'
#' [get_gower_dist()], [cluster::pam()]
#'
#' @examples
#' data("characters")
#'
#' #Reading example file as categorical data
#' Dmatrix <- get_gower_dist(characters)
#'
#' #Get silhouette widths for k=7
#' sw <- get_sil_widths(Dmatrix, max.k = 7)
#'
#' sw
#'
#' plot(sw, color = "red", size =2)

#' @export
get_sil_widths <- function(dist_mat, max.k = 10) {

  sil_width <- vapply(seq_len(max.k)[-1L], function(i) {
    cluster::pam(dist_mat, diss = TRUE, k = i)$silinfo$avg.width
  }, numeric(1L))

  sil_width_df <- data.frame(k = seq_len(max.k)[-1L], sil_width = sil_width)
  class(sil_width_df) <- c("sil_width_df", "data.frame")

  return(sil_width_df)
}

#' @exportS3Method plot sil_width_df
#' @rdname get_sil_widths
plot.sil_width_df <- function(x, ...) {
  ggplot(data = x) +
    geom_line(aes(x = .data$k, y = .data$sil_width), ...) +
    labs(x = "Number of clusters", y = "Silhouette Width") +
    theme_bw()
}
