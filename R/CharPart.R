#Implementing the procedure in CharPartitioning_Tetrapods.R

#Load data and compute Gower distance matrix
get_gower_dist <- function(file, numeric = FALSE) {
  if (is.matrix(file) || is.data.frame(file)) {
    Data_M <- as.data.frame(t(as.matrix(file)))
  }
  else if (is.character(file)) {
    Data_M <- as.data.frame(t(as.matrix(read.csv(file, colClasses = "character"))))
  }
  else stop("'file' must be a file path to a csv file or a data frame.", call. = FALSE)

  if (numeric) Data_M[] <- lapply(Data_M, as.numeric)

  Dmatrix <- StatMatch::gower.dist(Data_M, KR.corr = TRUE)

  return(Dmatrix)
}

#Compute silhouette widths for different numbers of clusters; optionally plot
get_sil_widths <- function(Dmatrix, max.k = 10, plot = TRUE) {

  sil_width <- vapply(2:max.k, function(i) {
    cluster::pam(Dmatrix, diss = TRUE, k = i)$silinfo$avg.width
  }, numeric(1L))

  sil_width_df <- data.frame(k = 2:max.k, sil_width = sil_width)
  class(sil_width_df) <- c("sil_width_df", "data.frame")

  if (plot) {
    print(
      plot.sil_width_df(sil_width_df)
    )
    invisible(sil_width_df)
  }
  else {
    return(sil_width_df)
  }
}

#Plot silhouette widths
plot.sil_width_df <- function(x, ...) {
  ggplot(x, aes(x = k, y = sil_width)) +
    geom_line() +
    labs(x = "Number of clusters", y = "Silhouette Width") +
    theme_bw()
}

#Perform clustering,
make_clusters <- function(Dmatrix, k, tsne = FALSE, plot = tsne, ...) {
  pam.fit <- cluster::pam(Dmatrix, diss = TRUE, k = k)

  df.Clusters <- data.frame(character_number = seq_along(pam.fit$clustering),
                            cluster = factor(pam.fit$clustering, levels = seq_len(k)))

  if (tsne) {
    tsne <- Rtsne::Rtsne(Dmatrix, is_distance = TRUE, theta = 0, dims = 2)
    df.Clusters$tSNE_Dim1 <- tsne$Y[,1]
    df.Clusters$tSNE_Dim2 <- tsne$Y[,2]
  }

  class(df.Clusters) <- c("cluster_df", "data.frame")

  if (plot) {
    print(
        plot.cluster_df(df.Clusters, ...)
    )
  }

  attr(df.Clusters, "pam.fit") <- pam.fit

  return(df.Clusters)
}

plot.cluster_df <- function(x, seed = 100, ...) {
  if (ncol(x) >= 4) {
    p <- ggplot(x, aes(x = tSNE_Dim1, y = tSNE_Dim2)) +
      geom_point(aes(color = cluster)) + suppressWarnings(
        ggrepel::geom_text_repel(aes(label = character_number, color = cluster),
                      show.legend = FALSE, ...)) +
      theme_bw()
  }
  else {
    pos <- position_jitter(seed = seed, width = .3)
    x$ypos <- 1
    p <- ggplot(x, aes(x = cluster, y = ypos)) +
      geom_point(aes(color = cluster), position = pos) +
      ggrepel::geom_text_repel(aes(label = character_number, color = cluster),
                      show.legend = FALSE, position = pos, ...) +
      theme_bw() + theme(axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         panel.grid = element_blank())
  }
  p
}

#Write cluster output to a Nexus file and returning the cluster
cluster_to_nexus <- function(clusters, file = "") {

  s <- data.frame(cluster = levels(clusters$cluster),
                  vector = vapply(levels(clusters$cluster), function(i) {
                    paste(clusters$character_number[clusters$cluster==i], collapse = " ")
                  }, character(1L)))

  out <- paste0("#NEXUS\n\n",
                paste(paste0("charset morph_p", s$cluster," = ", s$vector, ";"), collapse = "\n"),
                "\nset partition=all;\n")

  cat(out, file = file)

  return(invisible(out))
}
