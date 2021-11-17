#Implementing the procedure in CharPartitioning_Tetrapods.R

#Load data and compute Gower distance matrix
get_gower_dist <- function(file, numeric = FALSE) {
  if (length(file) == 1 && (is.matrix(file) || is.data.frame(file))) {
    Data_M <- as.data.frame(as.matrix(file))
  }
  else if (length(file) == 1 && is.character(file) && endsWith(file, ".nex")) {
    Data_M <- ape::read.nexus.data(file)
    Data_M <- as.data.frame(Data_M)
  }
  else stop("'file' must be a file path to a .nex file or a data frame.", call. = FALSE)

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
make_clusters <- function(Dmatrix, k, tsne = FALSE, plot = tsne, tsne_dim = 2, ...) {
  pam.fit <- cluster::pam(Dmatrix, diss = TRUE, k = k)

  df.Clusters <- data.frame(character_number = seq_along(pam.fit$clustering),
                            cluster = factor(pam.fit$clustering, levels = seq_len(k)))

  if (tsne) {
    tsne <- Rtsne::Rtsne(Dmatrix, is_distance = TRUE, theta = 0, dims = tsne_dim)
    for (i in seq_len(tsne_dim)) {
      df.Clusters[[paste0("tSNE_Dim", i)]] <- tsne$Y[,i]
    }
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

plot.cluster_df <- function(x, seed = 100, nrow = 1, ...) {
  names(x)[names(x) == "cluster"] <- "Cluster"

  if (sum(startsWith(names(x), "tSNE_Dim")) > 1) {

    names(x) <- gsub("tSNE_Dim", "tSNE Dimension ", names(x), fixed = TRUE)

    dim_combs <- combn(names(x)[startsWith(names(x), "tSNE Dimension")], 2, simplify = FALSE)

    plots <- lapply(seq_along(dim_combs), function(i) {
      ggplot(x, aes(x = .data[[dim_combs[[i]][1]]],
                    y = .data[[dim_combs[[i]][2]]])) +
        geom_point(aes(color = Cluster)) +
          ggrepel::geom_text_repel(aes(label = character_number, color = Cluster),
                                   show.legend = FALSE, ...) +
        theme_bw()
    })

    p <- Reduce("+", plots) + patchwork::plot_layout(nrow = nrow, guides = "collect")
  }
  else {
    pos <- position_jitter(seed = seed, width = .3)
    x$ypos <- 1
    p <- ggplot(x, aes(x = cluster, y = ypos)) +
      geom_point(aes(color = cluster), position = pos) +
      ggrepel::geom_text_repel(aes(label = character_number, color = cluster),
                      show.legend = FALSE, position = pos, ...) +
      theme_bw() + theme(axis.title.y = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks.y = element_blank(),
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
