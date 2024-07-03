library(igraph)
dat_fun <- function(){ worcs::load_data(to_envir = FALSE, verbose = FALSE)$embeddings}

do_umap <- function(mat_vect){
# Reduce dimensionality
# we will want a larger n_neighbors value – small values will focus more on very local structure and are more prone to producing fine grained cluster structure that may be more a result of patterns of noise in the data than actual clusters.
# Second set min_dist to a very low value. Since we actually want to pack points together densely (density is what we want after all) a low value will help, as well as making cleaner separations between clusters. In this case we will simply set min_dist to be 0.

# clusterable_embedding = umap.UMAP(
#   n_neighbors=30,
#   min_dist=0.0,
#   n_components=2,
#   random_state=42,
# )
  library(umap)
  library(reticulate)
  use_virtualenv("tmsrcoda")
  pypkgs <- py_list_packages()
  if(!"umap-learn" %in% pypkgs$package) reticulate::py_install("umap-learn")
  config_umap <- umap.defaults
  config_umap$n_neighbors <- 5
  config_umap$n_components <- 10
  config_umap$min_dist <- 0.0001
  config_umap$metric <- "cosine"
  mat_vect[1:4,1:4]
  umap(d = mat_vect, config = config_umap, method = "umap-learn")
}

do_hdbscan <- function(mat, algorithm = "best",
                       allow_single_cluster = FALSE,
                       alpha = 1,
                       prediction_data = TRUE,
                       approx_min_span_tree = TRUE,
                       gen_min_span_tree = FALSE,
                       leaf_size = 40L,
                       core_dist_n_jobs = 19L,
                       metric = "euclidean",
                       min_cluster_size = 10L,
                       min_samples = 1L,
                       cluster_selection_epsilon =  .1,
                       cluster_selection_method = "leaf"){
  reticulate::use_virtualenv("tmsrcoda")
  pypkgs <- reticulate::py_list_packages()
  if(!"hdbscan" %in% pypkgs$package) reticulate::py_install("hdbscan")
  hdbscan <- reticulate::import('hdbscan', delay_load = TRUE)

  clusterer <- hdbscan$HDBSCAN(algorithm = algorithm,
                               allow_single_cluster = allow_single_cluster,
                               alpha = alpha,
                               prediction_data = prediction_data,
                               approx_min_span_tree = approx_min_span_tree,
                               gen_min_span_tree = gen_min_span_tree,
                               leaf_size = leaf_size,
                               core_dist_n_jobs = core_dist_n_jobs,
                               metric = metric,
                               min_cluster_size = min_cluster_size,
                               min_samples = min_samples,
                               cluster_selection_epsilon =  cluster_selection_epsilon,
                               cluster_selection_method = cluster_selection_method)

  reticulate::py_main_thread_func(clusterer$fit(mat))

  return(list(
    labels = clusterer$labels_,
    probabilities = clusterer$probabilities_,
    cluster_persistance = clusterer$cluster_persistence_,
    exemplars = clusterer$exemplars_,
    outlier_scores = clusterer$outlier_scores_))

}

remove_outliers <- function(mat_umap){

    mns_umap <- colMeans(mat_umap)
    cov_umap <- cov(mat_umap)
    dist_mah <- mahalanobis(x = mat_umap , center = mns_umap, cov = cov_umap)
    plot(density(dist_mah))

    # Cutoff value for ditances from Chi-Sqaure Dist.
    # with p = 0.95 df = 2 which in ncol(air)
    cutoff <- qchisq(p = 0.95, df = ncol(mat_umap))

    colr <- rep("black", nrow(mat_umap))
    outliers <- dist_mah > cutoff
    colr[outliers] <- "red"

    p <- ggplot(data.frame(x = mat_umap[,1], y = mat_umap[,2], col = colr), aes(x = x, y = y, color = col)) +
      geom_point(alpha = .05) +
      theme_bw()+
      labs(x = "UMAP dimension 1", y = "UMAP dimension 2")+
      scale_color_manual(values = c("black" = "black", "red"=  "red"))+
      theme(legend.position = "none")
    ggsave("umap_outliers.png", p, device = "png", width = 210, height = 200, units = "mm")
    ggsave("umap_outliers.svg", p, device = "svg")
    plot3d(x = mat_umap[,1],
           y = mat_umap[,2],
           z = mat_umap[,3],
           type="s",
           size=1,
           col = colr)
    # Check the worst words
    # tail(df$word[unique_words][outliers][order(dist_mah[outliers], decreasing = TRUE)])
    # Remove outliers
    mat_umap <- mat_umap[!outliers, ]
    plot3d(x = mat_umap[,1],
           y = mat_umap[,2],
           z = mat_umap[,3],
           type="s",
           size=1)

}


write_exemplar_metadata <- function(df, res_umap, res_hdbscan){
  all_clusters <- sort(unique(res_hdbscan$labels))
  # exemplars <- lapply(all_clusters, function(thisclust){
  #   #thisclust = 106
  #   print(thisclust)
  #   exmpls <- res_hdbscan$exemplars[[thisclust]]
  #   df$word[apply(exmpls, 1, FUN = function(x){which(rowSums(res_umap$layout == rep(x, each = nrow(res_umap$layout)))==length(x))})]
  # })
  words_by_clust <- split(df$word, factor(res_hdbscan$labels))
  next_words <- split(df$word[-1], factor(res_hdbscan$labels[-length(res_hdbscan$labels)]))
  tab_words <- lapply(words_by_clust, function(w){
    #w = words_by_clust[[371]]
    sort(table(textstem::lemmatize_strings(w)), decreasing = TRUE)
    })
  tab_next_words <- lapply(next_words, function(w){
    #w = words_by_clust[[371]]
    out <- sort(table(textstem::lemmatize_strings(w)), decreasing = TRUE)
    out[1:3]/sum(out)
  })
  labs <- lapply(tab_words, function(x){
    sample(names(x)[x == max(x)], 1)
  })
  out <- vector("list", length = length(all_clusters))
  for(i in seq_along(out)){
    out[[i]] <- list(
      label = labs[[i]],
      words = paste0(paste(names(tab_words[[i]]), tab_words[[i]]), collapse = ", "),
      next_words = paste0(paste(names(tab_next_words[[i]]), formatC(tab_next_words[[i]], digits = 2, format = "f")), collapse = ", "),
      include = TRUE
    )
  }
  names(out) <- all_clusters
  out[[1]]$label <- "noise"
  out[[1]]$include <- "FALSE"
  yaml::write_yaml(out, "exemplars.yaml")
  return("exemplars.yaml")
}


## Construct frequencies
pretty_words <- function(x){
  out <- x
  out[tolower(out) == "dysregulation"] <- "Emotion regulation"
  out <- stringr::str_to_sentence(gsub("[_-]", " ", out))
  out[out == "Adhd cd"] <- "ADHD/CD"
  out[out == "Ses"] <- "SES"
  out[out == "Ptsd"] <- "PTSD"
  out[out == "Schizo"] <- "Schizophrenia"
  out <- gsub("^Par\\.", "Parenting ", out)
  out
}

plot_cluster_freq <- function(exmplrs, res_hdbscan, prune = .99){
  library(ggplot2)
  library(svglite)

  df_plot <- do.call(rbind, lapply(exmplrs[-1], function(i){
    data.frame(lab = i$label,
               freq = sum(as.integer(gsub("^.*\\s(\\d{1,})$", "\\1", strsplit(i$words, ", ")[[1]]))),
               include = i$include
    )
  }))
  df_plot <- df_plot[df_plot$include, ]
  df_plot <- tapply(df_plot$freq, factor(df_plot$lab), sum)
  df_plot <- as.data.frame.table(df_plot)

  fit <- MASS::fitdistr(df_plot$Freq, "negative binomial")
  # thres <- qnbinom(p, size=fit$estimate["size"], mu=fit$estimate["mu"])
  # return(cooc > thres)

  p = ggplot(df_plot, aes(x = Freq)) +
    geom_density() +
    labs(x = "Number of words per cluster", y = "Density") +
    theme_bw() +
    scale_x_continuous(limits = range(df_plot$Freq), expand = c(0,0))

  p <- p + geom_function(data = data.frame(x = as.integer(sort(unique(df_plot$Freq)))), aes(x = x),
                    fun = function(x){ dnbinom(x = round(x), size=fit$estimate["size"], mu=fit$estimate["mu"])}, colour = "red") #+
    #geom_vline(xintercept = qnbinom(p = prune, size=fit$estimate["size"], mu=fit$estimate["mu"]))
  ggsave("plot_dist_clust.png", p, device = "png", width = 210, height = 200, units = "mm")
  ggsave("plot_dist_clust.svg", p, device = "svg", width = 210, height = 200, units = "mm")
  names(df_plot) <- c("Construct", "Frequency")
  write.csv(df_plot, "cluster_freq.csv", row.names = FALSE)

  df_plot$Construct <- pretty_words(as.character(df_plot$Construct))
  df_plot <- df_plot[order(df_plot$Frequency, decreasing = TRUE), ]
  df_plot$Construct <- ordered(df_plot$Construct, levels = df_plot$Construct[order(df_plot$Frequency)])

  #cat_cols <- c(Outcome = "gray50", Indicator = "tomato", Cause = "gold", Protective = "forestgreen")
  #df_plot$cat <- ordered(df_plot$cat, levels = c("Outcome", "Indicator", "Cause", "Protective"))
  #range(df_plot$Frequency)
  p <- ggplot(df_plot, aes(y = Construct, x = Frequency)) +
    geom_segment(aes(x = 0, xend = Frequency,
                     y = Construct, yend = Construct
                     #, linetype = faded
    ), colour = "grey50"
    ) +
    geom_vline(xintercept = 0, colour = "grey50", linetype = 1) + xlab("Construct frequency") +
    geom_point(data = df_plot, colour = "black", fill = "black", shape = 21, size = 1.5) +
    # scale_colour_manual(values = c(Outcome = "gray50", Indicator = "tomato", Cause = "gold", Protective = "forestgreen"), guide = NULL)+
    # scale_fill_manual(values = c(Outcome = "gray50", Indicator = "tomato", Cause = "gold", Protective = "forestgreen")) +
    scale_x_sqrt() +
    #scale_linetype_manual(values = c("TRUE" = 2, "FALSE" = 1), guide = NULL) +
    theme_bw() + theme(panel.grid.major.x = element_blank(),
                       panel.grid.minor.x = element_blank(), axis.title.y = element_blank(),
                       legend.position = c(.70,.125),
                       legend.title = element_blank(),
                       axis.text.y = element_text(hjust=0, vjust = 0, size = 6))

  ggsave("plot_freq.png", p, device = "png", width = 210, height = 400, units = "mm")
  ggsave("plot_freq.svg", p, device = "svg", width = 210, height = 400, units = "mm")
  return("plot_freq")
}
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }

  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

add.vertex.shape("circle2", clip=igraph.shape.noclip,
                 plot=mycircle, parameters=list(vertex.frame.color=1,
                                                vertex.frame.width=1))



select_cooc <- function(cooc, p = .95){
  fit <- MASS::fitdistr(cooc, "negative binomial")
  thres <- qnbinom(p, size=fit$estimate["size"], mu=fit$estimate["mu"])
  out <- cooc > thres
  attr(out, "thres") <- thres
  return(out)
}


library(igraph)
library(ggplot2)
library(ggrepel)
# df <- sotu::sotu_meta
# set.seed(1)
# dat <- data.frame(construct = factor(df$president), doc = sample.int(10, size = length(df$party), replace = TRUE))
#
# select_cooc <- function(cooc, r = .95){
#
#   fit <- MASS::fitdistr(cooc, "negative binomial")
#   thres <- qnbinom(q, size=fit$estimate["size"], mu=fit$estimate["mu"])
#   return(cooc > thres)
# }

plot_graph <- function(df, exmplrs, res_hdbscan, prune = .99){
  include <- as.logical(sapply(exmplrs, `[[`, "include")[as.character(res_hdbscan$labels)])
  labels <- sapply(exmplrs, `[[`, "label")[as.character(res_hdbscan$labels)]
  dat <- data.frame(doc = as.integer(factor(df$doc_id)),
                    include,
                    construct = labels)
  # Drop terms marked for exclusion
  dat <- dat[dat$include, -which(names(dat)=="include")]
  # Replace constructs with labels because some are merged
  dat$construct <- factor(dat$construct)

  duplicates = FALSE
  # Dat must have factor called construct and integer ID called doc

  # Words only count once per document
  if(!duplicates) dat <- dat[!duplicated(dat), ]
  # Make table
  cluster_tab <- table(dat$construct)
  cluster_freq <- data.frame(construct = names(cluster_tab), frequency = as.integer(cluster_tab))

  # Make coocurrence matrix
  V <- crossprod(table(dat[, c("doc", "construct")]))
  diag(V) <- 0

  df_plot <- as.data.frame.table(V)
  names(df_plot) <- c("term1", "term2", "cooc")
  df_plot <- df_plot[as.vector(lower.tri(V)), ]
  write.csv(df_plot, "complete_coocurrences.csv", row.names = FALSE)
  node_centrality <- sort(sapply(levels(df_plot$term1), function(l){ sum(df_plot$cooc[df_plot$term1 == l | df_plot$term2 == l]) })) #sort(igraph::degree(g_centr, mode="all"))
  df_nodecentral <- data.frame(construct = names(node_centrality),
                               frequency = cluster_freq$frequency[match(names(node_centrality), cluster_freq$construct)],
                               centrality_degree = node_centrality
  )

  write.csv(df_nodecentral, "complete_node_centrality.csv", row.names = FALSE)

  # Now prune
  if(!is.null(prune)){
    select_links <- suppressWarnings(select_cooc(df_plot$cooc, p = prune))
    writeLines(paste("Cooccurrence graph threshold: ", attr(select_links, "thres")), "coocurrence_threshold.txt")
    df_plot <- df_plot[select_links, ]
  }
  write.csv(df_plot, "coocurrences.csv", row.names = FALSE)
  # Calculate centrality
  node_centrality <- sort(sapply(levels(df_plot$term1), function(l){ sum(df_plot$cooc[df_plot$term1 == l | df_plot$term2 == l]) })) #sort(igraph::degree(g_centr, mode="all"))
  df_nodecentral <- data.frame(construct = names(node_centrality),
                               frequency = cluster_freq$frequency[match(names(node_centrality), cluster_freq$construct)],
                               centrality_degree = node_centrality
                               )

  write.csv(df_nodecentral, "node_centrality.csv", row.names = FALSE)

  # Plot as list ------------------------------------------------------------
  df_tab <- df_plot
  df_tab$term1 <- as.character(df_tab$term1)
  df_tab$term2 <- as.character(df_tab$term2)
  df_tab[c("term1", "term2")] <- do.call(rbind, apply(df_tab[c("term1", "term2")], 1, function(r){r[order(nchar(r), decreasing = T)]}, simplify = F))
  maxchar <- sapply(df_tab[c("term1", "term2")], function(i){max(nchar(i))})+1L
  df_tab$joint <- sapply(1:nrow(df_tab), function(i){
    paste0(sprintf(paste0("%-", maxchar[1], "s"), df_tab$term1[i]), "<->",
           sprintf(paste0("%", maxchar[2], "s"), df_tab$term2[i]))
  })
  df_tab$joint <- ordered(df_tab$joint, levels = df_tab$joint[order(df_tab$cooc, decreasing = F)])
  df_tab <- df_tab[order(df_tab$cooc, decreasing = T), ]
  df_tab <- df_tab[df_tab$cooc > attr(select_links, "thres"), ]

  p <- ggplot(df_tab, aes(y = joint, x = cooc)) +
    geom_segment(aes(x = 0, xend = cooc,
                     y = joint, yend = joint
                     #, linetype = faded
    ), colour = "grey50"
    ) +
    geom_vline(xintercept = 0, colour = "grey50", linetype = 1) + xlab("Co-occurrence frequency") +
    geom_point(data = df_tab, colour = "black", fill = "black", shape = 21, size = 1.5) +
    # scale_colour_manual(values = c(Outcome = "gray50", Indicator = "tomato", Cause = "gold", Protective = "forestgreen"), guide = NULL)+
    # scale_fill_manual(values = c(Outcome = "gray50", Indicator = "tomato", Cause = "gold", Protective = "forestgreen")) +
    scale_x_sqrt(expand = c(0,0)) +
    #scale_linetype_manual(values = c("TRUE" = 2, "FALSE" = 1), guide = NULL) +
    theme_bw() + theme(panel.grid.major.x = element_blank(),
                       panel.grid.minor.x = element_blank(), axis.title.y = element_blank(),
                       legend.position = c(.70,.125),
                       legend.title = element_blank(),
                       axis.text.y = element_text(hjust=0, vjust = 0, size = 6, family = "mono"))

  ggsave("plot_cooc_freq.png", p, device = "png", width = 210, height = 400, units = "mm", dpi = 150)
  ggsave("plot_cooc_freq.svg", p, device = "svg", width = 210, height = 400, units = "mm")

  # Without cooperation
  df_tab <- df_tab[!(df_tab$term1 == "cooperation"| df_tab$term2 == "cooperation"), ]
  df_nocoop <- do.call(rbind, lapply(1:nrow(df_tab), function(i){
    constructs <- df_tab[i, c(1:2)]
    data.frame(constructs[order(node_centrality[unlist(constructs)], decreasing = TRUE)], cooc = df_tab$cooc[i])
  }))
  names(df_nocoop) <- c("Central construct", "Peripheral construct", "Co-occurrence")
  write.csv(df_nocoop, "cooc_no_coop.csv", row.names = FALSE)
  p <- ggplot(df_tab, aes(y = joint, x = cooc)) +
    geom_segment(aes(x = 0, xend = cooc,
                     y = joint, yend = joint
                     #, linetype = faded
    ), colour = "grey50"
    ) +
    geom_vline(xintercept = 0, colour = "grey50", linetype = 1) + xlab("Co-occurrence frequency") +
    geom_point(data = df_tab, colour = "black", fill = "black", shape = 21, size = 1.5) +
    # scale_colour_manual(values = c(Outcome = "gray50", Indicator = "tomato", Cause = "gold", Protective = "forestgreen"), guide = NULL)+
    # scale_fill_manual(values = c(Outcome = "gray50", Indicator = "tomato", Cause = "gold", Protective = "forestgreen")) +
    scale_x_sqrt(expand = c(0,0)) +
    #scale_linetype_manual(values = c("TRUE" = 2, "FALSE" = 1), guide = NULL) +
    theme_bw() + theme(panel.grid.major.x = element_blank(),
                       panel.grid.minor.x = element_blank(), axis.title.y = element_blank(),
                       legend.position = c(.70,.125),
                       legend.title = element_blank(),
                       axis.text.y = element_text(hjust=0, vjust = 0, size = 6, family = "mono"))

  ggsave("plot_cooc_freq_nocoop.png", p, device = "png", width = 210, height = 400, units = "mm", dpi = 150)
  ggsave("plot_cooc_freq_nocoop.svg", p, device = "svg", width = 210, height = 400, units = "mm")

  # Create network ----------------------------------------------------------
  sq<-function(x){
    x^2
  }
  #inverse square function (square root)
  isq<-function(x){
    print(paste("isq",x))  #debug statement
    x<-ifelse(x<0, 0, x)
    sqrt(x)
  }
  edg <- df_plot
  edg$width = edg$cooc

  vert <- data.frame(name = as.character(unique(c(edg$term1, edg$term2))))
  vert$size <- cluster_freq$frequency[match(vert$name, cluster_freq$construct)]


  vert$size <- scales::rescale(vert$size^2, c(3, 25))

  g <- igraph::graph_from_data_frame(edg, vertices = vert,
                                     directed = FALSE)
  # node_centrality <- sort(igraph::degree(g, mode="all"))
  # write.csv(data.frame(construct = names(node_centrality), centrality_degree = node_centrality), "node_centrality.csv", row.names = FALSE)
  # edge thickness
  E(g)$width <- scales::rescale(sqrt(E(g)$width), to = c(.5, 8))

  edge.start <- ends(g, es=E(g), names = FALSE)[,1]
  edge.end <- ends(g, es=E(g), names = FALSE)[,2]

  # Layout
  set.seed(1) #4
  wts <- rep(1, length(E(g)$width))
  wts[edg$term1 == "cooperation" | edg$term2 == "cooperation"] <- 10
  lo <- layout_with_fr(g, weights = wts)
  vert <- data.frame(lo, vert[, c("name", "size")])
  vert$fill <- ordered(vert$size)
  edg$width <- scales::rescale(edg$width, to = c(.5, 2))
  df_from <- vert[edge.start, c("X1", "X2")]
  names(df_from) <- c("xstart", "ystart")
  df_to <- vert[edge.end, c("X1", "X2")]
  names(df_to) <- c("xend", "yend")
  edg <- data.frame(df_from, df_to, edg)

  p <- ggplot(data = NULL) +
    # Edges
    geom_segment(data = edg, aes(x = xstart, y = ystart, xend = xend, yend = yend, linewidth = width, alpha = width))+
    geom_point(data = vert, aes(x = X1, y = X2, size = size, fill = size), colour="black",pch=21, alpha = 1) +
    scale_fill_gradient(high = "#FF0000", low = "#FFFF00") +
    ggrepel::geom_label_repel(data = vert, aes(x = X1, y = X2, label = name), max.overlaps = 25) +

    theme_void() +
    scale_linewidth(range = c(.5, 3)) +
    scale_alpha(range = c(.2, 1))+
    scale_size(range = c(2,15), trans = scales::trans_new("sq", sq, isq))+
    theme(legend.position = "none")

  p
  #ggsave("network.eps", p, device = "eps", width = 210, height = 210, units = "mm")
  ggsave("network.svg", p, device = "svg", width = 210, height = 210, units = "mm")
  ggsave("network.png", p, device = "png", width = 7, height = 7, units = "in")
  saveRDS(p, "network.RData")
  set.seed(64)
  shifter <- function(x, n = 1) {
    if (n == 0) x else c(tail(x, -n), head(x, n))
  }
  lo <- layout_in_circle(g, order = shifter(V(g), -3))

  vert[, c(1:2)] <- lo
  df_from <- vert[edge.start, c("X1", "X2")]
  names(df_from) <- c("xstart", "ystart")
  df_to <- vert[edge.end, c("X1", "X2")]
  names(df_to) <- c("xend", "yend")
  edg[,c("xstart", "ystart", "xend", "yend")] <- data.frame(df_from, df_to)

  p <- ggplot(data = NULL) +
    # Edges
    geom_segment(data = edg, aes(x = xstart, y = ystart, xend = xend, yend = yend, linewidth = width))+
    geom_point(data = vert, aes(x = X1, y = X2, size = size, fill = size), colour="black",pch=21, alpha = 1) +
    scale_fill_gradient(high = "#FF0000", low = "#FFFF00") +
    ggrepel::geom_label_repel(data = vert, aes(x = X1, y = X2, label = name), max.overlaps = 25) +

    theme_void() +
    scale_linewidth(range = c(1, 3)) +
    theme(legend.position = "none")
  ggsave("network_circle.svg", p, device = "svg", width = 7, height = 7, units = "in")
  ggsave("network_circle.png", p, device = "png", width = 7, height = 7, units = "in")

  return("network")
}

d_to_r <- function(d){
  d / sqrt((d^2 + 4))
}

plot_meta <- function(){
  library(ggplot2)
  library(ggrepel)
  cf <- read.csv("cluster_freq.csv", stringsAsFactors = F)
  df_meta <- read.csv("meta-analytic-effects.csv", stringsAsFactors = F)
  df_meta$Variable <- trimws(tolower(df_meta$V1))
  all(df_meta$Variable %in% tolower(cf$Construct))
  df_meta$Variable[!df_meta$Variable %in% tolower(cf$Construct)]

  df_meta$Frequency <- as.integer(cf$Frequency[match(df_meta$Variable, cf$Construct)])
  df_meta$effectsize <- factor(substr(df_meta$V2, 1,1), levels = c("d", "r"))
  df_meta$Value <- abs(as.numeric(gsub("^.+?=(.+)\\).?", "\\1", df_meta$V2)))
  df_meta$Value[df_meta$effectsize == "d"] <- d_to_r(df_meta$Value[df_meta$effectsize == "d"])
  file.remove("results_meta_analysis.txt")
  levels(df_meta$effectsize) <- paste0(levels(df_meta$effectsize), ", ", sapply(levels(df_meta$effectsize), function(lev){
    x <- df_meta[df_meta$effectsize == lev, ]
    tmp <- cor.test(x$Value, log(x$Frequency), use = "pairwise.complete.obs")
    n_obs <- sum(complete.cases(x[c("Value", "Frequency")]))
    write(paste0("Correlation ", lev, ": r ", worcs::report(tmp$estimate), ", p ", worcs::report(tmp$p.value), ", n = ", n_obs), "results_meta_analysis.txt", append = TRUE)
    formatC(cor(x$Value, log(x$Frequency), use = "pairwise.complete.obs"), digits = 2, format = "f")
  }))


  df_plot <- df_meta[, c("Variable", "Value", "Frequency", "effectsize")]
  df_plot <- na.omit(df_plot)

  p <- ggplot(df_plot, aes(x = Value, y = Frequency, shape = effectsize, color = effectsize, label = Variable)) +
    #geom_point() +
    ggrepel::geom_text_repel(max.overlaps = 20, size = 4) +
    #geom_text(aes(label = Variable, size = 1)) +
    theme_bw()+
    theme(legend.position = c(.9,.9), legend.title = element_blank())+
    scale_y_log10()
  ggsave("meta_analysis_plot.svg", p, device = "svg", width = 200, height = 180, units = "mm")
  return("meta_analysis_plot.svg")
}
# df_meta$Variable[df_meta$Variable %in% tolower(cf$Construct)]c

make_tables <- function(){
tmp <- read.csv("complete_node_centrality.csv", stringsAsFactors = F)
tmp <- tmp[order(tmp$frequency, decreasing = TRUE), ]
cls <- 3
rws <- ceiling(nrow(tmp) / cls)
out <- tmp[1:rws, ]
for(i in 2:cls){
  out <- cbind(out, tmp[(((i-1)*rws)+1):(i*rws), ])
}
names(out) <- rep(c("Construct", "Freq.", "Cent."), cls)

write.csv(out, "Table_X.csv", row.names = F, na = "")
return("Table_X.csv")
}
