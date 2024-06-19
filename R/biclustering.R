spectral_coclustering <- function(X, n_clusters, clust_seed=NULL, arpack=TRUE, kmeanspp=TRUE, ...) {
  #' Spectral co-clustering
  #'
  #' Runs spectral co-clustering (biclustering) on matrix X.
  #'
  #' @param X a numeric matrix or a data frame of numeric values. In the case of negative values the data is scaled.
  #' @param n_clusters positive integer >= 2. The amount of clusters to use.
  #' @param clust_seed NULL or a single value, optional seed number designation for reproducible results.
  #' @param arpack boolean specifying whether arpack-method is used for single value decomposition. By default,
  #' SVD is conducted using irlba::svdr() function which uses the robust arpack method (recommended). The alternative uses the base::svd() function which might be faster with large matrices.
  #' @param kmeanspp boolean specifying whether to use kmeans++ for initialization (recommended). By default the clustering is conducted using maotai::kmeanspp(), which produces more accurate results (recommended).
  #' If faster run times are desired, the alternative uses stats::kmeans().
  #' @param ... other arguments are passed on to stats::kmeans.
  #'
  #' @details
  #' Biclustering algorithms simultaneously cluster rows and columns of a data matrix.
  #' These clusters of rows and columns are known as biclusters, which each determine a submatrix of the original data matrix with some desired properties.
  #' For instance, given a matrix of shape (10, 10), one possible bicluster with three rows and two columns induces a submatrix of shape (3, 2).
  #'
  #' The spectral coclustering algorithm models a data matrix as a bipartite graph between columns and rows,
  #' where simultaneous clustering problem can be posed as a bipartite graph partitioning problem.
  #' This is solved using an appropriate scaling of the data matrix followed by singular value decomposition.
  #' The resulting graph yields a checkerboard structure, where each row and columns belong exactly one bicluster.
  #'
  #' @references
  #' https://doi.org/10.1145/502512.502550
  #'
  #' @returns Returns a list with entries:
  #' row_vec: a vector of cluster assignments for the rows
  #' col_vec: a vector of cluster assignments for the columns
  #' n_cluster: integer, the number of clusters
  #'
  #' @examples
  #' {#generate data and visualize it
  #' data <- generate_example_matrix()
  #' heatmap(data$original, Rowv = NA, Colv = NA)
  #'
  #' #run coclustering on the randomized data and visualize the reordered data
  #' bc <- spectral_coclustering(data$randomized, n_clusters = 4, arpack = TRUE)
  #' r.ind <- sort.int(bc$row_vec,index.return = TRUE)$ix
  #' c.ind <- sort.int(bc$col_vec,index.return = TRUE)$ix
  #' heatmap(
  #'   data$randomized[r.ind,c.ind], Rowv = NA, Colv = NA)
  #' }
  #'
  #' @importFrom stats kmeans
  #' @export
  cross <- as.matrix(X)
  if (any(!is.numeric(cross))) {
    print("Matrix is not numeric")
    return()
  }
  if (any(is.na(cross))) {
    print("Matrix contains missing values")
    return()
  }
  if (any(cross < 0)) cross <- cross + abs(min(cross))
  R <- diag(x=1/sqrt(rowSums(cross)))
  C <- diag(x=1/sqrt(colSums(cross)))
  An <- R %*% cross %*% C
  n_svd <- 1+ceiling(log2(n_clusters))
  if(is.numeric(clust_seed)) { set.seed(clust_seed) }

  if (arpack == TRUE) { SVD <- irlba::svdr(An, n_svd) }
  else { SVD <- svd(An,nu=n_svd,nv=n_svd) }

  Z <- rbind(R %*% SVD$u[,-1], C %*% SVD$v[,-1])
  if(kmeanspp) {
    clusters <- list("cluster"=maotai::kmeanspp(Z, k = n_clusters)) }
  else {
    clusters <- kmeans(Z, centers = n_clusters, ...)
  }
  row_vec <- clusters$cluster[1:dim(cross)[1]]
  names(row_vec) <- 1:length(row_vec)
  col_vec <- clusters$cluster[dim(cross)[1]+1:dim(cross)[2]]
  names(col_vec) <- 1:length(col_vec)
  results <- list(
    "row_vec" = row_vec,
    "col_vec" = col_vec,
    n_clusters = n_clusters
  )
  return(results)
}


plot_scores <- function(X, min_clusters=2, max_clusters=10, score_a="ave.within.cluster.ss", score_b="ch", score_seed=NULL, ...) {
  #' Plot clustering scores
  #'
  #' Runs the spectral co-clustering on matrix X using the selected cluster amounts, runs cluster scoring, and plots clustering scores for both dimensions.
  #'
  #' @param X a numeric matrix or a data frame of numeric values.
  #' @param min_clusters positive integer >= 2. The minimum amount of clusters to try.
  #' @param max_clusters positive integer. The maximum amount of clusters to try.
  #' @param score_a character string specifying the first score to be plotted. Defaults to "ave.within.cluster.ss" (average within-cluster sum of squares). Acceptable inputs are e.g. "avg.silwidth", "within.cluster.ss", "ch", and the other scores returned by fpc::cluster.stats.
  #' @param score_b character string specifying the second score to be plotted. Defaults to "ch" (calinski-harabasz score).
  #' @param score_seed NULL or a single value, optional seed number for reproducible results
  #' @param ... other arguments are passed on to jhbiclusters::spectral_coclustering.
  #'
  #' @details
  #' Determining the optimal number of biclusters in absence of *a priori* knowledge can be done by assessing each dimension (rows, columns)
  #' separately using convenient metrics. This function uses the fpc::cluster_stats() function to run scoring across cluster numbers min_clusters:max_clusters for both dimensions.
  #'
  #' @returns a list with length 2:
  #' row_scores: clustering scores along the rows, for cluster numbers min_clusters:max_clusters
  #' col_scores: clustering scores along the columns
  #'
  #' @examples
  #' {# generate data and visualize the clustering scores
  #' data <- generate_example_matrix()
  #' scores <- plot_scores(X = data$randomized, score_a = "ch", score_b = "avg.silwidth", score_seed = 0)
  #' }
  #'
  #' @import ggplot2
  #' @importFrom stats dist
  #' @export
  if ((!requireNamespace("fpc", quietly = TRUE)) | !requireNamespace("cluster", quietly = TRUE)) {
    stop(
      "Packages \"fpc\" and \"cluster\" must be installed to use this function.",
      call. = FALSE
    ) }
  Xr <- dist(X, "euclidean")
  Xc <- dist(t(X), "euclidean")

  scores <- c()
  scores_c <- c()

  if (is.numeric(score_seed)) {set.seed(score_seed)}
  for (i in (min_clusters:max_clusters)) {
    bc <- spectral_coclustering(X, n_clusters = i, ...)
    suppressWarnings({
    scores <- rbind(scores, t(fpc::cluster.stats(Xr, aggregateonly=TRUE, bc$row_vec)))
    scores_c <- rbind(scores_c, t(fpc::cluster.stats(Xc, aggregateonly=TRUE, bc$col_vec)))
    })
  }
  if ((!requireNamespace("ggplot2", quietly = TRUE)) | !requireNamespace("patchwork", quietly = TRUE)) {

    stop("Packages \"ggplot2\" and \"patchwork\" must be installed to make visualizations.")
    return(list(row_scores=scores, col_scores=scores_c))
  }
  p <- patchwork::wrap_plots(
    ggplot(data=data.frame(row=as.numeric(scores[,score_a]),col=as.numeric(scores_c[,score_a]),
                      cn=as.numeric(scores[,"cluster.number"])), aes(x=cn)) +
      geom_path(aes(y=row), color="green", linewidth = 2, alpha=0.6) +
      geom_path(aes(y=col), color="red", linewidth = 2, alpha=0.6) +
      ylab(paste0(score_a," (green=rows, red=columns)"))+xlab("n clusters") + theme_bw(),

    ggplot(data=data.frame(row=as.numeric(scores[,score_b]),col=as.numeric(scores_c[,score_b]),
                      cn=as.numeric(scores[,"cluster.number"])), aes(x=cn)) +
      geom_path(aes(y=row), color="darkgreen", linewidth = 2, alpha=0.6) +
      geom_path(aes(y=col), color="darkred", linewidth = 2, alpha=0.6) +
      ylab(paste0(score_b," (green=rows, red=columns)"))+xlab("n clusters") + theme_bw()
  )
  print(p)
  return(list(row_scores=scores, col_scores=scores_c))
}


generate_example_matrix <- function() {
  #' Generate an example matrix
  #'
  #' Generates a numeric matrix, simulating correlation coefficients in a dataset, and randomizes it.
  #'
  #' @returns A list with length 2:
  #' original:  A numeric 14 x 14 matrix with a checkerboard structure.
  #' randomized:  The randomized matrix.
  #'
  #' @importFrom stats runif
  #' @export
  m1 <- data.frame(matrix(runif(16, 0.5, 0.9), nrow = 4, dimnames = list(paste0("R0",0:3),paste0("C0",0:3))))
  m2 <- data.frame(matrix(runif(9, 0.4, 1), nrow=3, dimnames = list(paste0("R0",4:6),paste0("C0",4:6))))
  m3 <- data.frame(matrix(runif(9, 0.4, 1), nrow=3, dimnames = list(paste0("R0",7:9),paste0("C0",7:9))))
  m4 <- data.frame(matrix(runif(18, 0.5, 0.7), nrow = 3, dimnames = list(paste0("R",11:13),paste0("C",11:16))))

  mm <- merge(merge(merge(m1, m2, by.y=0, by.x=0, all=TRUE), m3, by.x="Row.names", by.y=0, all=TRUE), m4, by.x="Row.names", by.y=0, all=TRUE)
  rownames(mm) <- mm[,1]
  mm <- mm[,-1]
  mm[is.na(mm)] <- runif(length(mm[is.na(mm)]),-0.5,0.15)
  return(list("original"=as.matrix(mm), "randomized"=as.matrix(mm[sample(rownames(mm)),sample(colnames(mm))])))
}
