#' Shuffle Clusters to Create Null Distribution
#'
#'
#'
shuffle_clusters <- function(clusters_df) {

    ref_point <- clusters_df[sample(nrow(clusters_df), 1), ]
    distances <- sqrt((clusters_df$x - ref_point$x)^2 + (clusters_df$y - ref_point$y)^2)
    order_indices <- order(distances)

    shuffled_labels <- sample(clusters_df$Cluster[order_indices])
    shuffled_clusters_df <- clusters_df
    shuffled_clusters_df$cluster[order_indices] <- shuffled_labels

    return(shuffled_clusters_df)

}
