#' Shuffle Clusters to Create Null Distribution
#'
#'
#'
shuffle_clusters <- function(clusters_df) {

    # choose a random reference point
    ref_point <- clusters_df[sample(nrow(clusters_df), 1), ]

    # calculate distances from the reference point
    distances <- sqrt((clusters_df$x - ref_point$x)^2 + (clusters_df$y - ref_point$y)^2)

    # order based on distance to reference point
    order_indices <- order(distances)

    # shuffle cluster labels of the ordered data
    shuffled_labels <- sample(clusters_df$cluster[order_indices])

    shuffled_clusters_df <- clusters_df
    shuffled_clusters_df$cluster[order_indices] <- shuffled_labels

    return(shuffled_clusters_df)

}
