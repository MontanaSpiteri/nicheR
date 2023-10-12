#' Compute Test Statistic
#'
#'
#'
compute_test_statistic <- function(clusters_df, dist_tiles) {

    grouped <- split(1:nrow(dist_tiles), clusters_df$cluster)
    intra_cluster_distances <- sapply(grouped, function(indices) {
        mean(dist_tiles[indices, indices])

    })

    return(mean(intra_cluster_distances))
}
NULL

#' Compute P Value Statistics
#'
#'
#'
compute_pvalue <- function(clusters_df, dist_tiles, n_permutations) {

    original_statistic <- compute_test_statistic(clusters_df, dist_tiles)

    null_statistics <- replicate(n_permutations, {
        shuffled_clusters_df <- shuffle_clusters(clusters_df)
        compute_test_statistic(shuffled_clusters_df, dist_tiles)
    })

    p_value <- mean(null_statistics >= original_statistic)
    #compute p-value (proportion of times the null statistics exceed the original statistic)

    return(p_value)

}
