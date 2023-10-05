#' Average Clusters within Shift
#'
#' @param start_tile
#' @param clusters
#' @param tile_shift
#' @param tile_height
#'
#'
average_clusters_within_shift <- function(start_tile, clusters, tile_shift, tile_height){

    tile <- rbind(as.numeric(start_tile), as.numeric(start_tile)+tile_shift)

    clusters %>%
        dplyr::filter(tile[1,1]>=x & tile[2,1]<=x+tile_height &
                          tile[1,2]>=y & tile[2,2]<=y+tile_height) %>%
        dplyr::count(cluster, .drop=FALSE)

}
