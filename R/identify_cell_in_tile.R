#' Identify Cell in Tile
#'
#' @param df_spe_1
#' @param df_final_1
#' @param tile_shift
#' @param x_use
#' @param y_use
#' @param barcode
#' @param cores
#' @param prob
#'
#' @return
#' @export
#'
#' @examples
identify_cell_in_tile <- function(df_spe_1, df_final_1, tile_shift,
                                  x_use, y_use, barcode,
                                  cores=10, prob=FALSE) {

    df_tmp <- df_spe_1[, c(x_use, y_use)]
    df_tmp <- as.list(as.data.frame(t(df_tmp)))
    names(df_tmp) <- df_spe_1[, barcode]

    if(prob){
        cluster_tmp <- mclapply(df_tmp, function(x) {
            tmp <- find_tile(x, df_final_1, tile_shift)
            df_final_1[tmp, -c(1:4)]
        }, mc.cores=cores)
    } else {
        cluster_tmp <- mclapply(df_tmp, function(x) {
            tmp <- find_tile(x, df_final_1, tile_shift)
            df_final_1$Cluster[tmp]
        }, mc.cores=cores)

        names(cluster_tmp) <- names(df_tmp)
        cluster_tmp <- unlist(cluster_tmp)
    }

    return(cluster_tmp)
}
