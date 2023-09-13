#' Identify Cell in Tile (Embedded in find_niches)
#'
#' @param df_spe_1
#' @param df_final_1
#' @param tile_shift
#' @param cores
#' @param prob
#'
#' @return
#' @export
#'
#' @examples
identify_cell_in_tile <- function(df_spe_1, df_final_1, tile_shift, cores=10, prob=FALSE){

    df_tmp <- df_spe_1[, c("CenterX_local_px", "CenterY_local_px")]
    df_tmp <- as.list(as.data.frame(t(df_tmp)))
    names(df_tmp) <- df_spe_1$Barcode

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

        names(cluster_tmp)<- names(df_tmp)
        cluster_tmp <- unlist(cluster_tmp)

    }

    return(cluster_tmp)

}