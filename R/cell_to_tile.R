#' Map Tiles Back to Cells from find_niches Output
#'
#' @param df_final A dataframe object from the find_niches output.
#' @param spe A SpatialExperiment object.
#' @param tile_shiftA vector specifying tile shift.
#' @param cores A vector specifying the number of cores to implement (Default: 10).
#' @param prob A logical value indicating if the user in interested in calculating the probability (Default: FALSE).
#'
#' @return A list of cells mapped to tiles
#' @export
#'
#' @examples
cell_to_tile <- function(df_final, spe,  tile_shift, cores=10, prob=FALSE){

    index <- split(1:nrow(df_final), df_final$fov)
    df_final_split <- lapply(index, function(x) df_final[x,])

    index <- split(1:ncol(spe), spe$fov)
    spe_split <- lapply(index, function(x) as.data.frame(colData(spe)[x,
                                                                      c("CenterX_local_px", "CenterY_local_px", "Barcode")]))
    df_final_split <- df_final_split[names(spe_split)]

    all_cluster <- mapply(Y = df_final_split, X = spe_split,
                          function(X,Y) identify_cell_in_tile(X, Y, tile_shift, cores=10, prob=prob))


    return(all_cluster)

}
