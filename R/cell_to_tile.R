#' Map Tiles Back to Cells from find_niches Output
#'
#' @param df_final A dataframe object from find_niches output.
#' @param spe A SpatialExperiment object.
#' @param fov A character string specifying the name of the fields of view column in spe (Default: 'fov').
#' @param barcode A character string specifying the name of the column containing cell barcode information in spe (Default: 'Barcode'). Barcode structure: 'CellID_fov'.
#' @param data_type A character string indicating if spe contains polygon or dot/point information (Options: 'polygon' or 'dot').
#' @param x_coord A character string specifying the name of the column containing x spatial co-ordinates (Default: 'x_local_px').
#' @param y_coord A character string specifying the name of the column containing y spatial co-ordinates (Default: 'y_local_px').
#' @param x_centre A character string specifying the name of the column containing x centroid co-ordinates (Default: 'CenterX_local_px'). Applicable for data_type='polygon'.
#' @param y_centre A character string specifying the name of the column containing y centroid co-ordinates (Default: 'CenterY_local_px'). Applicable for data_type='polygon'.
#' @param tile_shift A vector specifying tile shift.
#' @param cores A vector specifying the number of cores to implement (Default: 10).
#' @param prob A logical value indicating if the user in interested in calculating the probability (Default: FALSE).
#'
#' @return A list of cells mapped to tiles
#' @export
#'
#' @examples
cell_to_tile <- function(df_final, spe, fov="fov", barcode="Barcode", data_type=c("polygon","dot"), x_coord="x_local_px", y_coord="y_local_px",
                         x_centre="CenterX_local_px", y_centre="CenterY_local_px", tile_shift, cores=10, prob=FALSE){

    data_type <- match.arg(data_type)

    if (data_type == "polygon") {
        x_use <- x_centre
        y_use <- y_centre
    } else { # dot
        x_use <- x_coord
        y_use <- y_coord
    }

    index <- split(1:nrow(df_final), df_final$fov)
    df_final_split <- lapply(index, function(x) df_final[x,])

    index <- split(1:ncol(spe), spe$fov)
    spe_split <- lapply(index, function(x) as.data.frame(colData(spe)[x,
                                                                      c(x_use, y_use, barcode)]))
    df_final_split <- df_final_split[names(spe_split)]

    all_cluster <- mapply(Y = df_final_split, X = spe_split,
                          function(X,Y) identify_cell_in_tile(X, Y, tile_shift, cores=10, prob=prob, x_use, y_use, barcode))


    return(all_cluster)

}
